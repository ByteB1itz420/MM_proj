"""Two-node lumped thermal model of an air-cooled data hall.

Implements PROJECT_SPEC.md §4-§5 (v0.2). Stage 1 only: physics + verification.
No controllers, optimisation or energy-saving calculations live here.

Units (SI throughout):
    temperatures      degC (differences in K); absolute K only inside cop()
    heat flows        W_th   (watts of heat)
    electrical power  W_e    (watts of electricity)
    energy            J internally; convert with J_to_kWh() for reporting
    time              s

State      x = [T_s, T_a]                  (degC)
Inputs     w = [Q_c, P_IT, T_out]          (W_th, W_e, degC)
Augmented  X = [T_s, T_a, E_Q, E_IT, E_env]
    E_Q   = integral of Q_c dt                     heat removed by cooling   (J_th)
    E_IT  = integral of P_IT dt                    IT heat released          (J)
    E_env = integral of UA_e (T_out - T_a) dt      net heat in via envelope  (J_th)
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

J_PER_KWH = 3.6e6  # 1 kWh = 1000 W * 3600 s
KELVIN_OFFSET = 273.15

N_X = 2  # thermal states
N_Z = 3  # energy accumulators
N_W = 3  # inputs


def J_to_kWh(energy_J):
    return np.asarray(energy_J) / J_PER_KWH


def kWh_to_J(energy_kWh):
    return np.asarray(energy_kWh) * J_PER_KWH


@dataclass(frozen=True)
class Params:
    """Model parameters. [STD] = standard constant, [ILL] = illustrative (not measured)."""

    rho: float = 1.2          # air density, kg/m^3                       [STD]
    c_p: float = 1005.0       # air specific heat, J/(kg K)               [STD]
    V: float = 3000.0         # hall air volume, m^3                      [ILL]
    C_s: float = 2.5e7        # lumped IT-equipment heat capacity, J/K    [ILL]
    UA_s: float = 1e6 / 15    # IT -> air conductance, W/K                [ILL] (15 K rise at 1 MW)
    UA_e: float = 2e3         # envelope conductance, W/K                 [ILL]
    Q_max: float = 1.5e6      # cooling capacity, W_th                    [ILL]
    eta_II: float = 0.45      # fraction of Carnot COP achieved           [ILL]
    dT_ev: float = 12.0       # room air -> evaporator approach, K        [ILL]
    dT_cd: float = 8.0        # condenser -> outdoor approach, K          [ILL]
    COP_cap: float = 10.0     # upper edge of the COP model's domain      [ILL]

    @property
    def C_a(self):
        """Air heat capacity rho * c_p * V, J/K."""
        return self.rho * self.c_p * self.V


# ---------------------------------------------------------------------------
# Equations (1)-(3)
# ---------------------------------------------------------------------------

def rhs(p, x, w):
    """Right-hand side of eqs. (1)-(2) written term by term, degC/s."""
    T_s, T_a = x
    Q_c, P_IT, T_out = w
    dT_s = (P_IT - p.UA_s * (T_s - T_a)) / p.C_s
    dT_a = (p.UA_s * (T_s - T_a) + p.UA_e * (T_out - T_a) - Q_c) / p.C_a
    return np.array([dT_s, dT_a])


def system_matrices(p):
    """Matrix form x' = A x + B u + E d of eq. (3); returns A, B, E and G = [B E]."""
    A = np.array([
        [-p.UA_s / p.C_s, p.UA_s / p.C_s],
        [p.UA_s / p.C_a, -(p.UA_s + p.UA_e) / p.C_a],
    ])
    B = np.array([[0.0], [-1.0 / p.C_a]])
    E = np.array([
        [1.0 / p.C_s, 0.0],
        [0.0, p.UA_e / p.C_a],
    ])
    return A, B, E, np.hstack([B, E])


def augmented_matrices(p):
    """X' = F X + H w for the state plus the three energy accumulators."""
    A, _, _, G = system_matrices(p)
    K = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, -p.UA_e]])   # E_env' depends on T_a
    L = np.diag([1.0, 1.0, p.UA_e])                          # E_Q' = Q_c, E_IT' = P_IT, E_env' += UA_e T_out
    F = np.block([[A, np.zeros((N_X, N_Z))], [K, np.zeros((N_Z, N_Z))]])
    H = np.vstack([G, L])
    return F, H


def stored_energy(p, x):
    """U = C_s T_s + C_a T_a, J relative to 0 degC. Works on (..., 2) arrays."""
    x = np.asarray(x)
    return p.C_s * x[..., 0] + p.C_a * x[..., 1]


# ---------------------------------------------------------------------------
# Exact zero-order-hold propagation
# ---------------------------------------------------------------------------

def zoh_step_matrices(p, h):
    """Phi, Gamma with X(t+h) = Phi X(t) + Gamma w, exact when w is constant on [t, t+h).

    Uses the augmented matrix exponential expm([[F, H], [0, 0]] h) = [[Phi, Gamma], [0, I]].
    """
    F, H = augmented_matrices(p)
    n, m = F.shape[0], H.shape[1]
    M = np.zeros((n + m, n + m))
    M[:n, :n] = F
    M[:n, n:] = H
    Md = expm(M * h)
    return Md[:n, :n], Md[:n, n:]


def _as_float_array(name, value):
    try:
        return np.asarray(value, dtype=float)
    except (TypeError, ValueError) as err:
        raise ValueError(f"{name} must be numeric: {err}") from None


def validate_simulation_inputs(x0, t_grid, W):
    """Check and convert simulator inputs; raise ValueError with a clear message if invalid.

    x0: shape (2,) finite. t_grid: 1-D, >= 2 finite points, strictly increasing.
    W: shape (len(t_grid) - 1, 3) finite (a single row (3,) is accepted when there is one interval).
    """
    x0 = _as_float_array("x0", x0)
    t_grid = _as_float_array("t_grid", t_grid)
    W = _as_float_array("W", W)
    if x0.shape != (N_X,):
        raise ValueError(f"x0 must have shape ({N_X},) = [T_s, T_a]; got {x0.shape}")
    if t_grid.ndim != 1 or t_grid.size < 2:
        raise ValueError(f"t_grid must be 1-D with at least 2 points; got shape {t_grid.shape}")
    if W.ndim == 1 and t_grid.size == 2:
        W = W[np.newaxis, :]
    if W.shape != (t_grid.size - 1, N_W):
        raise ValueError(f"W must have shape (len(t_grid) - 1, {N_W}) = ({t_grid.size - 1}, {N_W}) "
                         f"with columns [Q_c, P_IT, T_out]; got {W.shape}")
    for name, arr in (("x0", x0), ("t_grid", t_grid), ("W", W)):
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} contains non-finite values (NaN or inf)")
    steps = np.diff(t_grid)
    if np.any(steps <= 0):
        k = int(np.argmax(steps <= 0))
        raise ValueError(f"t_grid must be strictly increasing; t[{k + 1}] = {t_grid[k + 1]} <= t[{k}] = {t_grid[k]}")
    return x0, t_grid, W


def simulate_zoh(p, x0, t_grid, W):
    """Propagate exactly across breakpoints t_grid with W[k] held on [t_k, t_{k+1}).

    Returns X with shape (len(t_grid), 5); accumulators start at zero.
    Raises ValueError for invalid shapes, non-finite values or a non-increasing time grid.
    """
    x0, t_grid, W = validate_simulation_inputs(x0, t_grid, W)
    X = np.zeros((len(t_grid), N_X + N_Z))
    X[0, :N_X] = x0
    cache = {}
    for k, h in enumerate(np.diff(t_grid)):
        key = round(h, 9)
        if key not in cache:
            cache[key] = zoh_step_matrices(p, h)
        Phi, Gam = cache[key]
        X[k + 1] = Phi @ X[k] + Gam @ W[k]
    return X


def simulate_ivp(p, x0, t_grid, W, method="Radau", rtol=1e-10, atol=None, eval_step=None):
    """Independent reference: solve_ivp on the augmented ODE, restarted at every input change.

    Splitting at breakpoints means the integrator never steps across a discontinuity.
    If eval_step (s) is given, also returns (t_dense, X_dense) sampled every eval_step inside
    each interval, for plotting or point-by-point comparison.
    Raises ValueError for invalid shapes, non-finite values or a non-increasing time grid.
    """
    x0, t_grid, W = validate_simulation_inputs(x0, t_grid, W)
    if eval_step is not None and not (np.isfinite(eval_step) and eval_step > 0):
        raise ValueError(f"eval_step must be a positive finite number of seconds; got {eval_step}")
    F, H = augmented_matrices(p)
    if atol is None:
        atol = np.array([1e-10, 1e-10, 1e-2, 1e-2, 1e-2])
    X = np.zeros((len(t_grid), N_X + N_Z))
    X[0, :N_X] = x0
    dense_t, dense_X = [t_grid[:1]], [X[:1]]
    for k in range(len(t_grid) - 1):
        bw = H @ W[k]
        t_eval = None
        if eval_step:
            t_eval = np.arange(t_grid[k] + eval_step, t_grid[k + 1] - 1e-9 * eval_step, eval_step)
            t_eval = np.append(t_eval, t_grid[k + 1])
        sol = solve_ivp(lambda t, X_, bw=bw: F @ X_ + bw, (t_grid[k], t_grid[k + 1]), X[k],
                        method=method, rtol=rtol, atol=atol, jac=F, t_eval=t_eval)
        if not sol.success:
            raise RuntimeError(f"solve_ivp failed on interval {k}: {sol.message}")
        X[k + 1] = sol.y[:, -1]
        if eval_step:
            dense_t.append(sol.t)
            dense_X.append(sol.y.T)
    if eval_step:
        return X, np.concatenate(dense_t), np.vstack(dense_X)
    return X


def energy_balance_residual(p, X):
    """U(t_f) - U(0) - (E_IT + E_env - E_Q), J. Zero if eq. (4) is integrated consistently."""
    dU = stored_energy(p, X[-1, :N_X]) - stored_energy(p, X[0, :N_X])
    E_Q, E_IT, E_env = X[-1, N_X:] - X[0, N_X:]
    return dU - (E_IT + E_env - E_Q)


# ---------------------------------------------------------------------------
# Equilibria and modes
# ---------------------------------------------------------------------------

def equilibrium(p, Q_c, P_IT, T_out):
    """Analytic steady state (spec §8.1): T_a* = T_out + (P - Q)/UA_e, T_s* = T_a* + P/UA_s."""
    T_a = T_out + (P_IT - Q_c) / p.UA_e
    return np.array([T_a + P_IT / p.UA_s, T_a])


def holding_cooling(p, T_a_target, P_IT, T_out):
    """Constant Q_c (W_th) whose equilibrium has T_a = T_a_target."""
    return P_IT + p.UA_e * (T_out - T_a_target)


def eigenvalues(p):
    """Eigenvalues of A (1/s), ascending, via the symmetric similar matrix D A D^-1.

    D = diag(sqrt(C_s), sqrt(C_a)) makes A symmetric, which proves the eigenvalues are real.
    """
    off = p.UA_s / np.sqrt(p.C_s * p.C_a)
    S = np.array([[-p.UA_s / p.C_s, off], [off, -(p.UA_s + p.UA_e) / p.C_a]])
    return np.linalg.eigvalsh(S)


def time_constants(p):
    """(tau_fast, tau_slow) in s, tau = -1/lambda."""
    lam = eigenvalues(p)
    return -1.0 / lam[0], -1.0 / lam[1]


# ---------------------------------------------------------------------------
# COP model, eq. (6), with an explicit validity domain
# ---------------------------------------------------------------------------

class COPDomainError(ValueError):
    pass


def cop(p, T_a, T_out):
    """COP = eta_II * T_ev / (T_cd - T_ev), temperatures in K.

    Valid domain (each violation raises COPDomainError):
      1. T_a and T_out finite;
      2. evaporator absolute temperature T_ev > 0 K;
      3. positive lift T_cd - T_ev > 0;
      4. COP <= COP_cap, i.e. lift >= eta_II T_ev / COP_cap.
    Outside 3-4 (small or negative lift, i.e. cool outdoor air) the Carnot-fraction model does not
    describe the plant (free cooling is not modelled), so we raise instead of silently capping.
    Check 2 matters because with T_ev <= 0 K the cap test alone would pass and return COP <= 0.
    """
    try:
        T_a = np.asarray(T_a, dtype=float)
        T_out = np.asarray(T_out, dtype=float)
    except (TypeError, ValueError) as err:
        raise COPDomainError(f"temperatures must be numeric: {err}") from None
    if not (np.all(np.isfinite(T_a)) and np.all(np.isfinite(T_out))):
        raise COPDomainError("non-finite temperature (NaN or inf) passed to cop()")
    T_ev = T_a - p.dT_ev + KELVIN_OFFSET
    T_cd = T_out + p.dT_cd + KELVIN_OFFSET
    if np.any(T_ev <= 0):
        raise COPDomainError(
            f"evaporator absolute temperature must be > 0 K; got T_ev = {np.min(T_ev):.2f} K "
            f"(T_a = {np.min(T_a):.2f} degC)")
    lift = T_cd - T_ev
    if np.any(lift <= 0):
        raise COPDomainError(
            f"temperature lift T_cd - T_ev must be > 0 K; got {np.min(lift):.2f} K")
    lift_min = p.eta_II * T_ev / p.COP_cap
    if np.any(lift < lift_min):
        raise COPDomainError(
            f"COP model outside its domain: lift {np.min(lift):.3f} K < required {np.max(lift_min):.3f} K")
    return p.eta_II * T_ev / lift


def min_valid_T_out(p, T_a):
    """Lowest outdoor temperature (degC) for which cop() is inside its domain at room temperature T_a."""
    T_ev = np.asarray(T_a, dtype=float) - p.dT_ev + KELVIN_OFFSET
    return T_ev + p.eta_II * T_ev / p.COP_cap - p.dT_cd - KELVIN_OFFSET
