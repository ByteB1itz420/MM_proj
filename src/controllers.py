"""PI cooling controllers and closed-loop simulator (PROJECT_SPEC.md §7.1-7.2).

B1 FixedSetpointPI: constant setpoint. B2 LoadAdaptivePI: same PI law, setpoint recomputed from the
measured IT load at every update. No optimisation or MPC lives here.

Sign convention: error e = T_a - T_set. Positive error (air too warm) increases cooling Q_c.
Units: Kp in W_th/K, Ki in W_th/(K s), integrator state I in W_th, dt in s.

PI update order at each control time t_k = k*dt (implemented in pi_update):
    0. setpoint T_set,k = ctrl.setpoint(P_IT(t_k)) (constant for B1)
    1. measure T_a(t_k) from the exact plant state
    2. e_k = T_a(t_k) - T_set,k
    3. tentative integrator  I_tent = I_{k-1} + Ki * dt * e_k
    4. tentative command     v = Q_ff + Kp * e_k + I_tent
    5. conditional integration (anti-windup): if the update would push the command past a limit,
       i.e. (v > Q_max and e_k > 0) or (v < Q_min and e_k < 0), integrate only up to that limit:
           high: I_k = max(I_{k-1}, Q_max - Q_ff - Kp e_k)
           low:  I_k = min(I_{k-1}, Q_min - Q_ff - Kp e_k)
       so the integrator never moves past the boundary and never moves backwards. If the command was
       already beyond the limit, this leaves I_k = I_{k-1} (fully blocked). Otherwise I_k = I_tent;
       updates that move back out of saturation are always accepted.
       Invariant: whenever integration is limited, the applied Q_c is at that limit. (An earlier draft
       rejected the whole update instead; that could freeze the integrator just inside a limit and
       keep cooling on while the room was already too cold, so it was replaced.)
    6. Q_c,k = clip(Q_ff + Kp * e_k + I_k, Q_min, Q_max)
    7. hold Q_c,k constant on [t_k, t_k + dt) while the plant is propagated exactly (ZOH),
       splitting the interval at every disturbance change.
"""

from dataclasses import dataclass, field

import numpy as np

import thermal_model as tm


# ---------------------------------------------------------------------------
# Setpoint and gains
# ---------------------------------------------------------------------------

def setpoint_from_spec(p, P_IT_max, T_a_max=27.0, T_s_max=40.0, margin=1.0):
    """Spec §7.1: T_set = min(T_a,max, T_s,max - P_IT,max / UA_s) - m.

    T_s_max (40 degC) and the margin m are illustrative; m = 1 K is provisional (spec D3).
    """
    return min(T_a_max, T_s_max - P_IT_max / p.UA_s) - margin


@dataclass(frozen=True)
class PIGains:
    Kp: float  # W_th per K of error
    Ki: float  # W_th per (K s)
    dt: float  # control update interval, s

    def __post_init__(self):
        for name in ("Kp", "Ki", "dt"):
            v = getattr(self, name)
            if not np.isfinite(v):
                raise ValueError(f"{name} must be finite; got {v}")
        if self.Kp < 0 or self.Ki < 0:
            raise ValueError(f"gains must be non-negative; got Kp={self.Kp}, Ki={self.Ki}")
        if self.dt <= 0:
            raise ValueError(f"dt must be positive; got {self.dt}")

    @property
    def Ti(self):
        """Integral time Kp / Ki, s."""
        return self.Kp / self.Ki if self.Ki > 0 else np.inf


# Gains selected in notebook 02 by design_gains(p, dt, STAGE2_KP_GRID, STAGE2_KI_GRID, gain_margin=2.0)
# with the default Params. Stored exactly so later stages reuse them without retuning.
STAGE2_KP_GRID = np.concatenate([[0.0], np.geomspace(1e3, 3e5, 31)])
STAGE2_KI_GRID = np.geomspace(1.0, 1e4, 41)
STAGE2_GAINS = {
    60: PIGains(Kp=8096.1359720760265, Ki=630.9573444801937, dt=60.0),
    300: PIGains(Kp=0.0, Ki=199.52623149688807, dt=300.0),
}


def closed_loop_matrix(p, gains):
    """Unsaturated discrete closed-loop matrix M for z_k = [T_s(t_k), T_a(t_k), I_{k-1}].

    From the update order: u_k = Q_ff + (Kp + Ki dt)(T_a,k - T_set) + I_{k-1} and
    I_k = I_{k-1} + Ki dt (T_a,k - T_set), with x_{k+1} = Phi x_k + Gamma_u u_k + Gamma_d d_k.
    """
    Phi, Gam = tm.zoh_step_matrices(p, gains.dt)
    Phi, Gu = Phi[:2, :2], Gam[:2, :1]
    c = np.array([[0.0, 1.0]])
    kd = gains.Kp + gains.Ki * gains.dt
    return np.block([[Phi + kd * Gu @ c, Gu], [gains.Ki * gains.dt * c, np.ones((1, 1))]])


def closed_loop_poles(p, gains):
    return np.linalg.eigvals(closed_loop_matrix(p, gains))


def spectral_radius(p, gains):
    return float(np.max(np.abs(closed_loop_poles(p, gains))))


def design_gains(p, dt, Kp_grid, Ki_grid, gain_margin=2.0):
    """Grid search: minimise the unsaturated spectral radius subject to the loop staying
    stable when both gains are multiplied by gain_margin.

    Returns (best PIGains, info) with info = {"rho", "rho_margin", "Kp_edge", "Ki_edge", "n_feasible"};
    *_edge is "low"/"high" if the optimum lies on that grid boundary (the grid may be too narrow,
    unless the boundary is a physical limit such as Kp = 0), else None.
    """
    Kp_grid, Ki_grid = np.asarray(Kp_grid, dtype=float), np.asarray(Ki_grid, dtype=float)
    best, n_ok = None, 0
    for i, Kp in enumerate(Kp_grid):
        for j, Ki in enumerate(Ki_grid):
            rho = spectral_radius(p, PIGains(Kp, Ki, dt))
            rho_m = spectral_radius(p, PIGains(gain_margin * Kp, gain_margin * Ki, dt))
            if rho < 1 and rho_m < 1:
                n_ok += 1
                if best is None or rho < best[0]:
                    best = (rho, rho_m, i, j)
    if best is None:
        raise ValueError("no gain pair on the grid satisfies the stability and margin requirement")
    rho, rho_m, i, j = best
    edge = lambda idx, n: "low" if idx == 0 else ("high" if idx == n - 1 else None)
    return PIGains(Kp_grid[i], Ki_grid[j], dt), {"rho": rho, "rho_margin": rho_m, "n_feasible": n_ok,
                                                 "Kp_edge": edge(i, len(Kp_grid)), "Ki_edge": edge(j, len(Ki_grid))}


# ---------------------------------------------------------------------------
# PI law with conditional-integration anti-windup
# ---------------------------------------------------------------------------

SAT_TOL = 1e-9  # relative tolerance for "command at a limit" (clamped commands land on it up to round-off)


def pi_update(I_prev, e, gains, Q_ff, Q_min, Q_max, anti_windup=True):
    """One PI update (steps 3-6 of the module docstring).

    Returns (Q_c, I_new, info) where info holds the tentative values and flags for logging.
    """
    dI = gains.Ki * gains.dt * e
    I_tent = I_prev + dI
    v = Q_ff + gains.Kp * e + I_tent
    I_new = I_tent
    if anti_windup and v > Q_max and dI > 0:
        I_new = max(I_prev, Q_max - Q_ff - gains.Kp * e)
    elif anti_windup and v < Q_min and dI < 0:
        I_new = min(I_prev, Q_min - Q_ff - gains.Kp * e)
    limited = I_new != I_tent
    u_unsat = Q_ff + gains.Kp * e + I_new
    u = min(max(u_unsat, Q_min), Q_max)
    tol = SAT_TOL * max(abs(Q_max), 1.0)
    sat = 1 if u_unsat >= Q_max - tol else (-1 if u_unsat <= Q_min + tol else 0)
    return u, I_new, {"I_tent": I_tent, "v_tent": v, "u_unsat": u_unsat, "sat": sat, "limited": limited,
                      "blocked": limited and I_new == I_prev}


@dataclass
class FixedSetpointPI:
    gains: PIGains
    T_set: float
    Q_min: float = 0.0
    Q_max: float = 1.5e6
    Q_ff: float = 0.0
    anti_windup: bool = True
    I: float = 0.0  # integrator state, W_th

    def initialize(self, Q_init):
        """Bumpless start: with e = 0 the first command equals Q_init."""
        if not (self.Q_min <= Q_init <= self.Q_max):
            raise ValueError(f"initial cooling {Q_init:.4g} W outside [{self.Q_min}, {self.Q_max}]")
        self.I = Q_init - self.Q_ff

    def setpoint(self, P_IT):
        """B1: the setpoint does not depend on load."""
        return self.T_set

    def update(self, T_a, P_IT=None):
        self.T_set = self.setpoint(P_IT)
        e = T_a - self.T_set
        I_prev = self.I
        u, self.I, info = pi_update(I_prev, e, self.gains, self.Q_ff, self.Q_min, self.Q_max, self.anti_windup)
        info.update(e=e, I_prev=I_prev, I=self.I, u=u, T_set=self.T_set)
        return u, info


def adaptive_setpoint(P_IT, UA_s, T_a_max=27.0, T_s_max=40.0, margin=1.0, lo=18.0, hi=27.0):
    """B2 rule: T_set = clip(min(T_a_max, T_s_max - P_IT / UA_s) - margin, lo, hi), degC.

    Quasi-steady reasoning (spec Prop. 5): at steady state T_s - T_a = P_IT / UA_s. The margin (1 K)
    is provisional, not a safety guarantee. The clip only bounds the setpoint; under overload it
    cannot make temperatures safe.
    """
    if not (np.isfinite(P_IT) and P_IT >= 0):
        raise ValueError(f"P_IT must be finite and non-negative; got {P_IT}")
    return float(np.clip(min(T_a_max, T_s_max - P_IT / UA_s) - margin, lo, hi))


@dataclass
class LoadAdaptivePI(FixedSetpointPI):
    """B2: identical PI law, gains and anti-windup to B1; only the setpoint rule differs.

    T_set holds the setpoint currently in force (initialise it with adaptive_setpoint at the
    starting load). UA_s must be the plant's Params.UA_s, W/K.
    """
    UA_s: float = None
    T_a_max: float = 27.0
    T_s_max: float = 40.0
    margin: float = 1.0
    T_set_lo: float = 18.0
    T_set_hi: float = 27.0

    def __post_init__(self):
        if self.UA_s is None or not (np.isfinite(self.UA_s) and self.UA_s > 0):
            raise ValueError("LoadAdaptivePI needs UA_s (W/K) from Params")

    def setpoint(self, P_IT):
        if P_IT is None:
            raise ValueError("LoadAdaptivePI.update needs the measured IT load P_IT (W)")
        return adaptive_setpoint(P_IT, self.UA_s, self.T_a_max, self.T_s_max, self.margin,
                                 self.T_set_lo, self.T_set_hi)


# ---------------------------------------------------------------------------
# Disturbances
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PiecewiseDisturbance:
    """P_IT and T_out held constant on [times[i], times[i+1]); times[0] must be 0."""

    times: np.ndarray
    P_IT: np.ndarray
    T_out: np.ndarray

    def __post_init__(self):
        t = np.asarray(self.times, dtype=float)
        P = np.asarray(self.P_IT, dtype=float)
        To = np.asarray(self.T_out, dtype=float)
        if not (t.ndim == P.ndim == To.ndim == 1 and t.size == P.size == To.size and t.size >= 1):
            raise ValueError("times, P_IT and T_out must be 1-D arrays of equal length")
        if not (np.all(np.isfinite(t)) and np.all(np.isfinite(P)) and np.all(np.isfinite(To))):
            raise ValueError("disturbance values must be finite")
        if t[0] != 0.0 or np.any(np.diff(t) <= 0):
            raise ValueError("disturbance times must start at 0 and be strictly increasing")
        if np.any(P < 0):
            raise ValueError("P_IT must be non-negative")
        object.__setattr__(self, "times", t)
        object.__setattr__(self, "P_IT", P)
        object.__setattr__(self, "T_out", To)

    def at(self, t):
        i = int(np.searchsorted(self.times, t, side="right") - 1)
        return self.P_IT[i], self.T_out[i]

    def boundaries_in(self, a, b):
        """Disturbance change times strictly inside (a, b)."""
        return self.times[(self.times > a) & (self.times < b)]


# ---------------------------------------------------------------------------
# Closed-loop simulation with intersample inspection
# ---------------------------------------------------------------------------

@dataclass
class ClosedLoopResult:
    t: np.ndarray       # fine-grid times, s
    X: np.ndarray       # [T_s, T_a, E_Q, E_IT, E_env] at t
    Q: np.ndarray       # cooling applied on [t_i, t_{i+1}) (last entry repeats)
    P_IT: np.ndarray    # disturbance on [t_i, t_{i+1})
    T_out: np.ndarray
    log: dict           # per control update: arrays keyed by name
    status: str         # "completed" or "stopped_cop_domain"
    message: str = ""
    stop_time: float = np.nan
    meta: dict = field(default_factory=dict)
    SP: np.ndarray = None       # setpoint in force on [t_i, t_{i+1}), degC
    P_cool: np.ndarray = None   # cooling electricity Q_c / COP at the start of [t_i, t_{i+1}), W_e
    E_cool: np.ndarray = None   # cumulative cooling electricity at t_i, J_e (Simpson per segment)
    peak_P_cool: float = np.nan  # max of Q_c / COP over segment ends and midpoints, W_e


def _g_or_error(p, T_a, T_out):
    """(1/COP, None) inside the COP domain, else (None, error message)."""
    try:
        return 1.0 / float(tm.cop(p, T_a, T_out)), None
    except tm.COPDomainError as err:
        return None, str(err)


def simulate_closed_loop(p, ctrl, dist, x0, t_end, fine_dt=5.0, check_cop=True, inv_cop=None):
    """Simulate the PI loop on the exact ZOH plant.

    Each control interval [t_k, t_k + dt) is propagated exactly on a fine grid (every fine_dt
    seconds) that is also split at every disturbance boundary, so temperatures are available
    between updates. The controller receives T_a(t_k) and the measured load P_IT(t_k).

    If check_cop is True, cooling electricity is accounted on every segment where Q_c > 0:
    E_cool += Q_c * h/6 * (g(a) + 4 g(mid) + g(b)), g = 1/COP(T_a, T_out), with the midpoint
    state propagated exactly (Simpson's rule; Q_c and T_out are constant on the segment).
    If COP is outside its domain at a segment end (checked first, as in Stage 2) or midpoint,
    the run is stopped there and the error is recorded. Constant fan power is not included.

    inv_cop: optional replacement for g, called as inv_cop(p, T_a, T_out) -> (1/COP, error message or None).
    Default None uses thermal_model.cop. Added in Stage 4 only for the frozen-COP invariance test.
    """
    dt = ctrl.gains.dt
    x0 = np.asarray(x0, dtype=float)
    if x0.shape != (2,) or not np.all(np.isfinite(x0)):
        raise ValueError(f"x0 must be 2 finite values; got {x0}")
    if not (np.isfinite(fine_dt) and 0 < fine_dt <= dt):
        raise ValueError(f"fine_dt must be in (0, dt]; got {fine_dt}")
    n_upd = int(round(t_end / dt))
    if n_upd < 1 or abs(n_upd * dt - t_end) > 1e-9 * max(t_end, 1.0):
        raise ValueError(f"t_end ({t_end}) must be a positive multiple of dt ({dt})")

    g_of = _g_or_error if inv_cop is None else inv_cop
    cache = {}

    def step(h):
        key = round(h, 9)
        if key not in cache:
            cache[key] = tm.zoh_step_matrices(p, h)
        return cache[key]

    X = np.zeros(tm.N_X + tm.N_Z)
    X[:2] = x0
    T_set0 = ctrl.T_set
    ts, Xs, Qs, Ps, Tos, SPs, Pcs, Ecs = [0.0], [X.copy()], [], [], [], [], [], [0.0]
    log = {k: [] for k in ("t", "T_s", "T_a", "T_set", "e", "I_prev", "I_tent", "I", "v_tent", "u_unsat", "u",
                           "sat", "limited", "blocked")}
    status, message, stop_time = "completed", "", np.nan
    E_cool, peak = 0.0, 0.0 if check_cop else np.nan

    for k in range(n_upd):
        t_k = k * dt
        u, info = ctrl.update(X[1], dist.at(t_k)[0])
        for name, val in (("t", t_k), ("T_s", X[0]), ("T_a", X[1])):
            log[name].append(val)
        for name in ("T_set", "e", "I_prev", "I_tent", "I", "v_tent", "u_unsat", "u", "sat", "limited", "blocked"):
            log[name].append(info[name])

        t_next = t_k + dt
        grid = np.unique(np.concatenate([[t_k], np.arange(t_k, t_next, fine_dt)[1:],
                                         dist.boundaries_in(t_k, t_next), [t_next]]))
        seg = {key: [] for key in ("t", "X", "Q", "P", "To", "SP", "Pc", "Ec")}
        for a, b in zip(grid[:-1], grid[1:]):
            P, To = dist.at(a)
            w = np.array([u, P, To])
            X_a = X
            Phi, Gam = step(b - a)
            X = Phi @ X_a + Gam @ w
            Pc_a = 0.0 if check_cop else np.nan
            if check_cop and u > 0:
                g_a, err_a = g_of(p, X_a[1], To)
                g_b, err_b = g_of(p, X[1], To)
                if err_a or err_b:
                    status, message = "stopped_cop_domain", err_a or err_b
                    stop_time = a if err_a else b
                    Pc_a = np.nan
                else:
                    Phi_h, Gam_h = step((b - a) / 2)
                    X_m = Phi_h @ X_a + Gam_h @ w
                    g_m, err_m = g_of(p, X_m[1], To)
                    if err_m:
                        status, message, stop_time = "stopped_cop_domain", err_m, (a + b) / 2
                        b, X, Pc_a = stop_time, X_m, np.nan   # keep the invalid midpoint as the final point
                    else:
                        E_cool += u * (b - a) / 6.0 * (g_a + 4.0 * g_m + g_b)
                        peak = max(peak, u * g_a, u * g_m, u * g_b)
                        Pc_a = u * g_a
            if status != "completed" and stop_time == a:
                break
            for key, val in (("t", b), ("X", X.copy()), ("Q", u), ("P", P), ("To", To), ("SP", ctrl.T_set),
                             ("Pc", Pc_a), ("Ec", E_cool)):
                seg[key].append(val)
            if status != "completed":
                break
        ts += seg["t"]; Xs += seg["X"]; Qs += seg["Q"]; Ps += seg["P"]; Tos += seg["To"]
        SPs += seg["SP"]; Pcs += seg["Pc"]; Ecs += seg["Ec"]
        if status != "completed":
            break

    # input arrays align with ts: entry i applies on [t_i, t_{i+1}); last entry repeats
    for arr in (Qs, Ps, Tos, SPs, Pcs):
        arr.append(arr[-1] if arr else np.nan)
    return ClosedLoopResult(
        t=np.array(ts), X=np.array(Xs), Q=np.array(Qs), P_IT=np.array(Ps), T_out=np.array(Tos),
        log={k: np.array(v) for k, v in log.items()}, status=status, message=message, stop_time=stop_time,
        meta={"dt": dt, "fine_dt": fine_dt, "T_set": T_set0, "Kp": ctrl.gains.Kp, "Ki": ctrl.gains.Ki,
              "anti_windup": ctrl.anti_windup, "Q_max": ctrl.Q_max, "Q_min": ctrl.Q_min,
              "controller": type(ctrl).__name__},
        SP=np.array(SPs), P_cool=np.array(Pcs), E_cool=np.array(Ecs), peak_P_cool=peak)


def _cop_ok(p, T_a, T_out):
    try:
        tm.cop(p, T_a, T_out)
        return True
    except tm.COPDomainError:
        return False


# ---------------------------------------------------------------------------
# Metrics (all evaluated on the fine intersample grid)
# ---------------------------------------------------------------------------

def event_metrics(res, t_event, t_next, tol=0.1):
    """Disturbance-rejection metrics for the window [t_event, t_next] on the fine grid.

    The error is measured against the setpoint in force at each fine-grid point (constant for B1).
    settling_time: time after t_event after which |T_a - T_set| <= tol for the rest of the
    window (None if still outside tol at the window end). overshoot: largest excursion on the
    opposite side of the setpoint from the main (largest) deviation, occurring after it.
    """
    m = (res.t >= t_event) & (res.t <= t_next)
    if np.count_nonzero(m) < 2:
        return {"reached": False}
    t, Ta, Ts = res.t[m], res.X[m, 1], res.X[m, 0]
    if res.SP is not None:
        sp = res.SP[m].copy()
        # SP[i] is the setpoint applied on [t_i, t_{i+1}). A point at the window's end t_next belongs to
        # this window, so judge it against the setpoint in force on the segment that ends there, not the
        # (possibly new) setpoint chosen at t_next for the next window.
        if len(sp) > 1 and t[-1] >= t_next - 1e-9:
            sp[-1] = sp[-2]
    else:
        sp = res.meta["T_set"]
    e = Ta - sp
    i_peak = int(np.argmax(np.abs(e)))
    peak = e[i_peak]
    after = e[i_peak:]
    overshoot = float(max(0.0, np.max(-np.sign(peak) * after))) if peak != 0 else 0.0
    outside = np.nonzero(np.abs(e) > tol)[0]
    truncated = t[-1] < t_next - 1e-9
    if outside.size == 0:
        settling = 0.0
    elif outside[-1] == len(t) - 1 or truncated:
        settling = None
    else:
        settling = float(t[outside[-1] + 1] - t_event)
    last_hour = t >= t[-1] - 3600.0
    in_win = (res.log["t"] >= t_event) & (res.log["t"] < t_next)
    dt = res.meta["dt"]
    return {
        "reached": True, "truncated": bool(truncated),
        "peak_dev": float(peak), "overshoot": overshoot, "settling_time": settling,
        "final_err": float(e[-1]), "max_err_last_hour": float(np.max(np.abs(e[last_hour]))),
        "T_a_min": float(Ta.min()), "T_a_max": float(Ta.max()),
        "T_s_min": float(Ts.min()), "T_s_max": float(Ts.max()),
        "sat_high_time": float(np.sum(res.log["sat"][in_win] == 1) * dt),
        "sat_low_time": float(np.sum(res.log["sat"][in_win] == -1) * dt),
    }


def saturation_summary(res):
    """Time at each cooling limit (from the control log), longest episode, limited integrations.

    limited_updates: integration reduced by anti-windup (clamped at a limit or fully blocked);
    blocked_updates: the subset where the integrator did not change at all.
    """
    dt, sat = res.meta["dt"], res.log["sat"]

    def longest(flag):
        best = run = 0
        for s in sat:
            run = run + 1 if s == flag else 0
            best = max(best, run)
        return best * dt

    return {"high_time": float(np.sum(sat == 1) * dt), "low_time": float(np.sum(sat == -1) * dt),
            "longest_high": float(longest(1)), "longest_low": float(longest(-1)),
            "limited_updates": int(np.sum(res.log["limited"])), "blocked_updates": int(np.sum(res.log["blocked"])),
            "I_min": float(res.log["I"].min()), "I_max": float(res.log["I"].max())}


def limit_violations(res, T_a_band=(18.0, 27.0), T_s_max=40.0, tol=0.0):
    """Time and degree-hours outside the (illustrative) limits, trapezoid on the fine grid.

    A point counts as violating only if it is beyond a limit by more than tol (K), a declared
    numerical tolerance. max_* are the largest exceedances in K (0 if none).
    """
    t, Ta, Ts = res.t, res.X[:, 1], res.X[:, 0]

    def dh(excess):
        return float(np.trapezoid(np.maximum(excess, 0.0), t) / 3600.0)

    def dur(mask):
        return float(np.sum(np.diff(t)[mask[:-1] | mask[1:]]))

    over_a, under_a, over_s = Ta - T_a_band[1], T_a_band[0] - Ta, Ts - T_s_max
    return {
        "T_a_above_max_h": dur(over_a > tol) / 3600, "T_a_above_max_Kh": dh(over_a),
        "T_a_below_min_h": dur(under_a > tol) / 3600, "T_a_below_min_Kh": dh(under_a),
        "T_s_above_max_h": dur(over_s > tol) / 3600, "T_s_above_max_Kh": dh(over_s),
        "T_a_max": float(Ta.max()), "T_a_min": float(Ta.min()), "T_s_max": float(Ts.max()),
        "max_T_a_over": float(max(over_a.max(), 0.0)), "max_T_a_under": float(max(under_a.max(), 0.0)),
        "max_T_s_over": float(max(over_s.max(), 0.0)),
        "within_limits": bool(over_a.max() <= tol and under_a.max() <= tol and over_s.max() <= tol),
    }


def energy_summary(p, res):
    """Energy bookkeeping for one run (units in the keys). Fan power is not included.

    dU = C_s (T_s,end - T_s,start) + C_a (T_a,end - T_a,start): stored heat gained over the run.
    For a stopped run all values cover only the part before the stop (cooling electricity up to
    the last fully valid segment) and must not be read as a complete-horizon energy.
    """
    x0, x1 = res.X[0, :2], res.X[-1, :2]
    E_Q, E_IT = res.X[-1, 2] - res.X[0, 2], res.X[-1, 3] - res.X[0, 3]
    E_cool = res.E_cool[-1] if res.E_cool is not None else np.nan
    dU = p.C_s * (x1[0] - x0[0]) + p.C_a * (x1[1] - x0[1])
    return {
        "complete": res.status == "completed",
        "E_cool_kWh_e": float(tm.J_to_kWh(E_cool)), "E_Q_kWh_th": float(tm.J_to_kWh(E_Q)),
        "E_IT_kWh_e": float(tm.J_to_kWh(E_IT)), "cool_to_IT_ratio": float(E_cool / E_IT) if E_IT > 0 else np.nan,
        "peak_P_cool_kW_e": float(res.peak_P_cool / 1e3), "dU_J": float(dU), "dU_kWh_th": float(tm.J_to_kWh(dU)),
        "dT_s_end_K": float(x1[0] - x0[0]), "dT_a_end_K": float(x1[1] - x0[1]),
    }
