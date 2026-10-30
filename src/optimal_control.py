"""Best-found perfect-foresight cooling benchmark and LP lower bound (PROJECT_SPEC.md §7.4, Stage 4).

OFFLINE BENCHMARK, NOT A CONTROLLER: the optimiser is given the whole future disturbance trace
(P_IT and T_out), which no real controller knows. Results are "best-found", never "globally optimal".

Discretised problem P_d on the update grid t_k = k dt, k = 0..N (spec §7.4):
    variables   u_k (k = 0..N-1, cooling command held on [t_k, t_k+1)),  x_k = [T_s,k, T_a,k] (k = 0..N)
    dynamics    x_{k+1} = Phi x_k + Gamma [u_k, P_IT,k, T_out,k]     exact ZOH (thermal_model.zoh_step_matrices)
    boundary    x_0 = x_start,  x_N = x_end (closed cycle: x_end = x_start)
    bounds      0 <= u_k <= Q_max;  18 <= T_a,k <= 27;  T_s,k <= 40   at the sample points only
    objective   J_d = sum_k dt u_k g(T_a,k, T_out,k),  g = 1/COP      (left-endpoint rule, J_e; reported in kWh_e)
Optional variant (separately labelled when used): the same limits also imposed at intersample points
of every interval, which are exact linear functions of (x_k, u_k).

LP lower bound: replace g(T_a,k, T_out,k) by g(T_a,max, T_out,k). Because g is strictly decreasing in T_a
(dg/dT_a = -T_cd / (eta_II T_ev^2) < 0 on the COP domain) and every feasible point has u_k >= 0 and
T_a,k <= T_a,max, u_k g(T_a,k) >= u_k g(T_a,max) term by term, so J_LP(z) <= J_d(z) for every feasible z,
hence min J_LP <= min J_d. The bound concerns J_d only, not continuous-time electricity.

Decision vector layout (shared by the NLP and the LP):  z = [u_0..u_{N-1} (MW_th), T_s,0..T_s,N, T_a,0..T_a,N]
u is scaled to MW_th inside z only; every function that returns commands returns W_th.
"""

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

import thermal_model as tm
import controllers as ct

U_SCALE = 1e6  # W_th per unit of u in z


def inv_cop_expr(p, T_a, T_out):
    """1/COP written with arithmetic only, so it works on numpy arrays and CasADi symbols.

    Same formula as thermal_model.cop but WITHOUT domain checks; callers must stay inside the domain
    (the NLP does, because T_a is bounded to the band; see check_bound_coefficients).
    """
    T_ev = T_a - p.dT_ev + tm.KELVIN_OFFSET
    T_cd = T_out + p.dT_cd + tm.KELVIN_OFFSET
    return (T_cd - T_ev) / (p.eta_II * T_ev)


@dataclass
class DiscreteProblem:
    """P_d for one scenario, update interval and boundary protocol. Built by build_problem()."""
    p: tm.Params
    dt: float
    N: int
    P_IT: np.ndarray         # (N,) W, held on [t_k, t_k+1)
    T_out: np.ndarray        # (N,) degC
    x_start: np.ndarray      # (2,) [T_s, T_a]
    x_end: np.ndarray        # (2,)
    T_a_band: tuple
    T_s_max: float
    Q_max: float
    Phi: np.ndarray          # (2, 2)
    Gam: np.ndarray          # (2, 3) columns [Q_c, P_IT, T_out]
    A_eq: sp.csr_matrix      # dynamics + boundary rows
    b_eq: np.ndarray
    lb: np.ndarray           # variable bounds on z
    ub: np.ndarray
    A_in: sp.csr_matrix = None   # optional intersample rows: lb_in <= A_in z <= ub_in
    lb_in: np.ndarray = None
    ub_in: np.ndarray = None
    c_in: np.ndarray = None      # disturbance offset of each intersample row: state value = A_in z + c_in
    t_in: np.ndarray = None      # time of each intersample row, s
    which_in: np.ndarray = None  # 0 = T_s row, 1 = T_a row
    intersample_dt: float = None
    label: str = ""

    @property
    def n(self):
        return 3 * self.N + 2

    def split(self, z):
        """z -> (u in W_th, x of shape (N+1, 2))."""
        N = self.N
        z = np.asarray(z, dtype=float)
        return z[:N] * U_SCALE, np.column_stack([z[N:2 * N + 1], z[2 * N + 1:]])

    def join(self, u_W, x):
        x = np.asarray(x, dtype=float)
        return np.concatenate([np.asarray(u_W, dtype=float) / U_SCALE, x[:, 0], x[:, 1]])

    @property
    def t(self):
        return np.arange(self.N + 1) * self.dt


def build_problem(p, dist, x_start, t_end, dt, x_end=None, T_a_band=(18.0, 27.0), T_s_max=40.0,
                  intersample_dt=None, label=""):
    """Assemble P_d. dist is a controllers.PiecewiseDisturbance whose changes lie on the update grid.

    If intersample_dt is given (it must divide dt), the band and T_s limits are also imposed at
    t_k + j*intersample_dt, j = 1..dt/intersample_dt - 1, inside every interval (a separate variant).
    """
    N = int(round(t_end / dt))
    if N < 1 or abs(N * dt - t_end) > 1e-9 * t_end:
        raise ValueError(f"t_end ({t_end}) must be a positive multiple of dt ({dt})")
    off = dist.times[(dist.times > 0) & (dist.times < t_end)]
    if np.any(np.abs(off / dt - np.round(off / dt)) > 1e-9):
        raise ValueError("disturbance changes must lie on the update grid (ZOH holds them constant per interval)")
    tk = np.arange(N) * dt
    PT = np.array([dist.at(t) for t in tk])
    P_IT, T_out = PT[:, 0], PT[:, 1]
    x_start = np.asarray(x_start, dtype=float)
    x_end = x_start.copy() if x_end is None else np.asarray(x_end, dtype=float)
    Phi, Gam = tm.zoh_step_matrices(p, dt)
    Phi, Gam = Phi[:2, :2], Gam[:2, :]

    n = 3 * N + 2
    iu = np.arange(N)
    i_s = N + np.arange(N + 1)
    i_a = 2 * N + 1 + np.arange(N + 1)
    rows, cols, vals, b = [], [], [], []

    def add(r, c, v):
        rows.append(r); cols.append(c); vals.append(v)

    r = 0
    for j, idx in enumerate((i_s[0], i_a[0])):            # initial state
        add(r, idx, 1.0); b.append(x_start[j]); r += 1
    for k in range(N):                                     # x_{k+1} - Phi x_k - Gam_u u_k = Gam_d d_k
        for i, nxt in enumerate((i_s, i_a)):
            add(r, nxt[k + 1], 1.0)
            add(r, i_s[k], -Phi[i, 0])
            add(r, i_a[k], -Phi[i, 1])
            add(r, iu[k], -Gam[i, 0] * U_SCALE)
            b.append(Gam[i, 1] * P_IT[k] + Gam[i, 2] * T_out[k])
            r += 1
    for j, idx in enumerate((i_s[N], i_a[N])):            # terminal state
        add(r, idx, 1.0); b.append(x_end[j]); r += 1
    A_eq = sp.csr_matrix((vals, (rows, cols)), shape=(r, n))

    lb = np.full(n, -np.inf); ub = np.full(n, np.inf)
    lb[iu], ub[iu] = 0.0, p.Q_max / U_SCALE
    ub[i_s] = T_s_max
    lb[i_a], ub[i_a] = T_a_band

    prob = DiscreteProblem(p=p, dt=float(dt), N=N, P_IT=P_IT, T_out=T_out, x_start=x_start, x_end=x_end,
                           T_a_band=tuple(T_a_band), T_s_max=float(T_s_max), Q_max=p.Q_max, Phi=Phi, Gam=Gam,
                           A_eq=A_eq, b_eq=np.array(b), lb=lb, ub=ub, label=label)
    if intersample_dt is not None:
        m = int(round(dt / intersample_dt))
        if m < 2 or abs(m * intersample_dt - dt) > 1e-9 * dt:
            raise ValueError("intersample_dt must divide dt into at least 2 parts")
        sub = [tm.zoh_step_matrices(p, j * intersample_dt) for j in range(1, m)]
        rows, cols, vals, lo, hi, cs, ts, wh = [], [], [], [], [], [], [], []
        r = 0
        for k in range(N):
            for j, (Ph, Gm) in enumerate(sub, start=1):
                Ph, Gm = Ph[:2, :2], Gm[:2, :]
                for i, (low, high) in enumerate(((-np.inf, T_s_max), T_a_band)):
                    c = Gm[i, 1] * P_IT[k] + Gm[i, 2] * T_out[k]
                    for col, v in ((i_s[k], Ph[i, 0]), (i_a[k], Ph[i, 1]), (iu[k], Gm[i, 0] * U_SCALE)):
                        rows.append(r); cols.append(col); vals.append(v)
                    lo.append(low - c); hi.append(high - c)
                    cs.append(c); ts.append(k * dt + j * intersample_dt); wh.append(i)
                    r += 1
        prob.A_in = sp.csr_matrix((vals, (rows, cols)), shape=(r, n))
        prob.lb_in, prob.ub_in = np.array(lo), np.array(hi)
        prob.c_in, prob.t_in, prob.which_in = np.array(cs), np.array(ts, dtype=float), np.array(wh)
        prob.intersample_dt = float(intersample_dt)
    return prob


# ---------------------------------------------------------------------------
# Independent checks on a candidate point (do not use A_eq, so they test the assembly too)
# ---------------------------------------------------------------------------

def check_point(prob, u_W, x):
    """Residuals of (u, x) against P_d, computed directly from Phi, Gamma and the limits.

    Returns K for temperatures and W for commands. Intersample rows (variant) are checked via A_in.
    """
    u_W, x = np.asarray(u_W, dtype=float), np.asarray(x, dtype=float)
    W = np.column_stack([u_W, prob.P_IT, prob.T_out])
    dyn = x[1:] - (x[:-1] @ prob.Phi.T + W @ prob.Gam.T)
    out = {
        "dyn_K": float(np.max(np.abs(dyn))),
        "init_K": float(np.max(np.abs(x[0] - prob.x_start))),
        "term_K": float(np.max(np.abs(x[-1] - prob.x_end))),
        "u_viol_W": float(max(0.0, -u_W.min(), u_W.max() - prob.Q_max)),
        "Ta_viol_K": float(max(0.0, prob.T_a_band[0] - x[:, 1].min(), x[:, 1].max() - prob.T_a_band[1])),
        "Ts_viol_K": float(max(0.0, x[:, 0].max() - prob.T_s_max)),
        "inter_viol_K": 0.0,
    }
    if prob.A_in is not None:
        a = prob.A_in @ prob.join(u_W, x)
        out["inter_viol_K"] = float(max(0.0, np.max(prob.lb_in - a), np.max(a - prob.ub_in)))
    return out


def is_feasible(chk, tol_K=1e-6, tol_W=1e-3):
    return (max(chk["dyn_K"], chk["init_K"], chk["term_K"], chk["Ta_viol_K"], chk["Ts_viol_K"],
                chk["inter_viol_K"]) <= tol_K and chk["u_viol_W"] <= tol_W)


def J_d(prob, u_W, T_a, frozen_cop=None):
    """Discrete objective sum dt u_k / COP(T_a,k, T_out,k) in kWh_e, using thermal_model.cop (domain-checked)."""
    T_a = np.asarray(T_a, dtype=float)[:prob.N]
    g = 1.0 / frozen_cop if frozen_cop is not None else 1.0 / tm.cop(prob.p, T_a, prob.T_out)
    return float(tm.J_to_kWh(np.sum(prob.dt * np.asarray(u_W) * g)))


def bound_coefficients(prob, frozen_cop=None):
    """Per-interval LP cost coefficient g_LB,k = min over the band of 1/COP = 1/COP(T_a,max, T_out,k)."""
    if frozen_cop is not None:
        return np.full(prob.N, 1.0 / frozen_cop)
    return 1.0 / tm.cop(prob.p, np.full(prob.N, prob.T_a_band[1]), prob.T_out)


def check_bound_coefficients(prob, n_grid=9001):
    """Verify, on a dense grid over the whole air band at every outdoor temperature in the trace, that
    cop() is inside its domain and that 1/COP(T_a,max) is the minimum (the LP coefficient is a lower bound).

    Returns (ok, info). Also returns the analytic derivative sign check dg/dT_a = -T_cd/(eta T_ev^2).
    """
    T = np.linspace(prob.T_a_band[0], prob.T_a_band[1], n_grid)
    worst, ok, msg = np.inf, True, ""
    max_dg = -np.inf
    for To in np.unique(prob.T_out):
        try:
            g = 1.0 / tm.cop(prob.p, T, np.full_like(T, To))
        except tm.COPDomainError as err:
            return False, {"message": f"T_out {To}: {err}"}
        g_lb = 1.0 / tm.cop(prob.p, prob.T_a_band[1], To)
        worst = min(worst, float(np.min(g - g_lb)))
        T_ev = T - prob.p.dT_ev + tm.KELVIN_OFFSET
        dg = -(To + prob.p.dT_cd + tm.KELVIN_OFFSET) / (prob.p.eta_II * T_ev ** 2)
        max_dg = max(max_dg, float(dg.max()))
        if np.any(np.diff(g) >= 0):
            ok, msg = False, f"1/COP not strictly decreasing at T_out {To}"
    ok = ok and worst >= 0.0 and max_dg < 0
    return ok, {"min_g_minus_gLB": worst, "max_dg_dTa": max_dg, "T_out_values": np.unique(prob.T_out), "message": msg}


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------

@dataclass
class SolveResult:
    label: str
    solver: str
    status: str            # solver's own return status text
    solver_success: bool   # solver's own success flag
    time_s: float          # wall time of the solve call
    iterations: int
    objective: float       # solver-reported objective, kWh_e (J_d for NLP, J_LP for LP)
    z: np.ndarray
    u_W: np.ndarray
    x: np.ndarray
    check: dict            # check_point residuals
    feasible: bool         # independent feasibility (is_feasible)
    J_d: float = np.nan    # independently re-evaluated J_d at the returned point (NaN if outside COP domain)
    data: dict = field(default_factory=dict)  # the exact arrays passed to the solver (for the shared-data check)


def _finish(prob, label, solver, status, success, t, iters, obj, z, data, frozen_cop, tol_K, tol_W):
    u, x = prob.split(z)
    chk = check_point(prob, u, x)
    try:
        jd = J_d(prob, u, x[:, 1], frozen_cop)
    except tm.COPDomainError:
        jd = np.nan
    return SolveResult(label, solver, status, bool(success), t, int(iters), float(obj), np.asarray(z), u, x, chk,
                       is_feasible(chk, tol_K, tol_W) and np.isfinite(jd), jd, data)


def implied_T_s_lower(prob):
    """Finite lower bounds L_k <= T_s,k valid at every feasible point of P_d (T_s has no lower limit of its own).

    From T_s,k+1 = Phi_00 T_s,k + Phi_01 T_a,k + Gam_00 u_k + Gam_01 P_k + Gam_02 T_out,k with Phi_00, Phi_01 >= 0 and
    Gam_00 <= 0 (checked): the smallest value uses T_s,k = L_k, T_a,k = T_a,min and u_k = Q_max. Used only to make the
    dual certificate finite.
    """
    Phi, Gam = prob.Phi, prob.Gam
    if not (Phi[0, 0] >= 0 and Phi[0, 1] >= 0 and Gam[0, 0] <= 0):
        raise ValueError("sign pattern needed for the implied T_s bound does not hold")
    L = np.zeros(prob.N + 1)
    L[0] = prob.x_start[0]
    for k in range(prob.N):
        L[k + 1] = (Phi[0, 0] * L[k] + Phi[0, 1] * prob.T_a_band[0] + Gam[0, 0] * prob.Q_max
                    + Gam[0, 1] * prob.P_IT[k] + Gam[0, 2] * prob.T_out[k])
    return L


def dual_bound(prob, c, y_eq, y_ub=None, A_ub=None, b_ub=None):
    """Lagrangian (weak-duality) lower bound on min c'z over P_d's constraints, for ANY multipliers.

    With r = c - A_eq' y_eq - A_ub' y_ub (y_ub <= 0 enforced by clipping), every feasible z satisfies
    c'z >= b_eq' y_eq + b_ub' y_ub + sum_i min(r_i lo_i, r_i hi_i), with lo/hi the variable bounds (T_s lower
    bounds from implied_T_s_lower). This certifies the LP bound independently of the solver's status text.
    """
    lo, hi = prob.lb.copy(), prob.ub.copy()
    lo[prob.N:2 * prob.N + 1] = implied_T_s_lower(prob)
    if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
        raise ValueError("dual bound needs finite variable bounds")
    r = c - prob.A_eq.T @ y_eq
    val = float(prob.b_eq @ y_eq)
    if A_ub is not None:
        y_ub = np.minimum(np.asarray(y_ub, dtype=float), 0.0)
        r = r - A_ub.T @ y_ub
        val += float(b_ub @ y_ub)
    return val + float(np.sum(np.minimum(r * lo, r * hi)))


# HiGHS configurations tried for every LP, in this order. The SciPy default ("highs", dual simplex with
# presolve) returned status 4 / "Not Set" on these LPs during development (SciPy 1.18.1), so every
# configuration is run and logged; the bound is certified by dual_bound, not by any status text.
LP_METHODS = (("highs (default)", "highs", {}),
              ("highs-ipm", "highs-ipm", {}),
              ("highs-ds, no presolve", "highs-ds", {"presolve": False}))


def solve_lp(prob, frozen_cop=None, time_limit=600.0, tol_K=1e-6, tol_W=1e-3):
    """LP lower bound with HiGHS (scipy.optimize.linprog) on exactly prob's constraints.

    Returns (best, attempts): attempts lists every configuration in LP_METHODS with its status, time,
    objective and certified dual bound; best is the optimal attempt with the largest certified bound
    (None if no attempt is optimal). best.data["dual_bound"] is the bound to report.
    """
    c = np.zeros(prob.n)
    c[:prob.N] = prob.dt * U_SCALE * bound_coefficients(prob, frozen_cop) / tm.J_PER_KWH
    A_ub = b_ub = None
    if prob.A_in is not None:
        fin_hi, fin_lo = np.isfinite(prob.ub_in), np.isfinite(prob.lb_in)
        A_ub = sp.vstack([prob.A_in[fin_hi], -prob.A_in[fin_lo]]).tocsr()
        b_ub = np.concatenate([prob.ub_in[fin_hi], -prob.lb_in[fin_lo]])
    bounds = np.column_stack([prob.lb, prob.ub])
    attempts = []
    for name, method, extra in LP_METHODS:
        t0 = time.perf_counter()
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=prob.A_eq, b_eq=prob.b_eq, bounds=bounds, method=method,
                      options={"time_limit": time_limit, **extra})
        t = time.perf_counter() - t0
        ok = res.status == 0 and res.x is not None
        z = res.x if ok else np.full(prob.n, np.nan)
        data = {"A_eq": prob.A_eq, "b_eq": prob.b_eq, "lb": prob.lb, "ub": prob.ub, "A_in": prob.A_in, "c": c,
                "method": name}
        out = _finish(prob, prob.label, f"HiGHS via scipy linprog [{name}]", res.message, ok, t,
                      getattr(res, "nit", -1), res.fun if ok else np.nan, z, data, frozen_cop, tol_K, tol_W)
        if ok:
            out.data["dual_bound"] = dual_bound(prob, c, res.eqlin.marginals,
                                                None if A_ub is None else res.ineqlin.marginals, A_ub, b_ub)
            out.data["eq_residual"] = float(np.max(np.abs(res.eqlin.residual)))
        else:
            out.data["dual_bound"] = np.nan
        attempts.append(out)
    good = [a for a in attempts if a.solver_success and np.isfinite(a.data["dual_bound"])]
    best = max(good, key=lambda a: a.data["dual_bound"]) if good else None
    return best, attempts


def make_nlp(prob, frozen_cop=None):
    """CasADi NLP for P_d built from prob.A_eq / bounds (and A_in for the variant). Returns (nlp dict, data)."""
    import casadi as ca
    N = prob.N
    z = ca.SX.sym("z", prob.n)
    u = z[:N] * U_SCALE
    T_a = z[2 * N + 1:2 * N + 1 + N]
    g = (1.0 / frozen_cop) if frozen_cop is not None else inv_cop_expr(prob.p, T_a, ca.DM(prob.T_out))
    f = ca.sum1(prob.dt * u * g) / tm.J_PER_KWH
    A = prob.A_eq if prob.A_in is None else sp.vstack([prob.A_eq, prob.A_in]).tocsr()
    lbg = prob.b_eq if prob.A_in is None else np.concatenate([prob.b_eq, prob.lb_in])
    ubg = prob.b_eq if prob.A_in is None else np.concatenate([prob.b_eq, prob.ub_in])
    Acsc = A.tocsc()
    Acsc.sort_indices()
    # CasADi stores matrices in compressed-column form: the data must be in CSC order
    A_ca = ca.DM(ca.Sparsity(A.shape[0], A.shape[1], Acsc.indptr.tolist(), Acsc.indices.tolist()), Acsc.data)
    gexpr = ca.mtimes(A_ca, z)
    nlp = {"x": z, "f": f, "g": gexpr}
    data = {"A_eq": prob.A_eq, "b_eq": prob.b_eq, "lb": prob.lb, "ub": prob.ub, "A_in": prob.A_in,
            "A_all": A, "lbg": lbg, "ubg": ubg, "g_fun": ca.Function("constraints", [z], [gexpr]),
            "jac_fun": ca.Function("constraint_jacobian", [z], [ca.jacobian(gexpr, z)]), "f_fun": ca.Function("objective", [z], [f])}
    return nlp, data


def solve_nlp(prob, z0, frozen_cop=None, max_wall_time=120.0, max_iter=3000, tol=1e-9, label=None,
              tol_K=1e-6, tol_W=1e-3, nlp=None):
    """Local NLP solve with IPOPT (CasADi, exact sparse Hessian, MUMPS). Returns SolveResult.

    z0 is projected into the variable bounds by IPOPT itself (bound_push); it need not be feasible.
    """
    import casadi as ca
    if nlp is None:
        nlp = make_nlp(prob, frozen_cop)
    expr, data = nlp
    opts = {"print_time": False, "ipopt.print_level": 0, "ipopt.sb": "yes", "ipopt.tol": tol,
            "ipopt.constr_viol_tol": 1e-9, "ipopt.max_iter": max_iter, "ipopt.max_wall_time": max_wall_time,
            "ipopt.linear_solver": "mumps", "ipopt.bound_relax_factor": 0.0}
    S = ca.nlpsol("oc", "ipopt", expr, opts)
    lbx = np.where(np.isfinite(prob.lb), prob.lb, -np.inf)
    ubx = np.where(np.isfinite(prob.ub), prob.ub, np.inf)
    t0 = time.perf_counter()
    sol = S(x0=np.asarray(z0, dtype=float), lbx=lbx, ubx=ubx, lbg=data["lbg"], ubg=data["ubg"])
    t = time.perf_counter() - t0
    st = S.stats()
    out = _finish(prob, label or prob.label, f"IPOPT {ca.__version__} (CasADi, MUMPS)", st["return_status"],
                  st["success"], t, st.get("iter_count", -1), float(sol["f"]), np.array(sol["x"]).ravel(), data,
                  frozen_cop, tol_K, tol_W)
    return out


def perturbed_start(prob, u_W, seed, rel=0.10):
    """Reproducible perturbed start: u + rel*Q_max*N(0,1) (seeded), clipped to [0, Q_max]; states obtained
    by propagating the dynamics from x_start (so they satisfy the dynamics but generally not the limits
    or the terminal condition)."""
    rng = np.random.default_rng(seed)
    u = np.clip(np.asarray(u_W) + rel * prob.Q_max * rng.standard_normal(prob.N), 0.0, prob.Q_max)
    return prob.join(u, propagate(prob, u))


def propagate(prob, u_W):
    """States x_0..x_N under commands u_W from x_start with the exact ZOH matrices of prob."""
    x = np.zeros((prob.N + 1, 2))
    x[0] = prob.x_start
    for k in range(prob.N):
        x[k + 1] = prob.Phi @ x[k] + prob.Gam @ np.array([u_W[k], prob.P_IT[k], prob.T_out[k]])
    return x


def nearest_feasible_start(prob, u_ref_W, time_limit=600.0):
    """Feasible initialisation when the reference (PI) sequence is infeasible for P_d:
    the point of P_d's feasible set that minimises sum_k |u_k - u_ref,k| (an LP, solved with HiGHS).
    Returns (z, info); info records which LP_METHODS configuration succeeded and every attempt's status."""
    N, n = prob.N, prob.n
    # variables [z (n), s (N)] with s_k >= |u_k - u_ref,k|
    c = np.concatenate([np.zeros(n), np.ones(N)])
    E = sp.hstack([sp.eye(N), sp.csr_matrix((N, n - N))])
    A_ub = sp.vstack([sp.hstack([E, -sp.eye(N)]), sp.hstack([-E, -sp.eye(N)])]).tocsr()
    b_ub = np.concatenate([u_ref_W / U_SCALE, -u_ref_W / U_SCALE])
    if prob.A_in is not None:
        fin_hi, fin_lo = np.isfinite(prob.ub_in), np.isfinite(prob.lb_in)
        A_ub = sp.vstack([A_ub, sp.hstack([prob.A_in[fin_hi], sp.csr_matrix((fin_hi.sum(), N))]),
                          sp.hstack([-prob.A_in[fin_lo], sp.csr_matrix((fin_lo.sum(), N))])]).tocsr()
        b_ub = np.concatenate([b_ub, prob.ub_in[fin_hi], -prob.lb_in[fin_lo]])
    A_eq = sp.hstack([prob.A_eq, sp.csr_matrix((prob.A_eq.shape[0], N))]).tocsr()
    bounds = np.vstack([np.column_stack([prob.lb, prob.ub]), np.column_stack([np.zeros(N), np.full(N, np.inf)])])
    log = []
    for name, method, extra in LP_METHODS:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=prob.b_eq, bounds=bounds, method=method,
                      options={"time_limit": time_limit, **extra})
        log.append((name, res.status, res.message))
        if res.status == 0:
            return res.x[:n], {"method": name, "attempts": log, "l1_MW": float(res.fun)}
    raise RuntimeError(f"no feasible initialisation found: {log}")


# ---------------------------------------------------------------------------
# Replay through the Stage 3 simulator
# ---------------------------------------------------------------------------

@dataclass
class CommandReplay:
    """Open-loop replay of a fixed command sequence through controllers.simulate_closed_loop, so the
    optimised commands see exactly the plant, fine grid and electricity accounting used for the PI runs.
    Measurements are ignored. Commands must already be within [Q_min, Q_max] (to 1e-3 W)."""
    commands: np.ndarray   # W_th, one per update
    dt: float
    Q_min: float = 0.0
    Q_max: float = 1.5e6

    def __post_init__(self):
        c = np.asarray(self.commands, dtype=float)
        if np.any(c < self.Q_min - 1e-3) or np.any(c > self.Q_max + 1e-3) or not np.all(np.isfinite(c)):
            raise ValueError("replayed commands must be finite and inside [Q_min, Q_max]")
        self.commands = np.clip(c, self.Q_min, self.Q_max)
        self.max_clip_W = float(np.max(np.abs(self.commands - c)))
        self.gains = ct.PIGains(0.0, 0.0, float(self.dt))
        self.T_set = np.nan
        self.anti_windup = False
        self.k = 0

    def update(self, T_a, P_IT=None):
        u = float(self.commands[self.k])
        self.k += 1
        tol = ct.SAT_TOL * max(abs(self.Q_max), 1.0)
        sat = 1 if u >= self.Q_max - tol else (-1 if u <= self.Q_min + tol else 0)
        nan = np.nan
        return u, {"T_set": nan, "e": nan, "I_prev": nan, "I_tent": nan, "I": nan, "v_tent": u, "u_unsat": u,
                   "u": u, "sat": sat, "limited": False, "blocked": False}


def replay(p, commands_W, dt, dist, x0, t_end, fine_dt=5.0, inv_cop=None):
    """Replay commands through the plant on the Stage 3 grid. Returns a controllers.ClosedLoopResult."""
    c = CommandReplay(commands_W, dt, 0.0, p.Q_max)
    res = ct.simulate_closed_loop(p, c, dist, x0, t_end, fine_dt=fine_dt, inv_cop=inv_cop)
    res.meta["x0"], res.meta["max_clip_W"] = np.asarray(x0, dtype=float), c.max_clip_W
    return res
