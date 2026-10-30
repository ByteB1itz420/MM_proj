"""Export the Stage 4 best-found OC solutions to artifacts/ for the main notebook.

Provenance: notebook/04_perfect_foresight_benchmark.ipynb is NOT modified. Its own code cells for setup, PI runs,
problem building/verification, the 300 s and 60 s solves and the replay (code cells 1-7) are executed unchanged in a
fresh in-memory kernel (cwd notebook/), then one extra export cell saves the best-found arrays. The exported numbers
are then compared with the values printed in notebook 04's saved outputs (best-found J_d and certified J_LB to the
printed 5 decimals, replayed E_cool to the printed 0.1 kWh_e). Any mismatch aborts the export.

Run from the project root:  python scripts/export_oc_artifacts.py      (about 15-20 min; IPOPT/HiGHS solves)
"""
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parent.parent
NB4 = ROOT / "notebook" / "04_perfect_foresight_benchmark.ipynb"
OUT = ROOT / "artifacts"
N_CODE_CELLS = 7   # setup, PI runs, problems + checks, solve definitions, 300 s solves, 60 s solves, replay

EXPORT = r'''
import json, platform, datetime, hashlib
import numpy as np
arrays, meta = {}, {"cases": {}}
for key, c in CASES.items():
    sid, dt, kind, var = key
    tag = f"{sid}_{dt}_{kind}_{var}"
    b, prob = c["best"], c["prob"]
    arrays[tag + "_u_W"] = b.u_W
    arrays[tag + "_x"] = b.x
    arrays[tag + "_x_start"] = prob.x_start
    r5 = REP[key]
    meta["cases"][tag] = {
        "scenario": sid, "dt_s": dt, "boundary_protocol": kind, "variant": var, "N": prob.N,
        "intersample_dt_s": prob.intersample_dt,
        "best_start": b.label, "best_status": b.status, "best_iterations": b.iterations, "best_time_s": b.time_s,
        "J_d_best_kWh_e": b.J_d, "J_d_solver_kWh_e": b.objective,
        "J_LB_certified_kWh_e": c["lp"].data["dual_bound"], "J_LP_primal_kWh_e": c["lp"].objective,
        "LP_method": c["lp"].data["method"],
        "LP_attempts": [[a.data["method"], a.status, bool(a.solver_success)] for a in c["lp_att"]],
        "starts": [[s.label, s.status, float(s.J_d), bool(s.feasible), float(s.time_s)] for s in c["sols"]],
        "warm_start": c["warm"], "pi_feasible_for_problem": bool(c["pi_feas"]), "J_d_PI_kWh_e": c["J_pi"],
        "replay_E_cool_5s_kWh_e": float(r5.E_cool[-1] / 3.6e6),
        "replay_status": r5.status,
    }
import scipy, casadi, sys as _s
meta["provenance"] = {
    "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    "source_notebook": "notebook/04_perfect_foresight_benchmark.ipynb (code cells 1-7 executed unchanged)",
    "python": _s.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__, "casadi": casadi.__version__,
    "platform": platform.platform(),
    "tolerances": {"TOL_FEAS_K": TOL_FEAS_K, "TOL_FEAS_W": TOL_FEAS_W, "NLP_WALL_S": NLP_WALL_S, "NLP_MAX_ITER": NLP_MAX_ITER,
                   "SEEDS": list(SEEDS), "PERTURB": PERTURB, "ISC_DT": ISC_DT, "HORIZON_s": HORIZON},
    "units": {"u_W": "W_th (command held on [t_k, t_k+dt))", "x": "degC, columns [T_s, T_a] at t_k, k = 0..N",
              "x_start": "degC [T_s, T_a]"},
}
np.savez_compressed("../artifacts/oc_best_found.npz", **arrays)
json.dump(meta, open("../artifacts/oc_best_found.json", "w"), indent=1)
print("EXPORTED", len(meta["cases"]), "cases")
'''


def main():
    nb = nbformat.read(NB4, as_version=4)
    code = [c for c in nb.cells if c.cell_type == "code"]
    src_hash = hashlib.sha256("\n".join(c.source for c in code[:N_CODE_CELLS]).encode()).hexdigest()
    printed = "\n".join(o.get("text", "") for c in code for o in c.get("outputs", []))
    run_nb = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(c.source) for c in code[:N_CODE_CELLS]]
                                      + [nbformat.v4.new_code_cell(EXPORT)], metadata=nb.metadata)
    t0 = time.perf_counter()
    NotebookClient(run_nb, timeout=None, kernel_name="python3",
                   resources={"metadata": {"path": str(ROOT / "notebook")}}).execute()
    print(f"executed notebook 04 code cells 1-{N_CODE_CELLS} + export in {time.perf_counter() - t0:.0f} s")
    meta = json.loads((OUT / "oc_best_found.json").read_text())
    meta["provenance"]["source_cells_sha256"] = src_hash
    # cross-check against notebook 04's saved, printed outputs
    pj = {}
    for m in re.finditer(r"^--- (S1|S2c) (\d+)s (B1|B2)-protocol (OC(?:-ISC5)?):.*?\n(.*?)(?=^--- |\Z)", printed, re.M | re.S):
        body = m.group(5)
        jd = re.search(r"best-found J_d = ([\d.]+) kWh_e", body)
        lb = re.search(r"J_LB = ([\d.]+) kWh_e", body)
        pj[f"{m[1]}_{m[2]}_{m[3]}_{m[4]}"] = (jd and jd[1], lb and lb[1])
    summ = printed[printed.rfind("STAGE 4 (perfect-foresight benchmark) VERIFICATION SUMMARY"):]
    pe = {f"{m[1]}_{m[2]}_{m[3]}_{m[4]}": m[5] for m in
          re.finditer(r"^(S1|S2c) (\d+)s (B1|B2)-prot (OC(?:-ISC5)?)\s+([\d.]+)", summ, re.M)}
    bad = []
    for tag, c in meta["cases"].items():
        now = (f"{c['J_d_best_kWh_e']:.5f}", f"{c['J_LB_certified_kWh_e']:.5f}", f"{c['replay_E_cool_5s_kWh_e']:.1f}")
        ref = pj.get(tag, (None, None)) + (pe.get(tag),)
        ok = now == ref
        print(f"{'OK  ' if ok else 'DIFF'} {tag:22s} exported {now} | notebook 04 printed {ref}")
        if not ok:
            bad.append(tag)
    meta["provenance"]["crosscheck_vs_notebook04_printed"] = "all 16 cases identical to printed precision" if not bad else f"MISMATCH: {bad}"
    (OUT / "oc_best_found.json").write_text(json.dumps(meta, indent=1))
    if bad or len(meta["cases"]) != 16:
        sys.exit(f"ABORT: exported artifacts do not reproduce notebook 04 ({bad})")
    print("all 16 cases reproduce notebook 04's printed values; artifacts written to artifacts/")


if __name__ == "__main__":
    main()
