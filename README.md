# Dynamic Thermal Modelling and Energy-Efficient Cooling of an AI Data Hall

Mathematical Modelling (MAC F314), BITS Pilani K K Birla Goa Campus, 2026-27.

> **STATUS: REVIEW PACKAGE. NOT READY TO SUBMIT.** See [Outstanding items](#outstanding-items-before-submission) below. The compiled report PDF is available; the AI acknowledgement and slides are still missing.

## What to open first

1. **`notebook/main.ipynb`**: the single main notebook. It runs: problem → equations → fixed PI → adaptive PI → energy/safety comparison → compact OC benchmark → sensitivity → conclusions. It is stored already executed. Its last cell prints the check summary.
2. **`report/main.pdf`**: the compiled, 25-page report. Editable source: **`report/main.tex`**, with `report/refs.bib` and `report/figures/`.
3. **`AI_USE_LOG.md`**: a dated log of every AI use, for writing the acknowledgement.

## Team collaboration

Clone the project:

```bash
git clone https://github.com/ByteB1itz420/MM_proj.git
cd MM_proj
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the branch, commit, push and pull-request workflow. The repository owner must grant collaborator access before teammates can push. Private repositories also require access before cloning.

## How to run

The code was tested with Python 3.14 and the pinned versions in `requirements.txt`: numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2, casadi 3.8.2, nbconvert 7.17.2, nbformat 5.11.1 and ipykernel 7.4.0.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebook/main.ipynb    # about 3-4 min
```

- **Paths are relative.** Run notebooks from inside `notebook/`, as Jupyter does by default.
- **No external data.** The two synthetic scenarios (S1, S2c) are defined in code (`SCEN` in `main.ipynb`, `CMP` in notebooks 03–05).
- **To recompile the report** with Tectonic (built here with version 0.17.0):

  ```bash
  cd report && tectonic --keep-logs main.tex
  ```

  Install Tectonic on your own machine and ensure it is on `PATH`. On the original workstation it is installed at `~/.local/bin/tectonic`. Tectonic runs the bibliography and reference passes automatically. On Linux, packages and fonts are cached in `~/.cache/tectonic/`; new dependencies require internet on their first use. The output is `report/main.pdf`.

  With a traditional TeX distribution, the alternative is:

  ```bash
  cd report && pdflatex main && bibtex main && pdflatex main && pdflatex main
  ```

  Figures are referenced without an extension, so `pdflatex` uses `figures/*.pdf` and `latex`/`dvips` uses `figures/*.eps`. The appendix code listings are read from `../src/`.

## Loaded versus recomputed results in `main.ipynb`

| Result | Status |
|---|---|
| Model checks, PI gain design, all 28 PI runs × 2 grids (B1/B2 comparison and sensitivity), electricity accounting, tables, figures | **Recomputed** |
| OC benchmark: best-found command and state arrays | **Loaded** from `artifacts/oc_best_found.npz` (+ `.json` provenance) |
| OC feasibility, $J_d$, certified LP lower bounds, replays on the 5 s and 1 s grids, energy balance, terminal checks | **Recomputed** from the loaded arrays |

**Where the OC arrays come from.** They were produced by `scripts/export_oc_artifacts.py`.

- That script executes notebook 04's own solve and replay cells **unchanged** in a separate kernel, then saves the best-found arrays.
- It checks every exported value against the values printed in notebook 04: all 16 cases are identical to the printed precision.
- To re-optimise instead (about 15–20 min), set `RUN_FULL_OC = True` in `main.ipynb`, or run the script.

## Verification: where the full checks are

| Notebook | Content | Checks |
|---|---|---|
| `notebook/01_thermal_model_verification.ipynb` | model, exact ZOH vs Radau, energy balance, eigenvalues, COP domain, invalid inputs | 50/50 |
| `notebook/02_fixed_setpoint_pi.ipynb` | fixed-setpoint PI, gain design, anti-windup, intersample extrema | 36/36 |
| `notebook/03_load_adaptive_pi_energy.ipynb` | load-adaptive PI, electricity accounting, B1/B2 comparison | 52/52 |
| `notebook/04_perfect_foresight_benchmark.ipynb` | OC benchmark, LP certificate, replays, constant-COP invariance | 197/197 |
| `notebook/05_margin_setpoint_sensitivity.ipynb` | margin / setpoint sensitivity, strict 5 s + 1 s safety rule | 152/152 |
| `notebook/main.ipynb` | fast re-checks; reproduces notebooks 03 and 05; re-verifies the loaded OC solutions | 28/28 |

These checks are numerical **verification** (the code solves our equations correctly). They are **not** real-world validation.

## Scope of the claims

- **Inputs and parameters.** All parameters are illustrative. Only two synthetic step scenarios were tested; no real workload or weather data were used, and nothing was validated against a real facility.
- **Gains.** The PI gains are stability-tuned (minimum closed-loop spectral radius with a 2× gain margin), not energy-tuned.
- **Run protocol.** Every run starts at its own equilibrium and must return to it over the same 17 h horizon. Each OC problem uses the start/end state of the PI controller it is compared with.
- **When savings are claimed.** Only when both compared runs pass the terminal checks and stay within limits on both the 5 s and the 1 s grid.
- **S1 OC-ISC5.** It exceeds a limit on the 1 s grid, by up to 7.5e-4 K at 60 s and 4.8e-5 K at 300 s. These results are diagnostic only: not continuously safe, and not eligible.
- **Energy figures and the LP gaps.** Energy results use replay-integrated electricity, not the optimiser's discrete objective $J_d$. The LP gaps concern $J_d$ only.

## Package contents

```
README.md                  this file
CONTRIBUTING.md            team access, branch and pull-request workflow
requirements.txt           pinned dependencies
AI_USE_LOG.md              detailed AI-use log (facts unchanged; table formatting fixed)
PROJECT_SPEC.md            project specification and change log (supporting)
src/                       thermal_model.py, controllers.py, optimal_control.py
notebook/main.ipynb        MAIN NOTEBOOK
notebook/01-05_*.ipynb     supporting notebooks with the full verification
artifacts/                 oc_best_found.npz + oc_best_found.json (provenance)
scripts/                   export_oc_artifacts.py (regenerates the OC artifacts)
figures/main/              main-notebook figures (.png, .eps, .pdf); fig1 is the model diagram, drawn by code in main.ipynb
figures/                   supporting-notebook figures (fig01-fig14)
report/                    main.tex, refs.bib, figures/ (.eps, .pdf, .png incl. BITS logo from the template)
```

## Outstanding items before submission

1. **AI acknowledgement (team).** It is intentionally blank in the report (*Generative AI Usage Declaration*) and in `main.ipynb` (*AI acknowledgement*). Write it from `AI_USE_LOG.md`. The brief warns that dishonest acknowledgement can cost 10 marks.
   - Notebooks 01–05 each contain a short note, written earlier, saying that Claude drafted them. Check that it agrees with your declaration.
2. **Compiled report PDF — completed 30 Oct 2026.** `report/main.pdf` was compiled with Tectonic 0.17.0. All 25 pages were rendered and visually reviewed; title-page overflow, wide tables, long equations and appendix headings were corrected. The abstract was subsequently rewritten from the model, notebook outputs and conclusions, and the PDF rebuilt and checked. The final build has no overfull boxes, undefined citations or unresolved references. Recheck after any content changes.
3. **Report fields left blank on purpose** (never invented):
   - Date of submission;
   - Acknowledgments (personal);
   - signatures and dates.
4. **Sources.**
   - **Check every source.** Bibliographic details were checked against DOI or catalogue records, but a team member should open each source and confirm the statement it is cited for.
   - **18–27 °C air band.** It is labelled illustrative and **not cited**: the official ASHRAE document could not be accessed.
   - **Carnot COP formula.** Add the thermodynamics textbook the team uses.
5. **Appendix A** shows only core code excerpts. The template asks for "all relevant programming code", and the full code is in `src/` and `main.ipynb`. Decide whether to include more.
6. **Slides.** Not yet made: 10 min talk plus 5 min Q&A, and every member must explain their own part.
7. **Deadline.** The brief says the `.tex`, PDF and `.ipynb` must be submitted **at least one week before** the presentation. Confirm the agreed date with the instructor.
8. **Bonus.** No bonus marks are claimed or guaranteed. Under the brief, bonus marks depend on mathematical depth and credible validation, and this project has no real-data validation.

## Checked against the course brief

| Brief requirement | Status |
|---|---|
| Editable `.tex` | done (`report/main.tex`) |
| Compiled PDF | done (`report/main.pdf`; item 2) |
| References | done (`report/refs.bib`); team verification pending (item 4) |
| Figure files in `.eps` | done (`report/figures/*.eps`, `figures/main/*.eps`) |
| Report content: introduction, motivation, objectives, assumptions, formulation, diagrams, analysis (equilibria, stability, sensitivity), simulations, graphs, interpretation | done (Chapters 1–3) |
| Executable `.ipynb`, colourful Markdown, equations and algorithms before code, labelled plots, interpretations | done (`notebook/main.ipynb`) |
| Acknowledge all AI platforms with their roles | **pending** (item 1) |
| Cite external datasets and results | no external datasets used; methods and software cited |
| Presentation | **pending** (item 6) |
