# Working together

Start with `README.md`, `notebook/main.ipynb` and `report/main.pdf`. The report source is `report/main.tex`; the model and controllers are in `src/`.

## Clone and make a change

```bash
git clone https://github.com/ByteB1itz420/MM_proj.git
cd MM_proj
git switch main
git pull --ff-only
git switch -c your-name/short-description
```

Follow the environment and build instructions in the README. Make a focused change, review it, then commit and push your branch:

```bash
git status
git diff
git add path/to/changed-file
git commit -m "Describe what changed"
git push -u origin your-name/short-description
```

Open a pull request on GitHub from your branch into `main`, describe the change and how you checked it, and ask another team member to review it. Coordinate notebook edits before starting, since simultaneous edits to the same notebook can be difficult to merge. Pull the latest `main` before beginning the next change.

## Check the parts you change

- For report edits, rebuild `report/main.pdf` and inspect the changed pages, citations and cross-references. Commit the `.tex` and updated PDF together.
- For model or controller changes, run the relevant supporting notebook checks and the main notebook. Update affected figures, tables and report claims together. Preserve executed notebook outputs needed to review the results.
- If you regenerate the optimal-control solutions, use `scripts/export_oc_artifacts.py` and commit the `.npz` and `.json` provenance together.
- Record AI assistance in `AI_USE_LOG.md`. Keep claims scoped to the synthetic scenarios and distinguish numerical verification from real-facility validation.

The repository tracks source, executed notebooks, numerical artifacts, figures, reference documents and the compiled report. Local environments, caches, build logs and generated ZIP archives are excluded; teammates recreate the Python environment from `requirements.txt`.

## Access

Direct pushes require collaborator access. The repository owner can invite each teammate through **Settings → Collaborators**. Accept the GitHub invitation before cloning a private repository or pushing a branch. Authenticate using GitHub CLI (`gh auth login`), a credential manager or SSH; never commit access tokens or passwords.
