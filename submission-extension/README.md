# Submission extension — pending repository integration

Reproduction code for new manuscript Figure 5 and Supplementary Sections S14–S17. This directory intentionally leaves the public historical code in `revision-2026-10-07/`, `revision-2026-10-12/` and `interval-extension/` **unchanged**.

From the repository root, run:

```sh
python submission-extension/run_all.py
```

The four module scripts generate results within their own directories; logs and the aggregate validation are written to `submission-extension/outputs/`. All scripts use fixed illustrative design parameters and fixed seeds. The requirement is Python 3.12+ with NumPy, SciPy and Matplotlib. A single full-root runner still requires an explicit one-line `run_all.py` integration as described in `../integration_instructions.md` in the staging package.

The Figure 5 `€8,419,002.39` and `€9,190,903.82` values are *conditional expenditure authorities* under the selected model inputs, not actual UEFA contractual figures. Table S16 values are break-even **additional sponsor costs** compared with optimized ramps, not measured costs. The 720-case S17 grid is illustrative, not a statistical sample and its counts are not probabilities. LPs audit derived finite-dimensional relationships; they do not independently prove all infinite-dimensional continuous-contract results.

The public repo should not be edited until the manuscript title and final artifact set are frozen; at that point update the root README, map and manifest and confirm the entire root runner on the supported locked CPython 3.12 environment.
