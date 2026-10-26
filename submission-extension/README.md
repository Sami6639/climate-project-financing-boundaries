# Submission extension — integrated release

Reproduction code for new manuscript Figure 5 and Supplementary Sections S14–S17. This directory intentionally leaves the public historical code in `revision-2026-10-07/`, `revision-2026-10-12/` and `interval-extension/` **unchanged**.

From the repository root, run:

```sh
python submission-extension/run_all.py
```

The four module scripts generate results within their own directories; logs and the aggregate validation are written to `submission-extension/outputs/`. All scripts use fixed illustrative design parameters and fixed seeds. The requirement is Python 3.12+ with NumPy, SciPy and Matplotlib. The root `run_all.py` invokes these four analyses as its final reproduction step.

The Figure 5 `€8,419,002.39` and `€9,190,903.82` values are *conditional expenditure authorities* under the selected model inputs, not actual UEFA contractual figures. Table S16 values are break-even **additional sponsor costs** compared with optimized ramps, not measured costs. The 720-case S17 grid is illustrative, not a statistical sample and its counts are not probabilities. LPs audit derived finite-dimensional relationships; they do not independently prove all infinite-dimensional continuous-contract results.

The original historical results remain unchanged. The release manifest separately includes the extension sources, while the scientific CSV/JSON/PNG outputs are regenerated on execution. The full-root integration is checked in GitHub Actions on the locked CPython 3.12 environment.
