# Financing Climate Projects: The Boundaries of Financial Capacity

Data and original scientific code supporting the supplied manuscript and supplementary information. This package combines the preserved October 7 local baseline, October 12 independent comparisons, and the matching October 14 strengthened interval/financing extension. It does not contain the later Ghana case.

## Reproduce

Tested runtime: CPython 3.12.14 on Linux x86_64. Create a virtual environment and install the complete hash-locked dependencies:

    python3.12 -m venv .venv
    .venv/bin/python -m pip install --only-binary=:all: --require-hashes -r interval-extension/requirements-lock-linux-cp312.txt
    .venv/bin/python verify_manifest.py
    .venv/bin/python run_all.py

No network is used for scientific execution after installation. Do not use Python -O. The entry point runs each original scientific implementation without changing its mathematical model. Outputs and execution logs are regenerated in package-relative folders. `interval-extension/run_all.py --verify-only` checks existing interval outputs. Hash verification applies to the initial snapshot; regeneration changes runtime records and PDF metadata.

## What is data and what is assumed

UEFA documents supply the reported inventory, announced fund scale, reported investment and programme rules. Recipient reports supply selected project amounts and advance narratives. These are documentary anchors, not an estimated dataset of project cash requirements or emissions-reporting errors. Reporting distributions, dispersion, target rates, aggregate cash need, donor exposure limits, financing costs and resource prices are hypothetical scenarios. No causal investment or emissions estimate is made.

The institutional source register contains URLs and locators. Third-party full texts, manuscript files, editorial working material and private notes are excluded. Existing scientific outputs include CSV/JSON tables and generated figures. Historical baseline filenames retain their original numbering; use `exhibit_file_map.csv` for the current manuscript numbering. In particular, current Figure 3 is `interval-extension/theory/Figure_3_Three_Way_Frontier.png`, not the old baseline figure3.png.

## Verification

The matching original extension archive passed its 109-file SHA-256 manifest. A fresh isolated-directory rerun on October 15 reproduced all 38 non-PDF extension outputs byte-for-byte and passed all 41 output gates. Three generated PDFs differ only in document metadata; plotted output is unchanged. The earlier baseline passed its 27-file manifest and all 17 scientific output comparisons. The supplied manuscript's four embedded figures match their authoritative generated PNGs exactly. Main numerical tables and supplementary numerical comparisons were checked against the supplied Word files. See `verification_summary.json` and `exhibit_file_map.csv` for scope and boundaries.

Original software license scope is specified in LICENSE and THIRD_PARTY_DATA.md. Source documents remain with their owners. The historical baseline repository is https://github.com/Sami6639/uefa-euro2024-carbon-linked-finance at commit 05d5f580bfe9a1b216fe900983c414abede0d484. The interval extension was a separately supplied companion archive; this package does not claim it was part of that historical commit. No archival DOI or formal GitHub release is asserted.
