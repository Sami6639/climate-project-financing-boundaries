#!/usr/bin/env python3
"""Verify the October 2026 historical rerun with narrow numerical tolerances.

Historical original files and original rerun script are never modified. This
fallback is used only if the original rerun fails its byte-identical CSV check.
Non-Gaussian CSV/JSON outputs must still agree exactly and PNG pixels must match.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
from PIL import Image, ImageChops

ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / 'revision-2026-10-07' / 'results'
REPLAY = ROOT / 'revision-2026-10-12' / 'portable_rerun' / 'results'
LOGS = ROOT / 'revision-2026-10-12' / 'logs'
ALLOWED_NUMERIC = {'figure4_series.csv', 'gaussian_lp_audit.csv'}
ABS_TOL = 1e-10
REL_TOL = 1e-10


def check_csv(a: Path, b: Path) -> dict:
    with a.open(newline='', encoding='utf8') as f:
        rdr = csv.DictReader(f)
        fields_a, rows_a = rdr.fieldnames, list(rdr)
    with b.open(newline='', encoding='utf8') as f:
        rdr = csv.DictReader(f)
        fields_b, rows_b = rdr.fieldnames, list(rdr)
    assert fields_a == fields_b and len(rows_a) == len(rows_b), (a.name, 'CSV grain/columns')
    maximum_absolute = 0.0
    maximum_relative = 0.0
    n = 0
    for i, (old, new) in enumerate(zip(rows_a, rows_b)):
        for col in fields_a:
            x, y = old[col], new[col]
            if x == y:
                continue
            try:
                vx, vy = float(x), float(y)
            except (TypeError, ValueError):
                raise AssertionError((a.name, i, col, x, y))
            assert math.isfinite(vx) and math.isfinite(vy), (a.name, i, col, x, y)
            diff = abs(vx - vy)
            maximum_absolute = max(maximum_absolute, diff)
            maximum_relative = max(maximum_relative, diff / max(1, abs(vx), abs(vy)))
            assert diff <= ABS_TOL + REL_TOL * max(abs(vx), abs(vy)), (a.name, i, col, vx, vy, diff)
            n += 1
    return {'rows': len(rows_a), 'differing_numeric_cells': n,
            'maximum_absolute_difference': maximum_absolute,
            'maximum_scaled_difference': maximum_relative}


def main() -> None:
    assert PUBLISHED.is_dir() and REPLAY.is_dir(), 'Missing original/replayed results'
    for name in ['portable_original.log', 'portable_independent.log']:
        assert (LOGS / name).is_file() and (LOGS / name).stat().st_size > 0, name
    original_files = {p.name for p in PUBLISHED.iterdir() if p.is_file() and p.name != 'runtime.json'}
    replay_files = {p.name for p in REPLAY.iterdir() if p.is_file() and p.name != 'runtime.json'}
    assert original_files == replay_files, (sorted(original_files - replay_files), sorted(replay_files - original_files))
    results = {}
    for name in sorted(original_files):
        a, b = PUBLISHED / name, REPLAY / name
        if name in ALLOWED_NUMERIC:
            results[name] = check_csv(a, b)
        elif a.suffix == '.csv':
            assert a.read_bytes() == b.read_bytes(), (name, 'non-Gaussian CSV differs')
        elif a.suffix == '.json':
            assert json.loads(a.read_text()) == json.loads(b.read_text()), (name, 'JSON differs')
        elif a.suffix == '.png':
            with Image.open(a) as ai, Image.open(b) as bi:
                x, y = ai.convert('RGB'), bi.convert('RGB')
                assert x.size == y.size and ImageChops.difference(x, y).getbbox() is None, name
        else:
            assert a.read_bytes() == b.read_bytes(), name
    report = {'status': 'PASS', 'historical_files_checked': len(original_files),
              'non_Gaussian_results_unchanged': len(original_files) - len(ALLOWED_NUMERIC),
              'Gaussian_numeric_outputs': results,
              'absolute_tolerance': ABS_TOL, 'relative_tolerance': REL_TOL,
              'note': 'A narrow numerical check replaces bitwise equality only for two listed Gaussian CSV files.'}
    (LOGS / 'gaussian_numeric_replay_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS: historical results reproduced; Gaussian CSV differences bounded:', results)


if __name__ == '__main__':
    main()
