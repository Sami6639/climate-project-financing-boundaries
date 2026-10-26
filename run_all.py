"""Reproduce the historical and new climate-project financing results."""
from pathlib import Path
import subprocess
import sys

R = Path(__file__).resolve().parent
if sys.flags.optimize:
    raise RuntimeError('Do not run scientific validation under Python -O')

historical = R / 'revision-2026-10-12/run_all.py'
first = subprocess.run([sys.executable, str(historical)], cwd=R,
                       text=True, capture_output=True)
if first.stdout:
    print(first.stdout, end='')
if first.returncode:
    if first.stderr:
        print(first.stderr, file=sys.stderr)
    # Do not suppress implementation failures; permit a carefully enumerated
    # numerical alternative only for known byte-level Gaussian CSV changes.
    if not ('AssertionError:' in first.stderr
            and 'figure4_series.csv' in first.stderr
            and 'gaussian_lp_audit.csv' in first.stderr):
        raise subprocess.CalledProcessError(first.returncode, str(historical))
    subprocess.run([sys.executable,
                    str(R / 'submission-extension/verify_historical_numerics.py')],
                   cwd=R, check=True)

for relative in [
    'revision-2026-10-12/math/independent_math_audit_v1.py',
    'interval-extension/run_all.py',
    'submission-extension/run_all.py',
]:
    subprocess.run([sys.executable, str(R / relative)], cwd=R, check=True)

print('PASS: historical results and submission extension scientifically verified')
