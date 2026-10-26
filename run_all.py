"""Run the preserved local, independent, and interval/financing computations."""
from pathlib import Path
import subprocess,sys
R=Path(__file__).resolve().parent
if sys.flags.optimize: raise RuntimeError("Do not use Python -O")
for p in ["revision-2026-10-12/run_all.py", "revision-2026-10-12/math/independent_math_audit_v1.py", "interval-extension/run_all.py"]:
    subprocess.run([sys.executable,str(R/p)],cwd=R,check=True)
subprocess.run([sys.executable,str(R/'submission-extension/run_all.py')],cwd=R,check=True)
print('PASS: all historical and submission-extension scientific entry points completed')
