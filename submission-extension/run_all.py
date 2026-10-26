#!/usr/bin/env python3
"""Reproduce post-2026-10-15 manuscript extensions without changing historical code.

This command is designed to be invoked from the repository root runner after
revision-2026-10-07, revision-2026-10-12 and interval-extension scripts finish.
"""
from __future__ import annotations
import json, os, subprocess, sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if sys.flags.optimize:
    raise RuntimeError('Never run scientific validation under Python -O')

MODULES = [
  ('figure5', 'figure5_and_authority_regions.py'),
  ('s14', 'table_s14_sensitivity.py'),
  ('s16', 'reproduce_complexity_premium.py'),
  ('s17', 'multifactor_robustness.py'),
]

def main():
    out = ROOT / 'outputs'
    out.mkdir(exist_ok=True)
    for sub, script in MODULES:
        target = ROOT / sub / script
        env = os.environ.copy()
        env['PYTHONHASHSEED'] = '0'
        done = subprocess.run([sys.executable, str(target)], cwd=target.parent,
                              env=env, text=True, capture_output=True)
        log = (done.stdout or '') + ('\nSTDERR:\n' + done.stderr if done.stderr else '')
        (out / f'{sub}_run.log').write_text(log,encoding='utf8')
        if done.returncode:
            print(f'FAIL: {sub}: {done.returncode}\n{log[-4000:]}',file=sys.stderr)
            raise SystemExit(done.returncode)
        print(f'PASS {sub}: {script}')
    from verify_outputs import verify
    report = verify(ROOT)
    (out / 'submission_extension_verification.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print('PASS submission extension:',report['checks'])

if __name__=='__main__': main()
