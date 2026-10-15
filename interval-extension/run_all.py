"""Regenerate and audit every scientific result; exit nonzero on any failure.

Usage: python run_all.py
       python run_all.py --verify-only  (audit already-generated reports)
This is a scientific reproducibility entry point, not a software-product test.
"""
import argparse, importlib.metadata, json, os, subprocess, sys, tempfile
from pathlib import Path
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'validation'))
from audit_baseline_outputs import verify as verify_baseline, finite_tree
from audit_extension_outputs import verify as verify_extensions

def load_plan():
    return json.loads((ROOT / 'execution_plan.json').read_text())

def verify_outputs(plan):
    for name in plan['generated_outputs']:
        p = ROOT / name
        if not p.is_file() or p.stat().st_size == 0:
            raise RuntimeError(f'Missing/empty scientific output: {name}')
        if p.suffix == '.json':
            finite_tree(json.loads(p.read_text()), name)
    result = {'baseline':verify_baseline(ROOT),'extensions':verify_extensions(ROOT)}
    for item in plan.get('report_checks', []):
        data = json.loads((ROOT / item['path']).read_text())
        value = data
        for key in item['keys']:
            value = value[key]
        if item['operator'] == 'equals':
            passed = value == item['value']
        elif item['operator'] == 'maximum':
            passed = isinstance(value, (int,float)) and 0 <= value <= item['value']
        else:
            raise RuntimeError(f"Unknown check operator: {item['operator']}")
        if not passed:
            raise RuntimeError(f"Report check failed: {item['path']}:{item['keys']} = {value}")
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'):
        raise RuntimeError('Python optimization must be disabled so original assertion checks cannot be skipped.')
    plan = load_plan()
    (ROOT/'run_summary.json').write_text(json.dumps({'status':'RUNNING','mode':'verify-only' if args.verify_only else 'regenerate-and-verify'},indent=2)+'\n')
    if not args.verify_only:
        runtime = ROOT / '_runtime'; runtime.mkdir(exist_ok=True)
        env = os.environ.copy()
        for key in ('PYTHONPATH','PYTHONHOME','PYTHONOPTIMIZE'):
            env.pop(key,None)
        env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1')
        with tempfile.TemporaryDirectory(prefix='interval-mpl-') as mpl, tempfile.TemporaryDirectory(prefix='interval-cache-') as cache:
            env.update(MPLCONFIGDIR=mpl,XDG_CACHE_HOME=cache)
            with (runtime/'execution.log').open('w') as log:
                for script in plan['scripts']:
                    print('Running',script,flush=True)
                    log.write('\n=== '+script+' ===\n'); log.flush()
                    p = subprocess.run([sys.executable, '-B', str(ROOT/script)],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
                    log.flush()
                    if p.returncode:
                        raise RuntimeError(f'{script} failed with exit {p.returncode}. See _runtime/execution.log.')
    checks=verify_outputs(plan)
    versions={name:importlib.metadata.version(name) for name in plan['dependencies']}
    report={'status':'PASS','python':sys.version.split()[0],'dependency_versions':versions,'random_seeds':plan['random_seeds'],'tolerances':plan['tolerances'],'checks':checks,'completed_scripts':plan['scripts'],'generated_output_count':len(plan['generated_outputs'])}
    (ROOT/'run_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    print('PASS: all scientific outputs regenerated and audited.' if not args.verify_only else 'PASS: existing scientific outputs audited.')

if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        (ROOT/'run_summary.json').write_text(json.dumps({'status':'FAIL','error':str(exc)},indent=2)+'\n')
        print(f'FAIL: {exc}', file=sys.stderr)
        sys.exit(1)
