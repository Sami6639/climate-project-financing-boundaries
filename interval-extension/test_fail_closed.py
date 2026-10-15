"""Deliberately corrupt temporary report copies and prove entry-command failure.

The package's own scientific outputs are never changed by these negative tests.
Run after python run_all.py. Evidence is written under evidence/.
"""
from pathlib import Path
import json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent

def mutate(path,fn):
    d=json.loads(path.read_text());fn(d);path.write_text(json.dumps(d,indent=2))

def run():
    cases=[
      ('construction_residual','validation/stress_results.json',lambda d:d['exact_construction'][0].update(mass_error=1e-4)),
      ('missing_construction_case','validation/stress_results.json',lambda d:d['exact_construction'].pop()),
      ('LP_bracket_residual','validation/stress_results.json',lambda d:d['lp_cases'][0].update(theory_above_upper=1e-4)),
      ('infeasibility_failure','validation/inadmissible_cases.json',lambda d:d[0].update(actual_infeasible=False)),
      ('nonfinite_number','validation/stress_results.json',lambda d:d['exact_construction'][0].update(reserve_error=float('nan'))),
      ('cost_objective_mismatch','cost/cost_results.json',lambda d:d['scenario_rows'][0].update(lp_cost=-100.)),
      ('single_ramp_formula_failure','theory/single_ramp_self_check.json',lambda d:d.update(maximum_scaled_order_error=1e-3)),
      ('independent_ramp_integral_failure','validation/optimized_ramp_checks.json',lambda d:d['summary'].update(max_independent_integral_error=1e-3)),
      ('Gaussian_continuum_failure','validation/gaussian_interval_certificates.json',lambda d:d['records'][0].update(gaussian_continuum_marginal_lower_over_p=-1)),
      ('price_region_face_failure','cost/price_region_validation_report.json',lambda d:d.update(max_normalized_optimal_face_error=1e-2)),
      ('primary_figure_failure','validation/ramp_equality_and_figure_audit.json',lambda d:d.update(maximum_figure_error_eur=1.)),
    ]
    results=[]
    with tempfile.TemporaryDirectory(prefix='interval-failure-gates-') as tmp:
        p=Path(tmp)/'package'
        shutil.copytree(ROOT,p,ignore=shutil.ignore_patterns('.venv','__pycache__','*.pyc','_runtime'))
        command=[sys.executable,str(p/'run_all.py'),'--verify-only']
        positive=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        if positive.returncode!=0:
            raise RuntimeError('Unmodified positive control failed: '+positive.stdout)
        for name,path,fn in cases:
            f=p/path;original=f.read_bytes()
            mutate(f,fn)
            result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            f.write_bytes(original)
            if result.returncode==0 or 'FAIL:' not in result.stdout:
                raise RuntimeError(f'Entry command did not reject {name}')
            results.append({'case':name,'returncode':result.returncode,'message':result.stdout.strip(),'passed':True})
        optimized=subprocess.run([sys.executable,'-O',str(p/'run_all.py'),'--verify-only'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        if optimized.returncode==0:raise RuntimeError('Optimized Python was accepted')
        results.append({'case':'python_assertion_optimization','returncode':optimized.returncode,'message':optimized.stdout.strip(),'passed':True})
    report={'status':'PASS','positive_control_returncode':positive.returncode,'negative_cases':len(results),'results':results}
    (ROOT/'evidence').mkdir(exist_ok=True)
    (ROOT/'evidence/failure_gate_self_test.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report
if __name__=='__main__':run()
