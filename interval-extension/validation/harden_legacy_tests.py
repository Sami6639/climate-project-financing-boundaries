"""Retain every original test and add fail-hard numerical gates in stated units."""
from pathlib import Path
import importlib.util,json,sys
import numpy as np
OUT=Path(__file__).resolve().parent
LEGACY=OUT/'legacy'/'Interval_Financing_Reproducibility'
if not LEGACY.exists():LEGACY=OUT.parent # integrated package fallback
sys.path.insert(0,str(LEGACY/'validation'))
import stress_formula

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def run():
    stress_formula.run()
    result=json.loads((stress_formula.OUT/'stress_results.json').read_text());summary=result['summary']
    assert summary['construction_cases']==112 and summary['lp_cases']==44
    tolerances={'max_exact_error':2e-10,'max_bracket_violation':2e-8,'max_primal_violation':2e-8,'max_duality_gap':2e-8}
    for key,tol in tolerances.items():assert summary[key]<=tol,(key,summary[key],tol)
    cost=load(LEGACY/'cost'/'cost_comparison.py','legacy_cost_validation')
    rng=np.random.default_rng(20261014);max_abs=0.;max_rel=0.;max_cash=0.;max_gate=0.
    for _ in range(1000):
        H=float(rng.uniform(1,1e7));a=float(rng.uniform(1,1.4));V=float(rng.uniform(0,1e6));R=float(rng.uniform(0,2e6));D=float(rng.uniform(0,4e6));B=float(rng.uniform(0,1.5e7));U=float(rng.uniform(0,1.5e7));cb,cu,cs=rng.uniform(.01,2,3)
        row=cost.solve(H,B,U,R,D,a,V,float(cb),float(cu),float(cs))
        err=abs(row['weighted_resource_cost']-row['lp_cost']);max_abs=max(max_abs,err);max_rel=max(max_rel,err/max(1,abs(row['lp_cost'])))
        cash=abs(row['donor_advance']+row['external_nonrepayable_funds']-H);max_cash=max(max_cash,cash)
        gate=max(0,row['budget_deficit']-row['budget_upgrade']-a*row['external_nonrepayable_funds'],row['exposure_deficit']-row['exposure_upgrade']-a*row['external_nonrepayable_funds']);max_gate=max(max_gate,gate)
        assert err<=1e-6 and cash<=1e-6 and gate<=1e-6,(err,cash,gate)
    report=dict(original_test_counts=dict(construction=112,LP=44,cost=1000),seed=20261014,uniform_summary=summary,uniform_thresholds=tolerances,uniform_units='normalized money p*h; normalized slope p for marginal residuals',cost_summary=dict(max_absolute_objective_error_eur=max_abs,max_relative_objective_error=max_rel,max_initial_cash_error_eur=max_cash,max_financing_gate_violation_eur=max_gate,absolute_tolerance_eur=1e-6),all_passed=True)
    (OUT/'legacy_hard_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report
if __name__=='__main__':run()
