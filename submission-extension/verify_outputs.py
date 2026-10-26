#!/usr/bin/env python3
"""Fail-closed substantive validation of extension results, not empirical validation."""
import csv, json, math
from pathlib import Path

def jread(path): return json.loads(Path(path).read_text(encoding='utf8'))
def csvread(path):
    with open(path,newline='',encoding='utf8') as f:return list(csv.DictReader(f))
def near(actual,expected,tol=0.02):
    assert math.isfinite(actual) and abs(actual-expected)<=tol,(actual,expected)

def verify(root):
    root=Path(root)
    fig=jread(root/'figure5/outputs/figure5_authority_validation.json')
    s14=jread(root/'s14/outputs/sensitivity_audit.json')
    s16=jread(root/'s16/complexity_cost_validation.json')
    s17=jread(root/'s17/multifactor_audit.json')
    a14=csvread(root/'s14/outputs/one_factor_sensitivity.csv')
    a16=csvread(root/'s16/complexity_cost_break_even.csv')
    a17=csvread(root/'s17/multifactor_sensitivity_720.csv')
    summary=csvread(root/'s17/multifactor_summary.csv')
    assert fig['status']=='PASS' and fig['checks']==6004
    near(fig['Bcrit_eur'],8419002.387956878,1e-4)
    near(fig['Ucrit_eur'],9190903.820731005,1e-4)
    assert s14['status']=='PASS' and s14['case_count']==14 and s14['independent_lp_cases']==300
    assert len(a14)==14
    near(float(a14[0]['Hmax']),6638497.612043122,1e-5)
    assert s16['status']=='PASS' and s16['checks']==2454 and len(a16)==9
    pp={float(x['lambda_h']):float(x['break_even_incremental_donor_cost_eur']) for x in a16}
    near(pp[.75],17833.930380366743,1e-4)
    near(pp[1.25],494002.38795687817,1e-4)
    assert s17['status']=='PASS' and s17['grid_cases']==720 and s17['lp_checks']==720 and s17['ramp_numeric_checks']==240
    assert len(a17)==720 and len(summary)==4
    assert s17['full_only_meets']==11 and s17['ramp_meets']==162
    assert sum(int(x['full_only']) for x in summary)==11
    full_only=[x for x in a17 if x['full_only_meets_need']=='True']
    assert len(full_only)==11 and all(x['alpha']=='0.25' for x in full_only)
    assert Path(root/'figure5/outputs/Figure_5_Fixed_Cash_Financing_Regimes.png').stat().st_size>10000
    return {'status':'PASS','checks':{'figure5_LP_and_thresholds':fig['checks'],'s14_LP':s14['independent_lp_cases'],'s16':s16['checks'],'s17_LP':s17['lp_checks'],'s17_ramp':s17['ramp_numeric_checks']},'results':{'figure5_budget_eur':fig['Bcrit_eur'],'figure5_exposure_eur':fig['Ucrit_eur'],'s16_premium_075':pp[.75],'s16_premium_125':pp[1.25],'s17_full_only':len(full_only)},'limit':'Conditional scenario and numerical code verification; cannot validate actual climate finance parameters or exhaustively prove a continuous-contract theorem'}

if __name__=='__main__':
    import sys
    print(json.dumps(verify(Path(__file__).resolve().parent),indent=2))
