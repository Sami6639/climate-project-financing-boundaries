"""Same-project financing repair under a budget ceiling; no borrowing.
All prices are hypothetical opportunity costs per euro of authority/resources.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

OUT=Path(__file__).parent

def solve(Hreq,B,U,R,D,a=1.,V=0.,cb=.1,cu=.05,cs=.4):
    assert Hreq>=0 and a>0 and min(cb,cu,cs)>0
    A=a*Hreq+V+R-B
    E=a*Hreq+V+D-U
    candidates=sorted(set([0.,Hreq,max(0.,min(Hreq,A/a)),max(0.,min(Hreq,E/a))]))
    lp=linprog([cb,cu,cs],A_ub=[[-1.,0.,-a],[0.,-1.,-a]],b_ub=[-A,-E],bounds=[(0,None),(0,None),(0,Hreq)],method='highs')
    # Objective includes the external resource cost; explicit form for audit.
    rows=[]
    for s0 in candidates:
        b0=max(0.,A-a*s0);u0=max(0.,E-a*s0)
        rows.append((cb*b0+cu*u0+cs*s0,b0,u0,s0))
    J,b,u,s=min(rows)
    assert lp.success,lp.message
    assert np.isclose(J,lp.fun,rtol=1e-9,atol=1e-6),(J,lp.fun)
    assert b+a*s>=A-1e-6 and u+a*s>=E-1e-6
    return dict(Hreq=Hreq,B=B,U=U,R=R,D=D,a=a,V=V,cb=cb,cu=cu,cs=cs,budget_deficit=A,exposure_deficit=E,budget_upgrade=b,exposure_upgrade=u,external_nonrepayable_funds=s,donor_advance=Hreq-s,total_initial_project_cash=Hreq,weighted_resource_cost=J,lp_cost=float(lp.fun),lp_solution=lp.x.tolist(),candidate_s=candidates)

def validate_random():
    rng=np.random.default_rng(20261014)
    for _ in range(1000):
        H=float(rng.uniform(1,1e7));a=float(rng.uniform(1,1.4));V=float(rng.uniform(0,1e6));R=float(rng.uniform(0,2e6));D=float(rng.uniform(0,4e6));B=float(rng.uniform(0,1.5e7));U=float(rng.uniform(0,1.5e7));cb,cu,cs=rng.uniform(.01,2,3)
        solve(H,B,U,R,D,a,V,float(cb),float(cu),float(cs))
    return 1000

if __name__=='__main__':
    n=validate_random()
    cfgfile=OUT/'interval_inputs.json'
    rows=[]
    if cfgfile.exists():
        cfg=json.loads(cfgfile.read_text())
        for label,cb,cu,cs in [('Authority expansion',.25,.125,1.),('Mixed financing',.25,1.5,1.),('External resources',1.5,.25,1.)]:
            row=solve(**cfg,cb=cb,cu=cu,cs=cs);row['scenario']=label;rows.append(row)
    result={'random_positive_price_cases_passed':n,'scenario_rows':rows,'scope':'Budget ceiling. Known feasible reserve/span representation. External funds nonrepayable and available at initial cash date. Unit costs hypothetical.'}
    (OUT/'cost_results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
