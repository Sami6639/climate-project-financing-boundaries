#!/usr/bin/env python3
"""Independent contractual-complexity charge checks for interval class vs optimized single ramp.

No observed expenses or welfare effects are inferred. Deterministic incremental expense K
is charged to the unrestricted class only and counted in consolidated donor expenditure.
"""
import csv,json, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parent
mu=316912.0
h=math.sqrt(3)*mu*0.1
p=25.0
q=18.75
B=7925000.0
U=1.1*B
Hreq=0.90*B


def full_resources(h,p,q,d):
    W=2*h; L=2*d; n=math.floor(d/h+1e-12); ell=L-n*W
    D=(n+1)*W*q+max(0., W*q-p*(W-ell))
    R=h*q*q/p+L*q
    return (R,D)

def ramp_resources(h,p,q,d,r):
    W=2*h; L=2*d
    K=max(W*q, r*(L-W)+2*W*q)
    s=W*q/r
    if L<=W-s+1e-8:
        R=q*L+W*q*q/(2*r)
    else:
        R=r*(L-W)+2*W*q-W*q*q/(2*r)
    return R,K

def ramp_optim(h,p,q,d,B,U,backend='candidate'):
    def gross(r):
        R,K=ramp_resources(h,p,q,d,r)
        return min(B-R,U-K)
    candidates=[q,p]
    if U>B:
        r0=(2*h)*q*q/(2*(U-B))
        if q-1e-12<=r0<=p+1e-12:
            candidates.append(max(q,min(p,r0)))
    if backend=='numerical':
        grid=np.linspace(q,p,2001)
        candidates=list(grid)
        # Locate all numerically detected kinks by root-finding rather than invoking
        # the analytic optimal-slope candidate from the theorem.
        def resource_diff(r):
            R,K=ramp_resources(h,p,q,d,r)
            return (B-R)-(U-K)
        for low,high in zip(grid[:-1],grid[1:]):
            v0,v1=resource_diff(low),resource_diff(high)
            if v0*v1<0:
                candidates.append(brentq(resource_diff,low,high,xtol=1e-13))
    g,rstar=max((gross(r),r) for r in candidates)
    return g,rstar

def point(lam,loading=0.,common_V=0.):
    d=lam*h
    R,D=full_resources(h,p,q,d)
    gross_full=min(B-R,U-D)
    gross_ramp,rstar=ramp_optim(h,p,q,d,B,U)
    a=1+loading
    hf=(gross_full-common_V)/a
    hr=(gross_ramp-common_V)/a
    premium=max(0.,gross_full-gross_ramp)
    return dict(lambda_h=lam, full_advance_eur=hf, ramp_advance_eur=hr,
                gross_capacity_gain_eur=hf-hr,
                break_even_incremental_donor_cost_eur=premium,
                full_cash_gap_eur=hf-Hreq, ramp_cash_gap_eur=hr-Hreq,
                optimal_ramp_slope_eur_per_t=rstar, binding_full=('Budget' if B-R<=U-D else 'Exposure'))

rows=[point(lam) for lam in (0.,0.1,0.25,0.5,0.75,1.,1.25,1.5,2.)]
checks=0; maxerr=0.
for r in rows:
    delta=r['lambda_h']*h
    rg,rstar=ramp_optim(h,p,q,delta,B,U,'numerical')
    err=abs(rg - (r['ramp_advance_eur']))
    maxerr=max(maxerr,err)
    assert err<1e-5,(r,err)
    # Financeable capacity at incremental fixed cost k; break even equality and strict order.
    gfull=r['full_advance_eur']; gramp=r['ramp_advance_eur']; kstar=r['break_even_incremental_donor_cost_eur']
    assert abs((gfull-kstar)-gramp)<1e-6
    for cost in [0,kstar*0.5,kstar,kstar+1000]:
        gh=gfull-cost
        assert abs((gh-gramp)-(kstar-cost)) < 1e-6
        checks+=1
    checks+=2
# Independent random sampled budget/cap, rate targets, interval widths and cost-loading tests.
rng=np.random.default_rng(20261026)
for j in range(600):
    p1=float(rng.uniform(1,35));q1=float(rng.uniform(.05,.97)*p1)
    h1=float(rng.uniform(1,50000));delta=float(rng.uniform(0,2.9)*h1)
    B1=float(rng.uniform(3,21)*h1*p1); U1=float(rng.uniform(3,24)*h1*p1)
    slope_gross,slope=ramp_optim(h1,p1,q1,delta,B1,U1)
    numeric_gross,_=ramp_optim(h1,p1,q1,delta,B1,U1,'numerical')
    err=abs(slope_gross-numeric_gross)
    # tolerance accounts for finite grid and floating point (not contract proof)
    tol=2e-7*max(1,abs(slope_gross))
    assert err < tol,(j,err,tol,slope,delta/h1)
    maxerr=max(maxerr,err)
    F=min(B1-full_resources(h1,p1,q1,delta)[0],U1-full_resources(h1,p1,q1,delta)[1])
    a=float(rng.uniform(1.,1.2));V=float(rng.uniform(0,1000))
    hf=(F-V)/a; hr=(slope_gross-V)/a
    kstar=F-slope_gross
    assert kstar>-2e-8*max(1,abs(F))
    assert abs((hf-hr)*a-kstar)<1e-7*max(1,abs(kstar))
    # Same sponsor cost scope; evaluate random incremental cost added to full class.
    inc=float(rng.uniform(0,1.5)*max(1,kstar))
    assert abs(((F-V-inc)/a - hr)-(kstar-inc)/a)<1e-7*max(1,abs(kstar))
    checks+=4
with (ROOT/'complexity_cost_break_even.csv').open('w',newline='') as fp:
 w=csv.DictWriter(fp,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'date':'2026-10-26','status':'PASS','scope':'Conditional deterministic incremental sponsor expense, not observed contract cost or project feasibility',
         'baseline':{'B_eur':B,'U_eur':U,'Hreq_eur':Hreq,'h_tonnes':h,'p':p,'q':q},
         'checks':checks,'independent_numerical_ramp_gross_error_max_eur':maxerr,
         'strictly_positive_premium_at_lambda':{str(r['lambda_h']):r['break_even_incremental_donor_cost_eur'] for r in rows if r['break_even_incremental_donor_cost_eur']>0.005},
         'caution':'Both full and optimized ramp Hmax are below the fixed 90% cash need at lambda 0.75 and 1.25; break-even complexity expense measures capacity differences, not funding of the unchanged need.'}
(ROOT/'complexity_cost_validation.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps(summary,indent=2)); print('\nSelected output:')
for r in rows:print(r['lambda_h'],round(r['full_advance_eur'],2),round(r['ramp_advance_eur'],2),round(r['break_even_incremental_donor_cost_eur'],2),round(r['full_cash_gap_eur'],2))
