#!/usr/bin/env python3
"""Fixed illustrative full-factorial robustness audit for sharp interval vs optimized ramp.

This is a hypothetical stress grid, not a statistical sample of climate projects.
Independent levers are target fraction alpha=q/p, physical coverage delta as fraction
of fixed reference h0, reporting uncertainty as relative standard deviation sigma/mu,
and realized exposure limit U/B. Physical delta does not co-move when sigma changes.
"""
import csv, json, math, statistics
from pathlib import Path
from itertools import product
import numpy as np
from scipy.optimize import linprog, minimize_scalar, brentq
ROOT=Path(__file__).resolve().parent
mu=316912.;p=25.;B=7_925_000.;Hneed=.9*B
h0=math.sqrt(3)*.1*mu
ALPHAS=(.25,.5,.75,.9)
LAMBDAS=(0.,.1,.25,.5,.75,1.,1.25,1.5,2.) # delta/h0, not delta/h
SIGMAS=(.05,.075,.10,.125,.15)
CAPRATIOS=(1.05,1.1,1.2,1.3)
TOL=1e-7

def reqs(h,p,q,d):
    W=2*h;L=2*d
    # integer-period boundaries treated with a small relative tolerance
    ratio=L/W
    nearest=round(ratio)
    if abs(ratio-nearest)<1.e-12: ratio=float(nearest)
    n=math.floor(ratio)
    ell=L-n*W
    ell=max(0.,min(W,ell))
    D=(n+1)*W*q+max(0.,W*q-p*(W-ell))
    R=h*q*q/p+L*q
    return R,D

def rampR(h,p,q,d,r):
    W=2*h;L=2*d
    s=W*q/r
    rise=max(W*q,r*(L-W)+2*W*q)
    if L<=W-s+1e-8*max(1.,W):
        reserve=q*L+W*q*q/(2*r)
    else:
        reserve=r*(L-W)+2*W*q-W*q*q/(2*r)
    return reserve,rise

def analytic_ramp(h,p,q,d,B,U):
    # Endpoint/crossing result of Proposition R, checked independently below
    slopes=[q,p]
    if U>B:
        r0=h*q*q/(U-B)
        if q<=r0<=p and 2*d>=2*h-2*h*q/r0-1e-8:
            slopes.append(r0)
    return max((min(B-rampR(h,p,q,d,r)[0],U-rampR(h,p,q,d,r)[1]),r) for r in slopes)

def numerical_ramp(h,p,q,d,B,U):
    # Dense direct scan with golden-section searches on both sides of all kinks;
    # the numerical search does not call analytic candidate calculations.
    def f(r):
        rr,dd=rampR(h,p,q,d,r)
        return min(B-rr,U-dd)
    x=np.linspace(q,p,301)
    v=np.array([f(z) for z in x])
    idxs=[int(np.argmax(v)),0,len(x)-1]
    best=max(v)
    for i in set(idxs):
        lo=x[max(i-1,0)];hi=x[min(i+1,len(x)-1)]
        if hi>lo:
            s=minimize_scalar(lambda r:-f(r),bounds=(lo,hi),method='bounded',options={'xatol':1e-12,'maxiter':100})
            best=max(best,-s.fun)
    # exact crossing localized by independently bracketing the budget/cap equality
    def g(r):
        rr,dd=rampR(h,p,q,d,r)
        return (B-rr)-(U-dd)
    for a,b in zip(x[:-1],x[1:]):
        if g(a)*g(b)<0:
            rr=brentq(g,a,b,xtol=1e-12)
            best=max(best,f(rr))
    return best

rows=[];max_cap_err=0.;max_ramp_err=0.;lp_checks=0;ramp_checks=0
for alpha,lam,sigma,beta in product(ALPHAS,LAMBDAS,SIGMAS,CAPRATIOS):
    q=alpha*p;d=lam*h0;h=math.sqrt(3)*sigma*mu;U=beta*B
    R,D=reqs(h,p,q,d)
    b_resource=B-R;u_resource=U-D
    F=min(b_resource,u_resource)
    G,rstar=analytic_ramp(h,p,q,d,B,U)
    assert F+1e-7>=G,(alpha,lam,sigma,beta,F,G)
    budget_ok=b_resource>=Hneed-1e-7
    exposure_ok=u_resource>=Hneed-1e-7
    full_ok=budget_ok and exposure_ok
    ramp_ok=G>=Hneed-1e-7
    if ramp_ok: assert full_ok
    # distinguish negative gross bounds from a feasible zero-advance, no truncation
    positive_advance=F>=0
    premium=F-G
    region=('both_met' if full_ok else 'budget_missing' if not budget_ok and exposure_ok else 'exposure_missing' if budget_ok and not exposure_ok else 'both_missing')
    binding=('budget' if b_resource<u_resource-1e-7 else 'exposure' if u_resource<b_resource-1e-7 else 'tie')
    rows.append(dict(alpha=alpha,physical_delta_over_h0=lam,relative_sd=sigma,hard_cap_over_B=beta,delta_over_current_h=d/h,h_tonnes=h,delta_tonnes=d,q_eur_t=q,Rdelta_eur=R,Ddelta_eur=D,budget_capacity_eur=b_resource,exposure_capacity_eur=u_resource,sharp_advance_eur=F,ramp_advance_eur=G,ramp_optimal_slope=rstar,complexity_break_even_cost_eur=premium,budget_gate_pass=budget_ok,exposure_gate_pass=exposure_ok,fixed_cash_need_full_class=full_ok,fixed_cash_need_optimized_ramp=ramp_ok,full_only_meets_need=full_ok and not ramp_ok,positive_advance_possible=positive_advance,region=region,binding=binding))
    # Independent optimization of the financeability constraints (given the proved R/D result)
    c=linprog([-1.],A_ub=[[1.],[1.]],b_ub=[b_resource,u_resource],bounds=[(0,None)],method='highs')
    if F<-1e-6: assert not c.success
    else:
        assert c.success,(alpha,lam,sigma,beta,F,c.message)
        err=abs(c.x[0]-F);max_cap_err=max(err,max_cap_err)
        assert err<2e-6
    lp_checks+=1
    # Direct slope numerical optimization samples across the full joint parameter grid
    if (len(rows)-1)%3==0:
        number=numerical_ramp(h,p,q,d,B,U)
        err=abs(number-G);max_ramp_err=max(err,max_ramp_err)
        assert err<1e-4,(alpha,lam,sigma,beta,err)
        ramp_checks+=1

# aggregate factual result table, not probability claims
summary=[]
for a in ALPHAS:
    ss=[r for r in rows if r['alpha']==a]
    counts={k:sum(1 for x in ss if x['region']==k) for k in ('budget_missing','exposure_missing','both_missing','both_met')}
    summary.append(dict(alpha=a,n=len(ss),**counts,full_only=sum(x['full_only_meets_need'] for x in ss),positive_complexity_gain=sum(x['complexity_break_even_cost_eur']>0.005 and x['sharp_advance_eur']>=0 for x in ss),ramp_meets=sum(x['fixed_cash_need_optimized_ramp'] for x in ss)))

# predeclared exact baseline values, and numerical inequalities with interval monotonicity
cent=[r for r in rows if r['alpha']==.75 and r['physical_delta_over_h0']==.25 and r['relative_sd']==.10 and r['hard_cap_over_B']==1.10]
assert len(cent)==1
cent=cent[0]
assert abs(cent['sharp_advance_eur']-6638497.612043122)<1e-6
assert abs(cent['budget_capacity_eur']-6638497.612043122)<1e-6
assert abs(cent['exposure_capacity_eur']-6659096.179268995)<1e-6
assert cent['region']=='both_missing'
# monotonicity in U for all other parameters and monotonicity in physical delta
u_cases=0;delta_cases=0;sigma_cases=0;sigma_capacity_increases=0;sigma_span_decreases=0
for alpha,lam,sigma in product(ALPHAS,LAMBDAS,SIGMAS):
    subset=sorted((x for x in rows if x['alpha']==alpha and x['physical_delta_over_h0']==lam and x['relative_sd']==sigma),key=lambda x:x['hard_cap_over_B'])
    assert len(subset)==len(CAPRATIOS)
    for l,r in zip(subset,subset[1:]):
        assert r['sharp_advance_eur']+1e-7>=l['sharp_advance_eur']
        u_cases+=1
for alpha,sigma,beta in product(ALPHAS,SIGMAS,CAPRATIOS):
    subset=sorted((x for x in rows if x['alpha']==alpha and x['relative_sd']==sigma and x['hard_cap_over_B']==beta),key=lambda x:x['physical_delta_over_h0'])
    for l,r in zip(subset,subset[1:]):
        assert r['sharp_advance_eur']<=l['sharp_advance_eur']+1e-7
        delta_cases+=1
for alpha,lam,beta in product(ALPHAS,LAMBDAS,CAPRATIOS):
    subset=sorted((x for x in rows if x['alpha']==alpha and x['physical_delta_over_h0']==lam and x['hard_cap_over_B']==beta),key=lambda x:x['relative_sd'])
    for l,r in zip(subset,subset[1:]):
        # Mean reserve necessarily rises at fixed delta and q, but payment span need not.
        assert r['Rdelta_eur']>=l['Rdelta_eur']-1e-7
        if r['Ddelta_eur']<l['Ddelta_eur']-1e-6:sigma_span_decreases+=1
        if r['sharp_advance_eur']>l['sharp_advance_eur']+1e-6:sigma_capacity_increases+=1
        sigma_cases+=1

# Independent discretized-slope LP with exactly integrated uniform window weights.
# This tests the counterintuitive finite-span reversal without calling reqs/analytic_ramp.
def direct_slope_lp_normalized(L,q,B,U,N=64):
    W=2.;dx=1./N;cells=int(round((W+L)/dx))
    if abs(cells*dx-(W+L))>1e-9: raise ValueError('aligned grid required')
    lo=np.arange(cells)*dx;hi=lo+dx;z=np.arange(int(round(L/dx))+1)*dx
    over=np.maximum(0.,np.minimum(hi[None,:],(z+W)[:,None])-np.maximum(lo[None,:],z[:,None]))
    weights=np.empty(cells)
    for j,(a,b) in enumerate(zip(lo,hi)):
        if b<=L+1e-12: weights[j]=b-a
        else: weights[j]=((L+W)*(b-a)-(b*b-a*a)/2)/W
    mat=np.zeros((len(z)+2,cells+1));rhs=np.zeros(len(z)+2)
    mat[:len(z),:cells]=-over/W;rhs[:len(z)]=-q
    mat[-2,:cells]=weights;mat[-2,-1]=1.;rhs[-2]=B
    mat[-1,:cells]=dx;mat[-1,-1]=1.;rhs[-1]=U
    res=linprog(np.r_[np.zeros(cells),-1.],A_ub=mat,b_ub=rhs,bounds=[(0.,1.)]*cells+[(0.,None)],method='highs')
    assert res.success,res.message
    return float(res.x[-1]),float(np.max(mat@res.x-rhs))

independent_lp_cases=[]
for alpha,physical_lam,sigma,beta in [(.25,.75,.075,1.05),(.25,.75,.10,1.05),(.25,1.,.10,1.10)]:
    record=next(x for x in rows if x['alpha']==alpha and x['physical_delta_over_h0']==physical_lam and x['relative_sd']==sigma and x['hard_cap_over_B']==beta)
    scale=record['h_tonnes']*p
    lp_val,violation=direct_slope_lp_normalized(2*record['delta_over_current_h'],alpha,B/scale,beta*B/scale)
    err=abs(lp_val*scale-record['sharp_advance_eur'])
    assert err<1e-5 and violation<1e-8,(record,err,violation)
    independent_lp_cases.append({'alpha':alpha,'delta_over_h0':physical_lam,'sigma_rel':sigma,'U_over_B':beta,'lp_advance_eur':lp_val*scale,'error_eur':err,'primal_violation_normalized':violation})

columns=list(rows[0].keys())
with (ROOT/'multifactor_sensitivity_720.csv').open('w',newline='',encoding='utf8') as f:
    wr=csv.DictWriter(f,columns);wr.writeheader();wr.writerows(rows)
with (ROOT/'multifactor_summary.csv').open('w',newline='',encoding='utf8') as f:
    wr=csv.DictWriter(f,list(summary[0]));wr.writeheader();wr.writerows(summary)
result=dict(status='PASS',timestamp='2026-10-26',scope='Illustrative full-factorial stress grid; counts are NOT real-world frequencies or probabilities',
  design={'alpha_q_over_p':ALPHAS,'physical_delta_over_fixed_reference_h0':LAMBDAS,'relative_sd':SIGMAS,'hard_U_over_fixed_B':CAPRATIOS,'p':p,'B':B,'fixed_cash_need':Hneed,'reference_h0_tonnes':h0,'donor_cost_loading':0,'fixed_cost':0},
  grid_cases=len(rows),lp_checks=lp_checks,ramp_numeric_checks=ramp_checks,max_lp_error_eur=max_cap_err,max_numeric_ramp_error_eur=max_ramp_err,
  monotonicity_checks=dict(raising_U=u_cases,raising_fixed_delta=delta_cases,raising_sigma_compared=sigma_cases,sigma_capacity_increases=sigma_capacity_increases,sigma_span_decreases=sigma_span_decreases),
  region_counts={k:sum(1 for x in rows if x['region']==k) for k in ('budget_missing','exposure_missing','both_missing','both_met')},
  full_only_meets=sum(x['full_only_meets_need'] for x in rows),ramp_meets=sum(x['fixed_cash_need_optimized_ramp'] for x in rows),
  positive_complexity_gain_feasible_zero_advance=sum(x['complexity_break_even_cost_eur']>0.005 and x['sharp_advance_eur']>=0 for x in rows),
  combined_summary=summary,independent_slope_lp_checks=independent_lp_cases,baseline={k:cent[k] for k in ('sharp_advance_eur','ramp_advance_eur','budget_capacity_eur','exposure_capacity_eur','region','binding')},
  limitations=['Fixed cash requirement and other inputs are scenarios, not estimated parameters','Nonlinear co-variation or empirical reporting distribution not identified','Strictly illustrative parameter grid; no probability assigned to outcomes','For h changes physical delta is held fixed at its design level rather than assuming constant delta/h'])
(ROOT/'multifactor_audit.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
print(json.dumps({k:v for k,v in result.items() if k not in ('design','limitations')},indent=2))
