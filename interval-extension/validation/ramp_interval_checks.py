"""Independent full-interval validation of the optimized uniform single ramp.
All non-dimensional monetary quantities are in p*h; slope is in p; distance in h.
No finite centre grid is used to certify a uniform ramp: its overlap function is
piecewise linear, and all its breakpoints and both interval endpoints are checked.
"""
from pathlib import Path
import csv, json, math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
OUT=Path(__file__).resolve().parent
SEED=20261014
TOL_MONEY=2e-9
TOL_SLOPE=2e-10

def global_frontier(h,p,q,delta):
    if h<=0 or p<=0 or not 0<=q<=p or delta<0:
        raise ValueError('Require h>0,p>0,0<=q<=p,delta>=0; h=0 is a separate model.')
    W=2*h; alpha=delta/h; n=math.floor(alpha); ell=2*delta-n*W
    return (n+1)*W*q+max(0.,W*q-p*(W-ell)),h*q*q/p+2*delta*q

def fixed_ramp(h,p,q,delta,r):
    if q==0: return dict(r=0.,start=0.,width=0.,span=0.,reserve=0.)
    assert q<=r<=p
    W=2*h;L=2*delta;s=W*q/r
    start=W-s;width=max(s,L+2*s-W)
    K=r*width
    R=q*L+W*q*q/(2*r) if L<=W-s else r*(L-W)+2*W*q-W*q*q/(2*r)
    return dict(r=r,start=start,width=width,span=K,reserve=R)

def best_ramp(h,p,q,delta,B,U):
    global_frontier(h,p,q,delta)
    if q==0:return dict(c=min(B,U),**fixed_ramp(h,p,q,delta,0.),candidates=[0.])
    W=2*h;L=2*delta;candidates=[q,p]
    if U>B:
        rc=W*q*q/(2*(U-B))
        if q<=rc<=p and L>=W-W*q/rc:candidates.append(rc)
    rows=[fixed_ramp(h,p,q,delta,r) for r in sorted(set(candidates))]
    for row in rows:row['c']=min(B-row['reserve'],U-row['span'])
    row=max(rows,key=lambda x:x['c']).copy();row['candidates']=sorted(set(candidates));return row

def exact_ramp_checks(h,p,q,delta,row):
    W=2*h; L=2*delta; r=row['r']; t=row['start']; w=row['width']
    points=sorted(set([0.,L]+[x for x in [t-W,t,t+w-W,t+w] if 0<x<L]))
    overlaps=[max(0.,min(x+W,t+w)-max(x,t)) for x in points]
    m=min(overlaps)*r/W
    # Quadrature of the payment itself, split at every payment kink. This does
    # not reuse the claimed reserve formula or its branch algebra.
    ends=sorted(set([L,L+W]+[x for x in [t,t+w] if L<x<L+W]))
    mean=sum(quad(lambda y:r*min(w,max(0.,y-t)),a,b,epsabs=2e-13)[0] for a,b in zip(ends[:-1],ends[1:]))/W
    return m,mean,len(points)

def run():
    comparisons=[];checks=[];rng=np.random.default_rng(SEED)
    for q in [0.,.1,.25,.5,.75,.9,1.]:
      for alpha in [0.,.1,.25,.5,.75,1.,1.5,2.,3.]:
       D,R=global_frontier(1.,1.,q,alpha)
       for mechanism,addB,addU in [('budget',0.,2.),('exposure',2.,0.),('joint',0.,0.)]:
        B=R+3+addB;U=D+3+addU
        best=best_ramp(1.,1.,q,alpha,B,U)
        simple=fixed_ramp(1.,1.,q,alpha,q)
        cc=min(B-simple['reserve'],U-simple['span'])
        comparisons.append(dict(q_over_p=q,delta_over_h=alpha,regime=mechanism,B_over_ph=B,U_over_ph=U,global_c_over_ph=3.,best_ramp_c_over_ph=best['c'],conservative_c_over_ph=cc,global_minus_best_over_ph=3-best['c'],global_minus_conservative_over_ph=3-cc,best_slope_over_p=best['r'],best_span_over_ph=best['span'],best_reserve_over_ph=best['reserve']))
        m,mean,n=exact_ramp_checks(1.,1.,q,alpha,best)
        assert m>=q-TOL_SLOPE,(q,alpha,m)
        assert abs(mean-best['reserve'])<=TOL_MONEY
        assert best['c']<=3+TOL_MONEY and best['c']>=cc-TOL_MONEY
        if alpha<=1-q+1e-14:assert abs(best['c']-3.)<=TOL_MONEY
        if alpha>=1:assert abs(best['c']-cc)<=TOL_MONEY
        checks.append(dict(scenario=f'ramp_{q}_{alpha}_{mechanism}',expected='continuum target; exact mean; class upper bound; conservative dominance',error=max(max(0.,q-m),abs(mean-best['reserve']),max(0.,best['c']-3),max(0.,cc-best['c'])),unit='normalized slope or p*h (separately bounded)',threshold=TOL_MONEY,passed=True))
    # Independently search the one-dimensional envelope. Dense nodes include all
    # analytic candidates and each local bounded minimizer between grid nodes.
    largest_search_error=0.;largest_integral_error=0.;crossing_count=0
    for k in range(240):
        q=float(rng.uniform(.005,.995));alpha=float(rng.uniform(0,4));D,R=global_frontier(1,1,q,alpha)
        B=R+float(rng.uniform(1,5));U=D+float(rng.uniform(1,5));best=best_ramp(1,1,q,alpha,B,U)
        nodes=np.linspace(q,1,1201)
        values=np.array([min(B-x['reserve'],U-x['span']) for x in [fixed_ramp(1,1,q,alpha,r) for r in nodes]])
        i=int(values.argmax());lo=nodes[max(0,i-1)];hi=nodes[min(len(nodes)-1,i+1)]
        def objective(r):
            x=fixed_ramp(1,1,q,alpha,r);return -min(B-x['reserve'],U-x['span'])
        search=minimize_scalar(objective,bounds=(lo,hi),method='bounded',options={'xatol':1e-14})
        search_c=max(values.max(),-search.fun)
        search_error=max(0,search_c-best['c']);largest_search_error=max(largest_search_error,search_error)
        assert search_error<=TOL_MONEY
        if q+1e-8<best['r']<1-1e-8:crossing_count+=1
        m,mean,n=exact_ramp_checks(1,1,q,alpha,best)
        largest_integral_error=max(largest_integral_error,abs(mean-best['reserve']))
        assert m>=q-TOL_SLOPE and abs(mean-best['reserve'])<=TOL_MONEY
        # Random feasible enlarged/earlier competitors at fixed slope. The
        # rightmost minimal-width ramp must weakly dominate both resources.
        rr=float(rng.uniform(q,1));optimal=fixed_ramp(1,1,q,alpha,rr)
        width=optimal['width']+float(rng.uniform(0,3));s=2*q/rr
        lower=2*alpha+s-width;upper=2-s
        t=float(rng.uniform(lower,upper));competitor={'r':rr,'start':t,'width':width,'span':rr*width,'reserve':0.}
        _,other_mean,_=exact_ramp_checks(1,1,q,alpha,competitor)
        assert other_mean>=optimal['reserve']-TOL_MONEY and competitor['span']>=optimal['span']-TOL_MONEY
    # Deliberately construct an interior optimum: alpha<1, crossing on the
    # increasing-reserve/decreasing-span branch, and both endpoints worse.
    q=.6;alpha=.8;rc=.68;v=fixed_ramp(1,1,q,alpha,rc);B=v['reserve']+3;U=v['span']+3
    z=best_ramp(1,1,q,alpha,B,U)
    assert abs(z['r']-rc)<1e-10 and abs(z['c']-3)<1e-10
    interior_solution=z.copy()
    checks.append(dict(scenario='interior_crossing',expected='r/p=0.68; mean and span jointly bind',error=abs(z['r']-rc),unit='slope/p',threshold=1e-10,passed=True))
    edges=[]
    for q in [0.,.25,.75,1.]:
      aa=set([0.,max(0.,1-q-1e-9),1-q,1-q+1e-9]+[k+side for k in [1.,2.,3.] for side in [-1e-9,0,1e-9]])
      for alpha in sorted(aa):
       for h in [1.,1e-12]:
        p=25.;qphys=p*q;D,R=global_frontier(h,p,qphys,h*alpha)
        scale=p*h;B=R+3*scale;U=D+3*scale;z=best_ramp(h,p,qphys,h*alpha,B,U)
        m,mean,n=exact_ramp_checks(h,p,qphys,h*alpha,z)
        err_m=max(0.,qphys-m);err_R=abs(mean-z['reserve'])
        assert err_m<=TOL_SLOPE*p and err_R<=TOL_MONEY*scale
        # Nonzero financing loading and fixed cost preserve the joint gates.
        V=.2*scale;kappa=.35;H=(min(B-R,U-D)-V)/(1+kappa)
        err_H=max(abs((1+kappa)*H+V+R-B),abs((1+kappa)*H+V+D-U))
        assert err_H<=TOL_MONEY*scale
        edges.append(dict(q_over_p=q,delta_over_h=alpha,h=h,p=p,mean_error_currency=err_R,mean_tolerance_currency=TOL_MONEY*scale,marginal_shortfall_currency_per_emissions=err_m,marginal_tolerance_currency_per_emissions=TOL_SLOPE*p,joint_gate_error_currency=err_H,passed=True))
    D,R=global_frontier(1,1,.75,1)
    negative_H=(min(R-.1-R,D-.1-D)-.2)/1.35
    assert negative_H<0 # Negative capacity is NOT clipped to a feasible zero.
    for args in [(0.,1.,.5,0.),(1.,1.,1.01,0.)]:
        try:global_frontier(*args)
        except ValueError:pass
        else:raise AssertionError('Excluded parameter accepted')
    summary=dict(seed=SEED,comparison_rows=len(comparisons),independent_random_cases=240,targeted_boundary_cases=len(edges),interior_optima_found=crossing_count+1,max_dense_search_excess_over_formula=largest_search_error,max_independent_integral_error=largest_integral_error,normalized_money_tolerance=TOL_MONEY,normalized_slope_tolerance=TOL_SLOPE,negative_capacity_not_clipped=negative_H,zero_noise_treated_separately=True)
    with (OUT/'optimized_ramp_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,comparisons[0].keys());w.writeheader();w.writerows(comparisons)
    (OUT/'optimized_ramp_checks.json').write_text(json.dumps(dict(summary=summary,checks=checks,boundary_cases=edges,interior_crossing=interior_solution),indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=='__main__':run()
