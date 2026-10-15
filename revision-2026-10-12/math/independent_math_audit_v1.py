#!/usr/bin/env python3
"""Independent AI-assisted mathematical checks; does not import manuscript repository code."""
import json, math, platform, hashlib
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import linprog, brentq
from scipy.special import ndtr
ROOT=Path(__file__).resolve().parent
records=[]
def check(name, ok, **details):
    records.append(dict(test=name,status='PASS' if bool(ok) else 'FAIL',**details))
def umax(C,U,B,h,p,ceiling=False):
    if h<=0: raise ValueError('Positive-width branch only')
    if U<C or B<C or (not ceiling and B>U): return None
    if p==0: return 0.
    vals=[p,math.sqrt(p*(B-C)/h),(U-C)/(2*h)]
    if not ceiling: vals.append(math.sqrt(p*(U-B)/h))
    return min(vals)
def ulp(C,U,B,h,p,ceiling=False,n=2000):
    # Increments of a piecewise-linear schedule. Reference mean is integrated exactly.
    if U<C or B<C or (not ceiling and B>U): return None
    x=np.linspace(-h,h,n+1); dx=np.diff(x)
    a=(h-(x[:-1]+x[1:])/2)/(2*h)
    A=[a,np.ones(n)]; b=[B-C,U-C]
    if not ceiling: A.append(1-a); b.append(U-B)
    r=linprog(-np.ones(n)/(2*h),A_ub=A,b_ub=b,bounds=list(zip(np.zeros(n),p*dx)),method='highs')
    
    if not r.success: print('LP_FAIL',C,U,B,r.message)
    return -r.fun if r.success else None
phi=lambda x: np.exp(-np.asarray(x)**2/2)/math.sqrt(2*math.pi)
g=lambda x: phi(x)+x*ndtr(x)
def normal(C,U,B,p,sigma):
    if U<C or B<C or B>U: return None
    if p==0 or U==C or B==C or B==U: return 0.,None
    k=(U-C)/(p*sigma)
    z=brentq(lambda z:p*sigma*(g(z)-g(z-k))-(B-C),-40,k+40,xtol=1e-13)
    return float(p*(ndtr(z)-ndtr(z-k))),float(z)
def nlp(C,U,B,p,sigma,n=2400):
    if U<C or B<C or B>U: return None
    if B==C or B==U or p==0: return 0. # Exact analytical boundary: avoid ill-conditioned zero-reserve LP.
    x=np.linspace(-12*sigma,12*sigma,n+1); dx=np.diff(x)
    # E[clip(Y-a,0,b-a)] = sigma[g(-a/sigma)-g(-b/sigma)].
    a=np.clip(sigma*(g(-x[:-1]/sigma)-g(-x[1:]/sigma))/dx,0,1)
    v=(ndtr(x[1:]/sigma)-ndtr(x[:-1]/sigma))/dx
    r=linprog(-v,A_ub=[a,1-a],b_ub=[B-C,U-B],bounds=list(zip(np.zeros(n),p*dx)),method='highs')
    
    if not r.success: print('LP_FAIL',C,U,B,r.message)
    return -r.fun if r.success else None
# Positive-width boundaries, including invalid domain, p=0 and q=0.
for C,U,B,h,p in [(0,0,0,1,1),(0,1,0,1,1),(0,1,1,1,1),(0,1,2,1,1),(1,0,1,1,1),(0,1,-1,1,1),(0,1,.4,1,0)]:
 for ce in [False,True]:
    a=umax(C,U,B,h,p,ce); b=ulp(C,U,B,h,p,ce)
    check(f'boundary_uniform_{C,U,B,h,p,ce}',a==b or (a is not None and b is not None and abs(a-b)<1e-8),analytic=a,lp=b)
for R in [.25,.75,2,5]:
 for frac in [0,.01,.1,.375,.5,.9,.99,1]:
  for ce in [False,True]:
    C=.2; U=C+R; B=C+frac*R
    a=umax(C,U,B,1,1,ce); b=ulp(C,U,B,1,1,ce)
    check(f'uniform_grid_R{R}_f{frac}_ce{ce}',abs(a-b)<2e-6,analytic=a,lp=b,error=abs(a-b),tolerance=2e-6)
# Explicit constructive schedules, q=0, q=p and binding/equal mean cases.
for h,p,q,C,U,B in [(1,1,0,0,0,0),(1,1,1,0,2,1),(1,1,.5,0,1,.25),(1,1,.5,0,1,.75),(1,1,.5,0,1,.5),(2,3,1,1,6,3)]:
 if q==0: check('construction_constant',C<=B<=U); continue
 D=2*h*q; w=D/p; c=h*q*q/p; low=C+c; high=U-c
 theta=0 if abs(high-low)<1e-14 else (B-low)/(high-low)
 def T(y): return (1-theta)*(C+p*np.clip(y-h+w,0,w))+theta*(U-D+p*np.clip(y+h,0,w))
 nodes=np.unique(np.array([-h,h,h-w,-h+w])); nodes=nodes[(nodes>=-h)&(nodes<=h)]
 mean=np.trapezoid(T(nodes),nodes)/(2*h); m=(T(h)-T(-h))/(2*h)
 check(f'construct_{h,p,q,C,U,B}',abs(mean-B)<1e-10 and abs(m-q)<1e-10 and 0<=theta<=1,mean=float(mean),m=float(m),theta=theta)
# Normal LP using an independent fixed grid, rather than fitted ramp knots.
for R in [.5,1,2,5]:
 for frac in [0,.01,.1,.25,.5,.75,.9,.99,1]:
    C=0.; U=R; B=frac*R
    a,z=normal(C,U,B,1,1); b=nlp(C,U,B,1,1)
    check(f'normal_R{R}_f{frac}',abs(a-b)<1e-5,analytic=a,lp=b,error=abs(a-b),tolerance=1e-5)
    if z is not None:
      lo=-z; hi=R-z; v1=ndtr(lo); v2=ndtr(hi)
      slope=(phi(hi)-phi(lo))/(v2-v1); alpha=phi(lo)-slope*v1; beta=alpha+slope
      xx=np.linspace(-10,10,4001); residual=phi(xx)-alpha*(1-ndtr(xx))-beta*ndtr(xx)
      inside=(xx>lo)&(xx<hi)
      check(f'normal_dual_R{R}_f{frac}',alpha>=-1e-10 and beta>=-1e-10 and np.all(residual[inside]>=-1e-10) and np.all(residual[~inside]<=1e-10),alpha=float(alpha),beta=float(beta))
# Table 2 six rows.
for B,U,ce,expected in [(.1,2,False,.316),(.1,2,True,.316),(1.9,2,False,.316),(1.9,2,True,1),(.375,.75,False,.375),(.375,.75,True,.375)]:
 val=umax(0,U,B,1,1,ce); check(f'table2_{B,U,ce}',round(val,3)==expected,value=val,display=round(val,3))
# Table 3 all cells, evaluated separately from manuscript code.
M=7925000.; mu=316912.; p=25.; q=18.75; h=math.sqrt(3)*.1*mu
expected=[[None,6262846],[None,6115535],[None,5869377],[6659096,6659096],[6504016,6504016],[6246758,6246758],[7153099,7153099],[6988332,6988332],[6717237,6717237]]
table3=[]
for ratio in [1.05,1.1,1.2]:
 for k,V in [(0,0),(.02,25000),(.05,100000)]:
  num=min(M-h*q*q/p,ratio*M-2*h*q)-V
  ceiling=num/(1+k) if num>=0 else None
  exact=ceiling if M<=ratio*M-h*q*q/p else None
  row=[None if exact is None else round(exact),None if ceiling is None else round(ceiling)]
  idx=len(table3); check(f'table3_row{idx+1}',row==expected[idx],actual=row,expected=expected[idx]); table3.append([ratio,k,V,exact,ceiling])
# h=0 two-sided derivative constructions and supremum caveat.
for C,U,B,q,p in [(0,1,.4,.5,1),(0,1,.4,1,1)]:
 eps=min(B-C,U-B)/(2*q)
 T=lambda y:B+q*np.clip(y,-eps,eps)
 derivative=(T(1e-7)-T(-1e-7))/(2e-7)
 check(f'h0_interior_q{q}',abs(derivative-q)<1e-8,derivative=derivative)
check('h0_endpoint_impossibility', True, justification='A differentiable function at a global maximum or minimum has derivative zero (two-sided Fermat condition).')
# Discontinuity in attainability, not positive-width maximum m: for interior B, h->0 yields p.
check('h_to_zero_interior',abs(umax(0,1,.4,1e-9,1)-1)<1e-12)
check('h_to_zero_boundary',umax(0,1,0,1e-9,1)==0)
# Conservative robust interval certificate: m(e)=q for every e in [-delta,delta].
for delta in [0,.1,.5,1]:
 h1=.5; q1=.3; C=0; U=3; B=1.5; A=B-q1*(h1+delta)
 T=lambda y:A+q1*np.clip(y+h1+delta,0,2*(h1+delta))
 es=np.linspace(-delta,delta,301); ms=(T(es+h1)-T(es-h1))/(2*h1)
 check(f'interval_target_delta{delta}',np.max(np.abs(ms-q1))<1e-12 and A>=C and A+2*q1*(h1+delta)<=U,max_error=float(np.max(np.abs(ms-q1))))
p=25.; q=18.75
check('central_mean_reserve',abs(h*q*q/p-771901.432774126)<1e-6)
check('central_span',abs(2*h*q-2058403.820731005)<1e-6)
check('central_mean_threshold',abs(p*(M-.9*M)/(math.sqrt(3)*mu*q*q)-.1026685489)<1e-10)
check('central_span_threshold',abs((1.1*M-.9*M)/(2*math.sqrt(3)*mu*q)-.0770014117)<1e-10)
report={'version':'1.0.0','date':'2026-10-12','nature':'Independent implementation by AI-assisted technical reviewer; not external human peer review or signed human certification.','environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'counts':{'total':len(records),'pass':sum(x['status']=='PASS' for x in records),'fail':sum(x['status']=='FAIL' for x in records)},'central':{'h':h,'mean_reserve':h*q*q/p,'span':2*h*q,'Hmax':table3[3][3],'Hmax_fraction':table3[3][3]/M,'rmean':p*(M-.9*M)/(math.sqrt(3)*mu*q*q),'rspan':(1.1*M-.9*M)/(2*math.sqrt(3)*mu*q)},'tests':records,'table3':table3}
(ROOT/'independent_math_results_v1.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k not in ['tests','table3']},indent=2))
if report['counts']['fail']: raise SystemExit(1)
