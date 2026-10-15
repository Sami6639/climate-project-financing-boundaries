"""Independent node-value LP, tuned capped-ramp benchmark and decision thresholds.
No imports from the original implementation. Deterministic scenarios, not UEFA estimates.
"""
from pathlib import Path
import numpy as np, math, json, csv
from scipy.optimize import linprog, brentq
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.sparse import diags, hstack, csr_matrix
R=Path(__file__).parent

def save(name,rows):
 with (R/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
def frontier(B,U,C=0,h=1,p=1,ceiling=False):
 if h<=0 or p<0:raise ValueError('Use separate zero-width theorem branch')
 if B<C or U<C or (not ceiling and B>U):return float('nan')
 return min([p,math.sqrt(p*(B-C)/h),(U-C)/(2*h)]+([] if ceiling else [math.sqrt(p*(U-B)/h)]))
def rampmean(t,U,C=0,h=1,p=1):
 # Integral of C+p[(y-t)+ -(y-t-(U-C)/p)+] over uniform [-h,h].
 def hinge(a):return (max(h-a,0)**2-max(-h-a,0)**2)/(4*h)
 return C+p*(hinge(t)-hinge(t+(U-C)/p))
def tune(B,U,C=0,h=1,p=1,ceiling=False):
 target=min(B,(U+C)/2) if ceiling else B
 if target==C:return (h,0,target)
 if target==U:return (-h-(U-C)/p,0,target)
 t=brentq(lambda t:rampmean(t,U,C,h,p)-target,-h-(U-C)/p,h,xtol=1e-12)
 T=lambda x:C+np.clip(p*(x-t),0,U-C)
 return t,float((T(h)-T(-h))/(2*h)),rampmean(t,U,C,h,p)
def nodal_lp(B,U,C=0,h=1,p=1,ceiling=False,n=1024):
 # Decision variables are payment levels at n+1 nodes, not slopes.
 dx=2*h/n;D=diags([-np.ones(n),np.ones(n)],[0,1],shape=(n,n+1),format='csr')
 from scipy.sparse import vstack
 A=vstack([D,-D],format='csr');b=np.r_[np.full(n,p*dx),np.zeros(n)]
 mean=np.ones(n+1)/n;mean[[0,-1]]/=2
 c=np.zeros(n+1);c[0]=1/(2*h);c[-1]=-1/(2*h)
 if ceiling:A=vstack([A,csr_matrix(mean[None,:])]);b=np.r_[b,B];eq=None;beq=None
 else:eq=csr_matrix(mean[None,:]);beq=[B]
 z=linprog(c,A_ub=A,b_ub=b,A_eq=eq,b_eq=beq,bounds=[(C,U)]*(n+1),method='highs')
 return z.success,(-z.fun if z.success else None)
rows=[]
for ceiling in [False,True]:
 for U in [.25,.5,.75,1,1.5,2,4]:
  for a in [0,.005,.025,.1,.25,.5,.75,.9,.975,.995,1]:
   B=a*U;t,m,avg=tune(B,U,ceiling=ceiling);f=frontier(B,U,ceiling=ceiling)
   ok,l=nodal_lp(B,U,ceiling=ceiling)
   rows.append(dict(regime='Ceiling' if ceiling else 'Exact',cap=U,budget=B,threshold=t,ramp_mean=avg,ramp_marginal=m,formula_max=f,nodal_lp_max=l,ramp_formula_error=abs(m-f),formula_lp_gap=f-l))
save('uniform_tuned_ramp_benchmark.csv',rows)
# Gaussian ramp mean and local derivative computed by adaptive quadrature.
phi=lambda z:math.exp(-z*z/2)/math.sqrt(2*math.pi)
g=lambda z:phi(z)+z*ndtr(z)
grows=[]
for U in [.25,.5,1,2,4,8]:
 for a in [.01,.1,.25,.5,.9,.99]:
  B=a*U
  def expectation(t):
   return quad(lambda x:(x-t)*phi(x),t,t+U,epsabs=1e-11)[0]+U*ndtr(-(t+U))
  t=brentq(lambda t:expectation(t)-B,-40-U,40,xtol=1e-12)
  m=quad(phi,t,t+U,epsabs=1e-12)[0]
  z=brentq(lambda z:g(z)-g(z-U)-B,-40,40+U,xtol=1e-12)
  optimum=ndtr(z)-ndtr(z-U)
  grows.append(dict(cap=U,mean=B,ramp_threshold=t,ramp_mean=expectation(t),ramp_marginal=m,formula_max=optimum,absolute_gap=abs(m-optimum)))
save('gaussian_tuned_ramp_benchmark.csv',grows)
# Invalid-domain and exact-boundary feasibility checked against node LP.
bound=[]
for C,B,U in [(1,0,2),(1,1,2),(1,2,2),(1,3,2),(1,2,.5),(1,1,1),(0,0,0)]:
 for ceiling in [False,True]:
  if U<C:ok=False;l=None
  else:ok,l=nodal_lp(B,U,C,ceiling=ceiling,n=64)
  f=frontier(B,U,C,ceiling=ceiling)
  bound.append(dict(floor=C,budget=B,cap=U,ceiling=ceiling,formula_feasible=math.isfinite(f),lp_feasible=ok,formula=f,lp=l,passed=(math.isfinite(f)==ok and (not ok or abs(f-l)<1e-8))))
save('boundary_domain_audit.csv',bound)
B=7925000.;mu=316912.;p=25.;q=18.75;h=math.sqrt(3)*mu*.1;U=1.1*B;H=.9*B
reserve=h*q*q/p;span=2*h*q;Hmax=min(B-reserve,U-span)
t,m,avg=tune(B,U,Hmax,h,p)
central={'B_eur':B,'U_eur':U,'target_q':q,'H_test':H,'Hmax':Hmax,'ramp_threshold_report_tonnes':mu+t,'ramp_cap_start_report_tonnes':mu+t+(U-Hmax)/p,'ramp_expected_expenditure':avg,'ramp_marginal':m,'target_met_abs_error':abs(m-q),'cap_increase_for_90pct':H+span-U,'advance_reduction_from_90pct':H-Hmax,'cap_at_mean_limited_transition':B+span-reserve,'cap_ratio_at_mean_limited_transition':(B+span-reserve)/B,'rate_needed_for_mean_at_90pct':h*q*q/(B-H),'max_target_at_90pct':frontier(B,U,H,h,p),'sd_span_threshold':(U-H)/(2*math.sqrt(3)*mu*q),'sd_mean_threshold':(B-H)*p/(math.sqrt(3)*mu*q*q),'fixed_cost_max_at_80pct_no_loading':Hmax-.8*B,'loading_max_at_80pct_zero_fixed_cost':Hmax/(.8*B)-1}
(R/'new_analysis_summary.json').write_text(json.dumps(central,indent=2))
summary={'uniform_benchmark_cases':len(rows),'uniform_max_ramp_formula_gap':max(x['ramp_formula_error'] for x in rows),'independent_nodal_lp_cases':len(rows),'nodal_max_formula_gap':max(x['formula_lp_gap'] for x in rows),'nodal_min_formula_gap':min(x['formula_lp_gap'] for x in rows),'gaussian_quadrature_benchmark_cases':len(grows),'gaussian_max_ramp_formula_gap':max(x['absolute_gap'] for x in grows),'boundary_cases':len(bound),'boundary_pass':all(x['passed'] for x in bound),'central_ramp_target_error':abs(m-q)}
assert summary['uniform_max_ramp_formula_gap']<1e-8
assert summary['gaussian_max_ramp_formula_gap']<1e-8
assert summary['nodal_min_formula_gap']>-1e-7
assert summary['nodal_max_formula_gap']<3e-6
assert summary['boundary_pass']
assert abs(m-q)<1e-8
(R/'independent_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'validation':summary,'central':central},indent=2))
# Independent advance maximization: payment-node variables and a separate H variable.
# Scale euros by B and emissions by h to improve numerical conditioning.
from scipy.sparse import vstack
advance=[]
for ratio in [1.05,1.1,1.2]:
 for k,V in [(0,0),(.02,25000),(.05,100000)]:
  for ceiling in [False,True]:
   n=2048;dx=2/n;D=diags([-np.ones(n),np.ones(n)],[0,1],shape=(n,n+1),format='csr')
   zeros=csr_matrix((n,1));A=vstack([hstack([D,zeros]),hstack([-D,zeros])]);bb=np.r_[np.full(n,p*h/B*dx),np.zeros(n)]
   # T(-1) >= (1+k)H+V, and endpoint difference >= 2hq.
   floor=np.zeros(n+2);floor[0]=-1;floor[-1]=1+k
   marginal=np.zeros(n+2);marginal[0]=1;marginal[n]=-1
   A=vstack([A,csr_matrix(floor[None,:]),csr_matrix(marginal[None,:])]);bb=np.r_[bb,-V/B,-2*h*q/B]
   mean=np.r_[np.ones(n+1)/n,0];mean[0]/=2;mean[n]/=2
   if ceiling:A=vstack([A,csr_matrix(mean[None,:])]);bb=np.r_[bb,1];eq=None;beq=None
   else:eq=csr_matrix(mean[None,:]);beq=[1]
   cost=np.zeros(n+2);cost[-1]=-1
   z=linprog(cost,A_ub=A,b_ub=bb,A_eq=eq,b_eq=beq,bounds=[(0,ratio)]*(n+1)+[(0,None)],method='highs')
   feasible=ceiling or ratio*B-B>=reserve
   formula=(min(B-reserve,ratio*B-span)-V)/(1+k) if feasible else None
   lp=-z.fun*B if z.success else None
   passed=(z.success==feasible and (not feasible or abs(formula-lp)<1.0))
   advance.append(dict(cap_ratio=ratio,loading=k,fixed_cost_eur=V,ceiling=ceiling,formula_advance_eur=formula,node_lp_advance_eur=lp,gap_eur=None if not feasible else formula-lp,passed=passed))
save('advance_optimization_audit.csv',advance)
assert all(a['passed'] for a in advance)
print('Independent advance maximizations:18/18 pass; maximum euro gap:',max(a['gap_eur'] or 0 for a in advance))
