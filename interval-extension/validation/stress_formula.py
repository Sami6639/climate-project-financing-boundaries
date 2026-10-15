"""Checks independent of the symbolic proof, including exact window extrema.
All residuals reported in normalized money / slope units.
"""
import numpy as np,json
from validate_interval import solve,OUT

def theory(a,r):
 n=int(np.floor(a)); f=a-n
 D=2*((n+1)*r+max(0,r-1+f)); R=r*r+2*a*r
 return D,R

def exact_contract(a,r):
 end=2+2*a
 knots=sorted(set([0.,end]+[x for k in range(int(np.floor(a))+2) for x in (2*k+2-2*r,2*k+2) if 0<=x<=end]))
 l=np.array(knots[:-1]); z=np.array(knots[1:]); mid=(l+z)/2
 g=((np.mod(mid,2)>=2-2*r)&(r>0)).astype(float)
 centers=np.unique(np.clip(np.r_[0,2*a,l,z,l-2,z-2],0,2*a))
 windows=np.maximum(0,np.minimum(z[None,:],centers[:,None]+2)-np.maximum(l[None,:],centers[:,None]))@g
 # E T at uppercenter; survival weight=1 before final window.
 def primitive(t):
  boundary=2*a
  return np.where(t<=boundary,t,boundary+(end*(t-boundary)-(t*t-boundary*boundary)/2)/2)
 reserve=(primitive(z)-primitive(l))@g
 mass=(z-l)@g
 return mass,reserve,float(min(windows)),float(max(windows)),len(centers)

def run():
 exact=[]
 for a in [0,1e-8,.1,.25,.5,.75,.99999999,1,1.00000001,1.4,1.99999999,2,2.00000001,3.25,5,8.25]:
  for r in [0,.01,.25,.5,.75,.99,1]:
   D,R=theory(a,r); mass,reserve,wmin,wmax,n=exact_contract(a,r)
   exact.append(dict(a=a,r=r,mass_error=mass-D,reserve_error=reserve-R,window_min_error=wmin-2*r,window_max_error=wmax-2*r,breakpoints=n))
 lp=[]
 rng=np.random.default_rng(20261014)
 cases=[(a,r) for a in [0,.25,.999,1,1.001,2,2.75,4] for r in [0,.25,.75,1]]
 cases += [(float(rng.uniform(0,5)),float(rng.uniform(0,1))) for _ in range(12)]
 for i,(a,r) in enumerate(cases):
  D,R=theory(a,r)
  # Alternate binding mechanisms and tie cases.
  b=R+2+([0,.3,.7][i%3]); u=D+2+([.7,.3,0][i%3])
  c=min(b-R,u-D)
  low=solve(a,r,b,u,400); high=solve(a,r,b,u,400,True)
  lp.append(dict(a=a,r=r,b=b,u=u,theory_c=c,lower_c=low['c'],upper_c=high['c'],theory_below_lower=max(0,low['c']-c),theory_above_upper=max(0,c-high['c']),primal_violation=low['max_primal_violation'],duality_gap=low['primal_dual_gap']))
 out=dict(exact_construction=exact,lp_cases=lp,summary={'construction_cases':len(exact),'lp_cases':len(lp),'max_exact_error':max(abs(v) for row in exact for k,v in row.items() if 'error' in k),'max_bracket_violation':max(max(x['theory_below_lower'],x['theory_above_upper']) for x in lp),'max_primal_violation':max(x['primal_violation'] for x in lp),'max_duality_gap':max(x['duality_gap'] for x in lp)})
 (OUT/'stress_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out['summary'],indent=2))
if __name__=='__main__':run()
