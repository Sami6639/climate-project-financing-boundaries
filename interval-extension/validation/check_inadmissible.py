import json
from validate_interval import solve,OUT
from stress_formula import theory
out=[]
for a in [0,.1,1,3]:
 for r in [1.01,1.5]:
  x=solve(a,r,100,100,100)
  out.append({'a':a,'r':r,'expected':'infeasible because q>p','actual_infeasible':not x['success']})
for a,r in [(.25,.75),(1,1),(3,.75)]:
 D,R=theory(a,r)
 x=solve(a,r,R-.1,D+.5,400)
 out.append({'a':a,'r':r,'expected':'no nonnegative secured advance','unconstrained_floor':x['c'],'actual_no_nonnegative_advance':x['c']<0})
(OUT/'inadmissible_cases.json').write_text(json.dumps(out,indent=2))
assert all(x.get('actual_infeasible',x.get('actual_no_nonnegative_advance')) for x in out)
print('All 11 infeasibility checks passed.')
