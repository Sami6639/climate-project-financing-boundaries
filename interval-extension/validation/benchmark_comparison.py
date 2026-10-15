import numpy as np,json
from stress_formula import theory,OUT
h=np.sqrt(3)*.1*316912.;p=25.;q=18.75;B=7925000.;U=1.1*B;s=p*h
rows=[]
for a in [0,.1,.25,.5,.75,1,1.25,2,3]:
 D,R=theory(a,q/p);D*=s;R*=s
 Do=2*q*h*(1+a);Ro=q*h*(1+2*a)
 H=min(B-R,U-D);Ho=min(B-Ro,U-Do)
 rows.append(dict(delta_over_h=a,new_D=D,new_R=R,new_H=H,new_H_percent=100*H/B,old_D=Do,old_interval_R=Ro,old_H=Ho,old_H_percent=100*Ho/B,gain=H-Ho,new_nonnegative_feasible=bool(H>=0),old_nonnegative_feasible=bool(Ho>=0),new_binding='budget' if B-R<U-D else 'span',old_binding='budget' if B-Ro<U-Do else 'span'))
(OUT/'benchmark_comparison.json').write_text(json.dumps(rows,indent=2))
for x in rows:print(x)
