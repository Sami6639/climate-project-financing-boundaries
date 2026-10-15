"""Independent continuum-feasible piecewise-linear LP validator.
Units: t=(report-mu)/h; g=T'/p; money in p*h.
Every sliding-window breakpoint is included, so each returned piecewise-linear
contract satisfies the continuum constraints, not only an arbitrary sample grid.
A separate optimistic segment-mass LP bounds the unrestricted optimum above.
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix,vstack
import json, pathlib
OUT=pathlib.Path(__file__).parent

def solve(a,r,b,u,n=400,relaxed=False):
    knots=np.linspace(-1-a,1+a,n+1); left=knots[:-1]; right=knots[1:]; dx=right-left
    # Exact kink locations of moving integrals for piecewise constant g.
    centers=np.unique(np.clip(np.r_[-a,a,knots-1,knots+1],-a,a))
    overlap=np.maximum(0,np.minimum(right[None,:],centers[:,None]+1)-np.maximum(left[None,:],centers[:,None]-1))
    # Survival probability is one before the uppercenter's window.
    boundary=a-1; end=a+1
    def primitive(t):
        return np.where(t<=boundary,t,boundary+(end*(t-boundary)-(t*t-boundary*boundary)/2)/2)
    w=(primitive(right)-primitive(left))/dx
    if relaxed:
        # For arbitrary bounded slopes the full segment mass can optimistically
        # serve every intersecting window. Budget coefficient is min weight.
        overlap=np.where(overlap>1e-12,dx[None,:],0.)
        w=np.minimum(1.,(a+1-right)/2)
    A=np.vstack([np.r_[dx,1.],np.r_[dx*w,1.],np.c_[-overlap,np.zeros(len(centers))]])
    rhs=np.r_[u,b,np.full(len(centers),-2*r)]
    res=linprog(np.r_[np.zeros(n),-1.],A_ub=csr_matrix(A),b_ub=rhs,bounds=[(0,1)]*n+[(None,None)],method='highs')
    if not res.success: return dict(success=False,message=res.message)
    g=res.x[:-1]; c=res.x[-1]
    windows=overlap@g
    dual=float(rhs@res.ineqlin.marginals+np.ones(n)@res.upper.marginals[:-1])
    return dict(success=True,a=float(a),r=float(r),n=n,c=float(c),mass=float(dx@g),reserve=float((dx*w)@g),window_min=float(min(windows)),max_primal_violation=float(max(0,np.max(A@res.x-rhs))),dual_objective=dual,primal_dual_gap=float(abs(res.fun-dual)),knots=knots.tolist(),g=g.tolist())

def run():
    h=np.sqrt(3)*.1*316912.; p=25.; B=7925000.; U=1.1*B; scale=p*h
    records=[]
    for a in [0,.1,.25,.5,.75,1.]:
      for n in [200,400,800]:
        lo=solve(a,.75,B/scale,U/scale,n); hi=solve(a,.75,B/scale,U/scale,n,True)
        rec={k:v for k,v in lo.items() if k not in ('knots','g')}
        rec.update(advance_eur=lo['c']*scale,advance_percent=lo['c']*scale/B*100,upper_advance_eur=hi['c']*scale,error_bound_eur=(hi['c']-lo['c'])*scale)
        records.append(rec)
        if n==800: (OUT/f'contract_a{a}.json').write_text(json.dumps(lo,indent=2))
    (OUT/'convergence.json').write_text(json.dumps({'h':h,'money_scale':scale,'records':records},indent=2))
    print(json.dumps(records,indent=2))
if __name__=='__main__': run()
