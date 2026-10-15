"""Numerical interval-wide Gaussian bounds; no Gaussian iff claim.
Matched p=h=1, sigma=h/sqrt(3), q, delta, B and U. Both laws permit every
nondecreasing globally p-Lipschitz bounded payment under a ceiling at every centre.
The feasible side uses piecewise-constant slopes, with an analytic second-derivative
interpolation guard covering the entire centre interval. The upper side relaxes arbitrary
slope measures, including BOTH infinite tails, using density suprema and survival
infima on each report segment. A residual-corrected weak-duality bound is reported.
Floating coefficient/rounding safety allowances are explicit; these are numerical
certificates to the stated tolerances, not exact rational/interval-arithmetic proofs.
"""
from pathlib import Path
import csv,json,math
import numpy as np
from scipy.optimize import linprog
from scipy.special import ndtr
from scipy.sparse import csr_matrix
from ramp_interval_checks import global_frontier
OUT=Path(__file__).resolve().parent
SIGMA=1/math.sqrt(3)
PDF_CONSTANT=1/(SIGMA*math.sqrt(2*math.pi))
KERNEL_LIPSCHITZ=math.sqrt(2/math.pi)/SIGMA # p * integral |f'|
KERNEL_SECOND_DERIVATIVE=4*math.exp(-.5)/(math.sqrt(2*math.pi)*SIGMA**2) # p*integral |f''|
COEFFICIENT_SAFETY=2e-13
CERTIFICATE_SAFETY=2e-9 # normalized money p*h or normalized slope p, as stated
SOLVER_OPTIONS={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,'ipm_optimality_tolerance':1e-10}

def pdf(z):return np.exp(-.5*(np.asarray(z)/SIGMA)**2)*PDF_CONSTANT

def primitive_survival(y,e):
    z=(e-y)/SIGMA
    return (y-e)*ndtr(z)-SIGMA*np.exp(-.5*z*z)/math.sqrt(2*math.pi)

def report_grid(alpha,n):return np.linspace(-alpha-8*SIGMA,alpha+8*SIGMA,n+1)

def feasible_lower(alpha,q,B,U,n=600,max_mesh=.002):
    if not 0<=q<1:raise ValueError("Gaussian bounded payment requires q<p; q=p>0 is unattainable.")
    knots=report_grid(alpha,n);left=knots[:-1];right=knots[1:];dx=right-left
    count=max(1,int(math.ceil(2*alpha/max_mesh)))
    centers=np.linspace(-alpha,alpha,count+1) if alpha else np.array([0.])
    mesh=float(np.max(np.diff(centers))) if len(centers)>1 else 0.
    guard=KERNEL_SECOND_DERIVATIVE*mesh*mesh/8+20*COEFFICIENT_SAFETY
    kernel=ndtr((right[None,:]-centers[:,None])/SIGMA)-ndtr((left[None,:]-centers[:,None])/SIGMA)
    budget=primitive_survival(right,alpha)-primitive_survival(left,alpha)
    # A variable floor C is valid because test budgets ensure nonnegative C.
    A=np.vstack([np.r_[dx,1.],np.r_[budget,1.],np.c_[-kernel,np.zeros(len(centers))]])
    rhs=np.r_[U,B,np.full(len(centers),-q-guard-2e-8)]
    res=linprog(np.r_[np.zeros(n),-1.],A_ub=csr_matrix(A),b_ub=rhs,bounds=[(0,1)]*n+[(0,min(B,U))],method='highs',options=SOLVER_OPTIONS)
    if not res.success:raise AssertionError(res.message)
    # Explicit box repair never increases the slopes or available resources.
    g=np.clip(res.x[:-1],0,1)
    span=float(dx@g);reserve=float(budget@g)
    C=min(float(res.x[-1]),B-reserve,U-span)-CERTIFICATE_SAFETY
    computed_min=float(np.min(kernel@g))
    continuum_min_lower=computed_min-KERNEL_SECOND_DERIVATIVE*mesh*mesh/8-CERTIFICATE_SAFETY
    assert continuum_min_lower>=q,(alpha,q,continuum_min_lower)
    assert C>=0 and C+span<=U and C+reserve<=B
    return dict(C=C,span=span,reserve=reserve,min_grid_marginal=computed_min,continuum_marginal_lower=continuum_min_lower,grid_mesh=mesh,interpolation_guard=guard,second_derivative_bound=KERNEL_SECOND_DERIVATIVE,centers=len(centers),segments=n,primal_violation=float(max(0,np.max(A@res.x-rhs))),support_left=float(knots[0]),support_right=float(knots[-1]))

def unrestricted_upper(alpha,q,B,U,n=600,centers_n=101):
    knots=report_grid(alpha,n);left=knots[:-1];right=knots[1:];dx=right-left
    centers=np.linspace(-alpha,alpha,centers_n) if alpha else np.array([0.])
    closest=np.maximum(left[None,:],np.minimum(right[None,:],centers[:,None]))
    kernel=pdf(closest-centers[:,None])+COEFFICIENT_SAFETY
    tails=np.c_[pdf(left[0]-centers),pdf(right[-1]-centers)]+COEFFICIENT_SAFETY
    kernel=np.c_[tails[:,0],kernel,tails[:,1]]
    # Infimum of Gaussian upper-centre survival over each segment. Right tail
    # infimum is zero, left tail infimum is survival at the finite boundary.
    w=np.r_[ndtr((alpha-left[0])/SIGMA),ndtr((alpha-right)/SIGMA),0.]
    w=np.maximum(0,w-COEFFICIENT_SAFETY)
    count=n+2
    A=np.vstack([np.r_[np.ones(count),1.],np.r_[w,1.],np.c_[-kernel,np.zeros(len(centers))]])
    rhs=np.r_[U,B,np.full(len(centers),-q+COEFFICIENT_SAFETY)]
    ub=np.r_[U,np.minimum(dx,U),U,min(B,U)]
    obj=np.r_[np.zeros(count),-1.]
    res=linprog(obj,A_ub=csr_matrix(A),b_ub=rhs,bounds=list(zip(np.zeros(len(ub)),ub)),method='highs',options=SOLVER_OPTIONS)
    if not res.success:raise AssertionError(res.message)
    # For ANY y<=0, weak duality over boxes gives:
    # min c*x >= rhs*y + sum_j min(0,c_j-A_j*y)*ub_j.
    # This repairs dual reduced-cost residuals rather than trusting solver status.
    y=np.minimum(np.asarray(res.ineqlin.marginals),0.)
    reduced=obj-A.T@y
    correction=float(np.minimum(reduced,0.)@ub)
    dual_lower=float(rhs@y+correction)
    upper=-dual_lower+CERTIFICATE_SAFETY
    assert upper+1e-10>=res.x[-1]
    return dict(C_upper=upper,LP_C=float(res.x[-1]),weak_duality_correction=correction,dual_residual_bound=upper-float(res.x[-1]),segments=n,necessary_centers=len(centers),left_tail_mass=float(res.x[0]),right_tail_mass=float(res.x[-2]),primal_violation=float(max(0,np.max(A@res.x-rhs))))

def run():
    rows=[];records=[]
    # All pairs use the same two authority ceilings and no financing charges.
    # These are dimensionless scenarios, not empirical Gaussian calibration.
    B=8.;U=9.
    for q in [.25,.5,.75]:
      for alpha in [0.,.25,.5,1.,2.]:
        D,R=global_frontier(1,1,q,alpha);uniform=min(B-R,U-D)
        lower=feasible_lower(alpha,q,B,U);upper=unrestricted_upper(alpha,q,B,U)
        gap=upper['C_upper']-lower['C']
        assert gap>=-CERTIFICATE_SAFETY and gap<.075,(q,alpha,gap)
        difference_lo=lower['C']-uniform;difference_hi=upper['C_upper']-uniform
        order='Gaussian higher' if difference_lo>0 else ('Uniform higher' if difference_hi<0 else 'Not resolved by bounds')
        row=dict(q_over_p=q,delta_over_h=alpha,B_over_ph=B,U_over_ph=U,uniform_exact_capacity_over_ph=uniform,gaussian_capacity_lower_over_ph=lower['C'],gaussian_capacity_upper_over_ph=upper['C_upper'],gaussian_bound_gap_over_ph=gap,gaussian_minus_uniform_lower_over_ph=difference_lo,gaussian_minus_uniform_upper_over_ph=difference_hi,certified_order=order,gaussian_continuum_marginal_lower_over_p=lower['continuum_marginal_lower'])
        rows.append(row);records.append(dict(**row,feasible_certificate=lower,upper_certificate=upper))
        print(f'q/p={q:.2f}, delta/h={alpha:.2f}: Uniform={uniform:.8f}, Gaussian=[{lower["C"]:.8f},{upper["C_upper"]:.8f}], {order}',flush=True)
    with (OUT/'gaussian_interval_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
    summary=dict(cases=len(rows),variance=1/3,reporting_sigma=SIGMA,normalized_money_safety=CERTIFICATE_SAFETY,normalized_slope_safety=CERTIFICATE_SAFETY,coefficient_safety=COEFFICIENT_SAFETY,max_bound_gap=max(x['gaussian_bound_gap_over_ph'] for x in rows),orders={label:sum(x['certified_order']==label for x in rows) for label in ['Gaussian higher','Uniform higher','Not resolved by bounds']},solver_options=SOLVER_OPTIONS,finite_grid_is_not_continuum_certificate=True,tail_mass_relaxation_included=True,no_general_gaussian_iff_claim=True,support_margin_standard_deviations=8,normal_two_sided_probability_beyond_8_sigma=2*float(ndtr(-8)),tail_marginal_relaxation_bound=2*U*float(pdf(8*SIGMA)),gaussian_q_equal_p_unattainable=True)
    (OUT/'gaussian_interval_certificates.json').write_text(json.dumps(dict(summary=summary,records=records),indent=2))
    return summary
if __name__=='__main__':print(json.dumps(run(),indent=2))
