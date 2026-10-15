"""Complete fixed-cash financing price regions and independent LP validation.
Run: python cost/validate_price_regions.py
All monetary quantities are in one common currency; prices are hypothetical.
No plotting, publication or external network access. Failure raises/nonzero exit.
"""
from pathlib import Path
import csv, importlib.util, json, math, platform
import numpy as np
import scipy
from scipy.optimize import linprog

OUT = Path(__file__).resolve().parent
SEED = 20261014
LP_OPTIONS = {'primal_feasibility_tolerance': 1e-9,
              'dual_feasibility_tolerance': 1e-9,
              'ipm_optimality_tolerance': 1e-10}
OBJ_TOL = 1e-8       # objective divided by the common cash scale
FEAS_TOL = 1e-8      # primal residual divided by the common cash scale
FACE_TOL = 1e-7      # optimal-s endpoint divided by the common cash scale
DERIV_TOL = 2e-5     # relative derivative error with denominator max(1, |derivative|)


def price_region(A, E, H, alpha, cb, cu, cs):
    """Return the COMPLETE optimal s interval, not just one breakpoint.
    Price equalities are mathematical equalities, never rounded price bands.
    """
    assert H >= 0 and alpha > 0 and min(cb, cu, cs) > 0
    if H == 0:
        return 0., 0., 'zero_cash_requirement'
    if max(A,E) <= 0:
        return 0., 0., 'no_positive_deficit'
    if min(A,E) <= 0:
        d, c = (A,cb) if A > 0 else (E,cu)
        t = min(H,d/alpha)
        threshold = alpha*c
        if cs > threshold: return 0.,0.,'one_gate_authority'
        if cs == threshold: return 0.,t,'one_gate_price_tie'
        return t,t,'one_gate_substitution'
    if A == E:
        t = min(H,A/alpha)
        threshold = alpha*(cb+cu)
        if cs > threshold: return 0.,0.,'equal_gates_authority'
        if cs == threshold: return 0.,t,'equal_gates_price_tie'
        return t,t,'equal_gates_substitution'
    small, large = min(A,E), max(A,E)
    c_large = cb if A > E else cu
    t1,t2 = min(H,small/alpha), min(H,large/alpha)
    upper, lower = alpha*(cb+cu), alpha*c_large
    if cs > upper: return 0.,0.,'two_gates_authority'
    if cs == upper: return 0.,t1,'upper_price_tie'
    if cs > lower: return t1,t1,'mixed_or_cash_cap'
    if cs == lower: return t1,t2,'lower_price_tie'
    return t2,t2,'substitution_or_cash_cap'


def intervention(A,E,H,alpha,cb,cu,cs,s):
    b,u = max(0.,A-alpha*s), max(0.,E-alpha*s)
    cash = H-s+s
    # Binding flags are diagnostic only; price region selection above is exact.
    tol = 2e-9*max(1.,abs(A),abs(E),H)
    return dict(budget_upgrade=b,exposure_upgrade=u,external_nonrepayable_funds=s,
                donor_advance=H-s,total_initial_project_cash=cash,
                residual_cash_gap=max(0.,H-cash),
                residual_budget_gate_violation=max(0.,A-b-alpha*s),
                residual_exposure_gate_violation=max(0.,E-u-alpha*s),
                budget_slack=max(0.,alpha*s-A),exposure_slack=max(0.,alpha*s-E),
                budget_gate_binding=abs(b+alpha*s-A)<=tol,
                exposure_gate_binding=abs(u+alpha*s-E)<=tol,
                external_lower_bound_binding=s<=tol,
                external_upper_bound_binding=abs(s-H)<=tol,
                weighted_resource_cost=cb*b+cu*u+cs*s)


def lp_problem(A,E,H,alpha,cb,cu,cs):
    scale=max(1.,abs(A),abs(E),H)
    weights=np.array([cb,cu,cs],float)
    mat=np.array([[-1.,0.,-alpha],[0.,-1.,-alpha]])
    rhs=np.array([-A,-E])/scale
    bounds=[(0,None),(0,None),(0,H/scale)]
    result=linprog(weights,A_ub=mat,b_ub=rhs,bounds=bounds,
                   method='highs',options=LP_OPTIONS)
    assert result.success,result.message
    return result,scale,weights,mat,rhs,bounds


def check_case(label,A,E,H,alpha,cb,cu,cs,check_face=False,expected=None):
    lo,hi,regime=price_region(A,E,H,alpha,cb,cu,cs)
    r,scale,c,mat,rhs,bounds=lp_problem(A,E,H,alpha,cb,cu,cs)
    row=intervention(A,E,H,alpha,cb,cu,cs,lo)
    error=abs(row['weighted_resource_cost']/scale-r.fun)
    x=np.array([row['budget_upgrade'],row['exposure_upgrade'],lo])/scale
    primal=max(0.,float(np.max(mat@x-rhs)),float(-np.min(x)),lo/scale-H/scale)
    assert error<=OBJ_TOL,(label,'objective',error)
    assert primal<=FEAS_TOL,(label,'feasibility',primal)
    # Every point in the reported interval must be feasible and optimal.
    for s in [hi,(lo+hi)/2]:
        test=intervention(A,E,H,alpha,cb,cu,cs,s)
        assert abs(test['weighted_resource_cost']/scale-r.fun)<=OBJ_TOL,(label,'interval')
    face_error=0.
    if check_face:
        # Independent optimal-face LPs use no analytic breakpoint constraints.
        eq=np.array([c]); val=np.array([r.fun])
        min_s=linprog([0,0,1],A_ub=mat,b_ub=rhs,A_eq=eq,b_eq=val,bounds=bounds,method='highs',options=LP_OPTIONS)
        max_s=linprog([0,0,-1],A_ub=mat,b_ub=rhs,A_eq=eq,b_eq=val,bounds=bounds,method='highs',options=LP_OPTIONS)
        assert min_s.success and max_s.success,(label,'optimal-face LP failed')
        face_error=max(abs(min_s.x[2]-lo/scale),abs(max_s.x[2]-hi/scale))
        assert face_error<=FACE_TOL,(label,'optimal face',face_error,(lo,hi),min_s.x,max_s.x)
    if expected is not None:
        assert np.allclose([lo,hi],expected,rtol=0,atol=1e-12),(label,'expected s interval',(lo,hi),expected)
    return dict(label=label,A=A,E=E,Hreq=H,alpha=alpha,cB=cb,cU=cu,cS=cs,
                regime=regime,optimal_s_min=lo,optimal_s_max=hi,expected_result=f'LP minimum and maximum optimal s must equal [{lo:.17g}, {hi:.17g}]',**row,
                normalized_objective_error=error,normalized_primal_violation=primal,
                normalized_optimal_face_error=face_error,
                objective_pass_tolerance=OBJ_TOL,feasibility_pass_tolerance=FEAS_TOL,
                face_pass_tolerance=FACE_TOL if check_face else None,
                unit='currency; objective errors divided by max(1,|A|,|E|,Hreq)',passed=True)


def reserves(h,p,q,delta):
    assert h>0 and p>=0 and 0<=q<=p and delta>=0
    if q==0: return 0.,0.
    n=math.floor(delta/h)
    remainder=delta-n*h
    R=h*q*q/p+2*delta*q
    D=(n+1)*2*h*q+max(0.,2*h*q-2*p*(h-remainder))
    return R,D


def model_case(label,H=7132500.,B=7925000.,U=8717500.,kappa=0.,V=0.,
               h=math.sqrt(3)*.1*316912,p=25.,q=18.75,delta_ratio=.25,
               cb=.25,cu=.125,cs=1.,check_face=False):
    delta=delta_ratio*h
    R,D=reserves(h,p,q,delta)
    alpha=1+kappa
    A,E=alpha*H+V+R-B,alpha*H+V+D-U
    row=check_case(label,A,E,H,alpha,cb,cu,cs,check_face)
    row.update(B=B,U=U,kappa=kappa,V=V,h=h,p=p,q=q,delta=delta,delta_over_h=delta_ratio,R=R,D=D)
    return row


def purpose_cases():
    rows=[]
    # Each deficit shape targets a different clipping, degeneracy or active-set risk.
    shapes=[('both_negative',-3.,-1.,10.),('zero_and_negative',0.,-1.,10.),
            ('both_zero',0.,0.,10.),('budget_only',3.,-1.,10.),
            ('exposure_only',-1.,3.,10.),('zero_budget_positive_exposure',0.,3.,10.),
            ('positive_budget_zero_exposure',3.,0.,10.),
            ('equal_positive',3.,3.,10.),('budget_larger',6.,2.,10.),
            ('exposure_larger',2.,6.,10.),('cap_before_both_breaks',12.,16.,10.),
            ('cap_between_breaks',4.,16.,10.),('cap_at_first_break',10.,16.,10.),
            ('cap_at_last_break',4.,10.,10.),('cap_equal_breaks',10.,10.,10.),
            ('one_gate_beyond_cap',16.,-1.,10.),('zero_cash',3.,7.,0.)]
    # Distinct prices expose the reversed-ordering error; binary fractions ensure exact ties.
    cb,cu,alpha=.25,.75,1.
    for shape,A,E,H in shapes:
        for pname,cs in [('cheap',.125),('budget_price',.25),('middle',.5),
                         ('exposure_price',.75),('sum_price',1.),('expensive',1.25)]:
            rows.append(check_case(shape+'__'+pname,A,E,H,alpha,cb,cu,cs,True))
    # Predeclared expected intervals independently expose missing ties/cap truncation.
    known=[('known_upper_tie',6,2,10,1,.25,.75,1,(0,2)),
           ('known_lower_tie',6,2,10,1,.25,.75,.25,(2,6)),
           ('known_reverse_lower_tie',2,6,10,1,.25,.75,.75,(2,6)),
           ('known_equal_tie',3,3,10,1,.25,.75,1,(0,3)),
           ('known_single_tie',-2,3,10,1,.25,.75,.75,(0,3)),
           ('known_cap_tie',12,16,10,1,.25,.75,1,(0,10)),
           ('known_collapsed_lower_tie',12,16,10,1,.25,.75,.75,(10,10)),
           ('known_partial_lower_tie',4,16,10,1,.25,.75,.75,(4,10)),
           ('known_loading_upper_tie',6,2,10,2,.25,.75,2,(0,1)),
           ('known_loading_lower_tie',6,2,10,2,.25,.75,.5,(1,3))]
    for label,A,E,H,alpha,cb,cu,cs,exp in known:
        rows.append(check_case(label,A,E,H,alpha,cb,cu,cs,True,exp))
    # Two sides of each strict price threshold; equality is not a rounding band.
    for t in [.25,1.]:
        for side in [-1,1]:
            rows.append(check_case(f'price_threshold_{t}_side_{side}',6,2,10,1,.25,.75,t+side*1e-5,True))
    # Very close signed-deficit and equal-deficit boundaries.
    for eps in [-1e-5,0.,1e-5]:
        rows.append(check_case(f'signed_boundary_{eps}',eps,3,10,1,.25,.75,.9,True))
        rows.append(check_case(f'ordering_boundary_{eps}',3+eps,3,10,1,.25,.75,.5,True))
    # Specific interval kinks with h positive and fixed, plus zero/full targets.
    for ratio in [0,.25-1e-8,.25,.25+1e-8,1-1e-8,1,1+1e-8,
                  1.25-1e-8,1.25,1.25+1e-8,2-1e-8,2,2+1e-8]:
        for name,cb,cu,cs in [('authority',.25,.125,1),('mixed',.25,1.5,1),('external',1.5,.25,1)]:
            rows.append(model_case(f'interval_kink_{ratio}_{name}',delta_ratio=ratio,cb=cb,cu=cu,cs=cs))
    for q in [0,25]:
        for ratio in [0,.25,1,1+1e-8,2]:
            rows.append(model_case(f'target_endpoint_q{q}_r{ratio}',q=q,delta_ratio=ratio,kappa=.05,V=100000))
    rows.append(model_case('zero_rate_cap_zero_target',q=0,p=0,kappa=.05,V=100000))
    rows.append(model_case('small_positive_width',h=1e-10,delta_ratio=.25,kappa=.25,V=50000,cb=.25,cu=1.5))
    # At s=H both donor gates may still need authority; external cash is not permission to drop the target.
    rows.append(model_case('all_external_still_both_authority',B=.1,U=.2,H=2,kappa=.2,V=3,h=1,p=2,q=1,delta_ratio=.5,cb=2,cu=3,cs=.1))
    return rows


def prior_random_cases():
    spec=importlib.util.spec_from_file_location('legacy_cost',OUT/'legacy'/'cost_comparison.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    assert mod.validate_random()==1000
    rng=np.random.default_rng(SEED)
    rows=[]
    for i in range(1000):
        H=float(rng.uniform(1,1e7));alpha=float(rng.uniform(1,1.4));V=float(rng.uniform(0,1e6))
        R=float(rng.uniform(0,2e6));D=float(rng.uniform(0,4e6))
        B=float(rng.uniform(0,1.5e7));U=float(rng.uniform(0,1.5e7));cb,cu,cs=rng.uniform(.01,2,3)
        A,E=alpha*H+V+R-B,alpha*H+V+D-U
        rows.append(check_case(f'original_seed_case_{i:04d}',A,E,H,alpha,float(cb),float(cu),float(cs)))
    return rows


def sensitivities():
    rows=[]
    price_sets=[('P1',.25,.125,1),('P2',.25,1.5,1),('P3',1.5,.25,1)]
    for name,cb,cu,cs in price_sets:
        for par,values in [('kappa',[0,.02,.05,.1,.2]),('V',[0,25000,100000,250000,1000000]),
                           ('delta_ratio',[0,.1,.25,.5,.75,1,1.25,2])]:
            for value in values:
                row=model_case(f'{name}_{par}_{value}',cb=cb,cu=cu,cs=cs,**{par:value})
                row.update(sensitivity_parameter=par,sensitivity_value=value,price_scenario=name)
                rows.append(row)
    # Vary one price at a time including each baseline threshold.
    for par,values in [('cb',[.125,.25,.5,.75,1,1.25,1.5]),
                       ('cu',[.125,.5,.75,1,1.5]),('cs',[.125,.25,.375,.5,1,1.75,2])]:
        for value in values:
            row=model_case(f'price_{par}_{value}',**{par:value},check_face=True)
            row.update(sensitivity_parameter=par,sensitivity_value=value,price_scenario='one_at_a_time_from_authority')
            rows.append(row)
    return rows


def derivative_checks():
    rows=[]
    # Stay away from interval/price/active-set kinks; dual-envelope derivatives are exact there.
    for label,cb,cu,cs in [('P1',.25,.125,1),('P2',.25,1.5,1),('P3',1.5,.25,1)]:
      for dr in [.2,.3,.7,1.1]:
        base=dict(kappa=.05,V=25000.,delta_ratio=dr,cb=cb,cu=cu,cs=cs)
        row=model_case('base',**base)
        res,scale,weights,mat,rhs,bounds=lp_problem(row['A'],row['E'],row['Hreq'],row['alpha'],cb,cu,cs)
        lb,lu=-res.ineqlin.marginals
        frac=dr-math.floor(dr)
        dD=0 if frac<1-row['q']/row['p'] else 2*row['p']
        expected={'V':lb+lu,'kappa':(row['Hreq']-row['external_nonrepayable_funds'])*(lb+lu),
                  'delta_ratio':row['h']*(2*row['q']*lb+dD*lu),
                  'cb':row['budget_upgrade'],'cu':row['exposure_upgrade'],'cs':row['external_nonrepayable_funds']}
        for par in expected:
            step=1e-4 if par=='V' else 1e-6
            plus=base.copy();minus=base.copy();plus[par]+=step;minus[par]-=step
            vplus=model_case('plus',**plus)['weighted_resource_cost']
            vminus=model_case('minus',**minus)['weighted_resource_cost']
            measured=(vplus-vminus)/(2*step)
            err=abs(measured-expected[par])/max(1.,abs(expected[par]))
            assert err<=DERIV_TOL,(label,dr,par,err,measured,expected[par])
            rows.append(dict(price_scenario=label,delta_over_h=dr,parameter=par,
                             exact_envelope_derivative=expected[par],finite_difference_derivative=measured,
                             relative_error=err,tolerance=DERIV_TOL,step=step,passed=True))
    return rows


def dump_csv(name,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def run():
    legacy=prior_random_cases();purpose=purpose_cases();sens=sensitivities();deriv=derivative_checks()
    main=[model_case(name,cb=cb,cu=cu,cs=cs,check_face=True) for name,cb,cu,cs in
          [('authority',.25,.125,1),('mixed',.25,1.5,1),('external',1.5,.25,1)]]
    dump_csv('original_1000_price_region_checks.csv',legacy)
    dump_csv('purposeful_price_region_checks.csv',purpose)
    dump_csv('cost_sensitivity.csv',sens)
    dump_csv('cost_derivative_checks.csv',deriv)
    dump_csv('main_price_examples.csv',main)
    phase=[]
    # Coordinates x=alpha*cB/cS, y=alpha*cU/cS. This scenario has A>E>0 and no cash cap truncation.
    for x in np.linspace(.05,1.75,35):
      for y in np.linspace(.05,1.75,35):
        r=model_case('phase',cb=float(x),cu=float(y),cs=1.)
        phase.append({k:r[k] for k in ['cB','cU','cS','regime','optimal_s_min','optimal_s_max','budget_upgrade','exposure_upgrade','external_nonrepayable_funds','residual_cash_gap']})
    dump_csv('price_phase_data.csv',phase)
    all_rows=legacy+purpose+sens+main+phase
    # Phase rows are a display grid, not counted as independent scientific tests.
    checked=legacy+purpose+sens+main
    report=dict(status='PASS',seed=SEED,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        solver='scipy.optimize.linprog(method="highs")',solver_options=LP_OPTIONS,
        legacy_original_random_checks_retained=1000,original_same_cases_new_price_rule_checks=len(legacy),
        purposeful_edge_cases=len(purpose),sensitivity_rows=len(sens),derivative_checks=len(deriv),main_rows=len(main),
        normalized_objective_tolerance=OBJ_TOL,normalized_feasibility_tolerance=FEAS_TOL,
        normalized_optimal_face_tolerance=FACE_TOL,derivative_relative_tolerance=DERIV_TOL,
        max_normalized_objective_error=max(r['normalized_objective_error'] for r in checked),
        max_normalized_primal_violation=max(r['normalized_primal_violation'] for r in checked),
        max_normalized_optimal_face_error=max(r['normalized_optimal_face_error'] for r in checked),
        max_derivative_relative_error=max(r['relative_error'] for r in deriv),
        interpretation='Conditional hypothetical scarcity valuations; unchanged project cash requirement; no observed expenditure or welfare claim.',
        main_rows_detail=main)
    (OUT/'price_region_validation_report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k!='main_rows_detail'},indent=2))

if __name__=='__main__': run()
