"""Stage 6 independent one-at-a-time sensitivity, fixed physical interval unless explicitly changed.
Report samples are hypothetical; changes are not estimates.
"""
import math, csv, json, random
from pathlib import Path
from scipy.optimize import linprog
P=Path(__file__).resolve().parent/'outputs';P.mkdir(exist_ok=True)
B=7925000.0; U=1.1*B; mu=316912.0; p=25.0; q=18.75
h=math.sqrt(3)*.1*mu; delta=.25*h; Hreq=.9*B

def calculate(v):
    B,U,h,delta,p,q,kappa,V=[v[k] for k in ('B','U','h','delta','p','q','kappa','V')]
    assert h>0 and p>0 and q>0 and 0<=q<=p and delta>=0 and B>=0 and U>=0 and kappa>=0 and V>=0
    W=2*h; L=2*delta
    n=math.floor((L+1e-10*W)/W)
    ell=max(0.0,L-n*W)
    R=h*q*q/p+L*q
    D=(n+1)*W*q + max(0,W*q-p*(W-ell))
    capacity_budget=(B-R-V)/(1+kappa)
    capacity_cap=(U-D-V)/(1+kappa)
    capacity=min(capacity_budget, capacity_cap)
    deficit_budget=(1+kappa)*Hreq+V+R-B
    deficit_exposure=(1+kappa)*Hreq+V+D-U
    binding=('budget' if capacity_budget<capacity_cap-1e-7 else 'exposure' if capacity_cap<capacity_budget-1e-7 else 'tie')
    return dict(R=R,D=D,cap_B=capacity_budget,cap_U=capacity_cap,Hmax=capacity,shortfall=max(0,Hreq-capacity),binding=binding,expected_budget_deficit=deficit_budget,exposure_deficit=deficit_exposure,feasible=capacity>=Hreq-1e-7)

base=dict(B=B,U=U,h=h,delta=delta,p=p,q=q,kappa=0.,V=0.)
choices=[
('Base: δ/h₀=0.25',{}),
('Reporting SD: 7.5%',{'h':math.sqrt(3)*.075*mu}),
('Reporting SD: 12.5%',{'h':math.sqrt(3)*.125*mu}),
('Target q: €12/t',{'q':12.}),
('Target q: €15/t',{'q':15.}),
('Target q: €20/t',{'q':20.}),
('Rate cap p: €30/t',{'p':30.}),
('Mean budget: 110% B₀',{'B':1.10*B}),
('Hard cap: 120% B₀',{'U':1.20*B}),
('Advance loading: 2%',{'kappa':.02}),
('Verification cost: €25,000',{'V':25000.}),
('Interval δ: 0',{'delta':0.}),
('Interval δ: 0.5h₀',{'delta':.5*h}),
('Interval δ: 1.25h₀',{'delta':1.25*h}),
]
res=[]
for title,changed in choices:
    case=base.copy();case.update(changed)
    d=calculate(case)
    res.append(dict(scenario=title,**d,delta_over_current_h=case['delta']/case['h'],h_tonnes=case['h'],p=case['p'],q=case['q'],B=case['B'],U=case['U'],kappa=case['kappa'],V=case['V']))

# independent LP for financing given the two sharp resource requirements
rng=random.Random(20261026)
maxerr=0
for k in range(300):
    c=dict(B=rng.uniform(8.,30.),U=rng.uniform(5.,30.),h=rng.uniform(.1,2.),delta=rng.uniform(0,3.),p=rng.uniform(.5,4.),q=0,kappa=rng.uniform(0,1.),V=rng.uniform(0,.5))
    c['q']=rng.uniform(.01,c['p'])
    z=calculate(c)
    # LP over (H,c), where c is guaranteed cost-floor baseline. Max H subject to c>=(1+k)*H+V,
    # c+R<=B, c+D<=U; simply ensures joint feasibility of R and D construction.
    lp=linprog(c=[-1,0],A_ub=[[1+c['kappa'],-1],[0,1],[0,1]],b_ub=[-c['V'],c['B']-z['R'],c['U']-z['D']],bounds=[(0,None),(0,None)],method='highs')
    if z['Hmax'] < -1e-7:
        # all inputs nonnegative, negative bound => no feasible H; no LP certificate
        assert not lp.success, (k,z,lp.message)
    else:
        assert lp.success,(k,z,lp.message)
        err=abs(z['Hmax']-lp.x[0])
        maxerr=max(maxerr,err)
        assert err<1e-6,(k,err)

# exact check with known base results and two reference variations
assert abs(res[0]['Hmax']-6638497.612043122)<1e-6
assert abs(res[0]['shortfall']-494002.38795687817)<1e-6
assert abs(res[-3]['Hmax']-6659096.179268995)<1e-6
assert abs(res[-2]['Hmax']-5972961.5723586595)<1e-6
assert abs(res[-1]['Hmax']-4580093.791312116)<1e-6

with open(P/'one_factor_sensitivity.csv','w',newline='',encoding='utf8') as f:
  writer=csv.DictWriter(f,fieldnames=res[0].keys());writer.writeheader();writer.writerows(res)
(P/'sensitivity_audit.json').write_text(json.dumps({'status':'PASS','case_count':len(res),'independent_lp_cases':300,'maximum_lp_gap':maxerr,'physical_halfwidth_held_constant_except_named_delta_scenarios':delta,'one_factor_policy':'vary exactly one of B,U,h,delta,p,q,kappa,V; keep initial project cash need fixed','interpretation':'scenario comparisons only; uniform fixed-width reporting model, not empirically estimated; no new Gaussian theorem'},indent=2),encoding='utf8')
print('PASS',len(res),'scenarios; 300 LP checks; max gap',maxerr)
for row in res:
 print(f"{row['scenario']:30s} Hmax {row['Hmax']:12,.2f} shortfall {row['shortfall']:11,.2f} bind {row['binding']:8s} feasible {row['feasible']}")
