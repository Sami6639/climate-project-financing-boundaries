"""Independent finite-candidate equality-regime and manuscript-figure audit."""
from pathlib import Path
import csv,json,math,sys
from ramp_interval_checks import global_frontier,best_ramp,fixed_ramp
OUT=Path(__file__).resolve().parent;THEORY=OUT.parent/'theory'
sys.path.insert(0,str(THEORY))
from single_ramp_optimum import solve_single_ramp,zero_gap_region

def run():
    count=0;max_equal_error=0.;max_formula_error=0.;min_strict_gap=float('inf')
    for q in [0.,.01,.1,.25,.5,.75,.9,.99,1.]:
      boundary=1-q
      widths=sorted(set([0.,.1,.25,.5,.75,1.,1.5,2.,2.5,3.,boundary]+[max(0.,boundary+x) for x in [-1e-4,1e-4]]+[k+x for k in [1.,2.,3.] for x in [-1e-4,1e-4]]))
      for alpha in widths:
       for budget_gap in sorted(set([-.2,0.,.05,.5,1.2,q*q-1e-4,q*q,q*q+1e-4,q-1e-4,q,q+1e-4])):
        B=20.;U=B+budget_gap;D,R=global_frontier(1.,1.,q,alpha);C=min(B-R,U-D);opt=best_ramp(1.,1.,q,alpha,B,U);gap=C-opt['c']
        other=solve_single_ramp(h=1.,delta=alpha,p=1.,q=q,B=B,U=U)
        implementation_error=abs(other['c']-opt['c']);max_formula_error=max(max_formula_error,implementation_error)
        assert implementation_error<2e-10
        predicted=(q in [0,1] or alpha<=1-q+1e-13 or (alpha<1 and budget_gap<=q*q+1e-13) or (alpha>=1 and alpha==int(alpha) and budget_gap<=q+1e-13))
        assert predicted==zero_gap_region(1.,alpha,1.,q,B,U),(q,alpha,budget_gap)
        if predicted:
            max_equal_error=max(max_equal_error,abs(gap));assert abs(gap)<2e-10
        else:
            min_strict_gap=min(min_strict_gap,gap);assert gap>1e-10,(q,alpha,budget_gap,gap)
        count+=1
    mu=316912.;h=math.sqrt(3)*.1*mu;B=7925000.;U=1.1*B;p=25.;q=18.75
    max_currency_error=0.;figure_rows=0
    with (THEORY/'Figure_3_three_way_data.csv').open() as f:
      for r in csv.DictReader(f):
        alpha=float(r['delta_over_h']);D,R=global_frontier(h,p,q,h*alpha);C=min(B-R,U-D);opt=best_ramp(h,p,q,h*alpha,B,U);con=fixed_ramp(h,p,q,h*alpha,q);cons=min(B-con['reserve'],U-con['span'])
        expected=dict(global_H=C,ramp_H=opt['c'],conservative_H=cons,global_minus_ramp=C-opt['c'],global_minus_conservative=C-cons)
        for key,val in expected.items():max_currency_error=max(max_currency_error,abs(float(r[key])-val))
        figure_rows+=1
    assert figure_rows==901 and max_currency_error<=1e-6
    report=dict(equality_cases=count,max_predicted_equality_error_normalized_money=max_equal_error,minimum_tested_strict_gap_normalized_money=min_strict_gap,max_independent_implementation_discrepancy_normalized_money=max_formula_error,normalized_money_tolerance=2e-10,figure_rows=figure_rows,maximum_figure_error_eur=max_currency_error,figure_tolerance_eur=1e-6,authoritative_figure='theory/Figure_3_Three_Way_Frontier.png',authoritative_figure_code='theory/make_three_way_figure.py',all_passed=True)
    (OUT/'ramp_equality_and_figure_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return report
if __name__=='__main__':run()
