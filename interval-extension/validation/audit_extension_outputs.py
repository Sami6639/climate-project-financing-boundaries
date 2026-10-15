"""Independent report-level gates for new extension outputs.

Scientific modules perform their own detailed numerical checks. These gates
also reject stale/missing cases, failed report rows and non-finite output.
"""
import csv
from pathlib import Path
from audit_baseline_outputs import read, require
ROOT=Path(__file__).resolve().parents[1]

def csv_rows(root,name):
    with (Path(root)/name).open(newline='') as f:
        return list(csv.DictReader(f))

def verify(root=ROOT):
    ramp=read(root,'theory/single_ramp_self_check.json')
    require(ramp['status']=='PASS' and not ramp['failures'],'Single-ramp formula audit failed')
    require(ramp['row_count']==517,'Expected all 517 deterministic ramp comparison cases')
    require(ramp['maximum_scaled_order_error']<=2e-10,'Single-ramp formula ordering error exceeds tolerance')
    independent=read(root,'validation/optimized_ramp_checks.json')
    s=independent['summary']
    require(s['independent_random_cases']==240,'Expected 240 independent random ramp cases')
    require(s['comparison_rows']==189,'Expected 189 independent ramp comparison rows')
    require(s['targeted_boundary_cases']==94 and len(independent['boundary_cases'])==94,'Expected 94 ramp boundary cases')
    require(s['max_dense_search_excess_over_formula']<=2e-9,'Independent ramp search mismatch')
    require(s['max_independent_integral_error']<=2e-9,'Independent ramp integration mismatch')
    require(s['negative_capacity_not_clipped']<0,'Negative financing capacity was silently clipped')
    for row in independent['checks']:
        require(row['passed'] is True and abs(row['error'])<=row['threshold'],'Independent ramp report contains a failed case')
    for row in independent['boundary_cases']:
        require(row['passed'] is True,'Boundary report contains a failed case')
        require(abs(row['mean_error_currency'])<=row['mean_tolerance_currency'],'Boundary expected payment error')
        require(row['marginal_shortfall_currency_per_emissions']<=row['marginal_tolerance_currency_per_emissions'],'Boundary marginal constraint failure')
    gaussian=read(root,'validation/gaussian_interval_certificates.json')
    require(gaussian['summary']['cases']==15 and len(gaussian['records'])==15,'Expected 15 Gaussian interval certificates')
    for row in gaussian['records']:
        lo=row['gaussian_capacity_lower_over_ph'];hi=row['gaussian_capacity_upper_over_ph']
        require(0<=hi-lo<.075,'Gaussian numerical enclosure invalid or wider than target')
        require(row['gaussian_continuum_marginal_lower_over_p']>=row['q_over_p'],'Gaussian continuum marginal certificate failed')
        f=row['feasible_certificate']
        require(f['C']>=0 and f['C']+f['span']<=row['U_over_ph']+2e-9 and f['C']+f['reserve']<=row['B_over_ph']+2e-9,'Gaussian feasible resource certificate failed')
    cost=read(root,'cost/price_region_validation_report.json')
    require(cost['status']=='PASS','Price-region report failed')
    for key,expected in [('legacy_original_random_checks_retained',1000),('original_same_cases_new_price_rule_checks',1000),('purposeful_edge_cases',174),('sensitivity_rows',73),('derivative_checks',72),('main_rows',3)]:
        require(cost[key]==expected,f'Price-region count mismatch: {key}')
    for key,tol in [('max_normalized_objective_error',1e-8),('max_normalized_primal_violation',1e-8),('max_normalized_optimal_face_error',1e-7),('max_derivative_relative_error',2e-5)]:
        require(0<=cost[key]<=tol,f'Price-region tolerance violation: {key}')
    for name,expected in [('original_1000_price_region_checks.csv',1000),('purposeful_price_region_checks.csv',174),('cost_sensitivity.csv',73),('cost_derivative_checks.csv',72),('main_price_examples.csv',3)]:
        rows=csv_rows(root,'cost/'+name)
        require(len(rows)==expected,'Price-region CSV case count mismatch: '+name)
        for row in rows:
            require(row.get('passed')=='True','Failed price-region CSV row: '+name)
    return {'status':'PASS','deterministic_ramp_cases':517,'independent_ramp_random_cases':240,'independent_ramp_comparison_rows':189,'ramp_boundary_cases':94,'gaussian_numerical_certificates':15,'gaussian_max_enclosure_width':gaussian['summary']['max_bound_gap'],'cost_original_cases_retested':1000,'cost_targeted_edge_cases':174,'cost_sensitivity_rows':73,'cost_derivative_checks':72}

if __name__=='__main__':
    import json
    print(json.dumps(verify(),indent=2))
