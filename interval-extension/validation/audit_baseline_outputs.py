"""Fail-closed checks of the preserved original scientific benchmark suite.

Reports are checked individually, not merely by their stored summary. This
module raises RuntimeError on a count, non-finite value, or tolerance failure.
All geometric and LP tolerances refer to normalized money/slope units.
"""
from pathlib import Path
import json, math

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_ATOL = 1e-10
COST_ATOL = 1e-6
COST_RTOL = 1e-9

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def finite_tree(value, where='root'):
    if isinstance(value, float):
        require(math.isfinite(value), f'Non-finite number in {where}')
    elif isinstance(value, dict):
        for k, v in value.items():
            finite_tree(v, where + '.' + k)
    elif isinstance(value, list):
        for k, v in enumerate(value):
            finite_tree(v, f'{where}[{k}]')

def read(root, name):
    value = json.loads((Path(root) / name).read_text())
    finite_tree(value, name)
    return value

def verify(root=ROOT):
    stress = read(root, 'validation/stress_results.json')
    exact = stress['exact_construction']; lp = stress['lp_cases']
    require(len(exact) == 112, 'Expected all 112 construction cases')
    require(len(lp) == 44, 'Expected all 44 LP-bracket cases')
    for i, row in enumerate(exact):
        for key in ('mass_error','reserve_error','window_min_error','window_max_error'):
            require(abs(row[key]) <= GEOMETRY_ATOL, f'Construction case {i}: {key}={row[key]}')
    for i, row in enumerate(lp):
        for key in ('theory_below_lower','theory_above_upper','primal_violation','duality_gap'):
            require(0 <= row[key] <= GEOMETRY_ATOL, f'LP case {i}: {key}={row[key]}')
        require(row['lower_c'] <= row['theory_c'] + GEOMETRY_ATOL, f'LP case {i}: lower bound exceeds theory')
        require(row['theory_c'] <= row['upper_c'] + GEOMETRY_ATOL, f'LP case {i}: upper bound excludes theory')
    inadmissible = read(root, 'validation/inadmissible_cases.json')
    require(len(inadmissible) == 11, 'Expected all 11 infeasibility checks')
    for i,row in enumerate(inadmissible):
        require(row.get('actual_infeasible',row.get('actual_no_nonnegative_advance')) is True, f'Infeasibility case {i} failed')
    convergence = read(root, 'validation/convergence.json')
    require(len(convergence['records']) == 18, 'Expected all 18 convergence LP runs')
    scale = convergence['money_scale']
    for i,row in enumerate(convergence['records']):
        require(row['success'] is True, f'Convergence LP {i} unsuccessful')
        for key in ('max_primal_violation','primal_dual_gap'):
            require(0 <= row[key] <= GEOMETRY_ATOL, f'Convergence LP {i}: {key}={row[key]}')
        require(row['error_bound_eur']/scale >= -GEOMETRY_ATOL, f'Convergence LP {i}: inverted bracket')
        require(row['window_min'] >= 2*row['r']-GEOMETRY_ATOL, f'Convergence LP {i}: window constraint failed')
    cost = read(root, 'cost/cost_results.json')
    require(cost['random_positive_price_cases_passed'] == 1000, 'Expected all 1,000 seeded cost checks')
    require(len(cost['scenario_rows']) == 3, 'Expected all three baseline hypothetical price scenarios')
    for i,row in enumerate(cost['scenario_rows']):
        require(abs(row['weighted_resource_cost']-row['lp_cost']) <= COST_ATOL+COST_RTOL*abs(row['lp_cost']), f'Cost scenario {i}: independent LP mismatch')
        require(row['budget_upgrade']+row['a']*row['external_nonrepayable_funds'] >= row['budget_deficit']-COST_ATOL, f'Cost scenario {i}: budget deficit')
        require(row['exposure_upgrade']+row['a']*row['external_nonrepayable_funds'] >= row['exposure_deficit']-COST_ATOL, f'Cost scenario {i}: exposure deficit')
    return {'status':'PASS','construction_cases':112,'lp_bracket_cases':44,'infeasibility_cases':11,'convergence_LP_cases':18,'seeded_cost_cases':1000,'geometry_absolute_tolerance':GEOMETRY_ATOL,'cost_absolute_tolerance':COST_ATOL,'cost_relative_tolerance':COST_RTOL,'seed':20261014,'max_exact_error':max(abs(r[k]) for r in exact for k in ('mass_error','reserve_error','window_min_error','window_max_error')),'max_bracket_violation':max(max(r['theory_below_lower'],r['theory_above_upper']) for r in lp),'max_primal_violation':max(r['primal_violation'] for r in lp),'max_duality_gap':max(r['duality_gap'] for r in lp)}

if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
