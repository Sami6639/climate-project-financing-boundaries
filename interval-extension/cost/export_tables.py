import json,csv
from pathlib import Path
P=Path(__file__).parent
rows=json.loads((P/'cost_results.json').read_text())['scenario_rows']
fields=['scenario','cb','cu','cs','budget_upgrade','exposure_upgrade','external_nonrepayable_funds','donor_advance','total_initial_project_cash','weighted_resource_cost']
with (P/'cost_table.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r[k] for k in fields} for r in rows)
with (P/'all_candidate_comparisons.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['scenario','candidate_external_funds','budget_upgrade','exposure_upgrade','weighted_resource_cost'])
 for r in rows:
  for s in r['candidate_s']:
   b=max(0,r['budget_deficit']-r['a']*s);u=max(0,r['exposure_deficit']-r['a']*s)
   w.writerow([r['scenario'],s,b,u,r['cb']*b+r['cu']*u+r['cs']*s])
