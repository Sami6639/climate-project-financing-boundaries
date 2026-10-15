"""One entry point; any exceeded tolerance raises and exits unsuccessfully."""
from pathlib import Path
import json,os,tempfile
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='interval-mpl-'))
os.environ.setdefault('XDG_CACHE_HOME',tempfile.mkdtemp(prefix='interval-cache-'))
import harden_legacy_tests,ramp_interval_checks,gaussian_interval_bounds,audit_ramp_equality_and_figure
P=Path(__file__).resolve().parent
if __name__=='__main__':
    report={}
    for name,module in [('preserved_legacy',harden_legacy_tests),('optimized_ramp',ramp_interval_checks),('equality_and_figure',audit_ramp_equality_and_figure),('gaussian_interval',gaussian_interval_bounds)]:
        print('Running',name,flush=True);report[name]=module.run()
    report['all_passed']=True
    (P/'strengthened_validation_summary.json').write_text(json.dumps(report,indent=2))
    print('All strengthened validation gates passed.',flush=True)
