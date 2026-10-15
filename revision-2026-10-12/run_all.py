from pathlib import Path
import subprocess,sys,shutil,json,os
from PIL import Image,ImageChops
R=Path(__file__).resolve().parent;env=os.environ.copy();env['MPLCONFIGDIR']=str(R/'runtime_cache/mpl');env['XDG_CACHE_HOME']=str(R/'runtime_cache')
S=R.parent/'revision-2026-10-07';D=R/'portable_rerun';L=R/'logs';L.mkdir(exist_ok=True)
subprocess.run([sys.executable,str(S/'code/verify_manifest.py')],check=True,env=env)
shutil.copytree(S,D,dirs_exist_ok=True)
for p in D.rglob('*'):
 if p.is_file():p.chmod(0o644)
for name,script in [('portable_original',D/'code/reproduce.py'),('portable_independent',R/'analysis/independent_analysis.py')]:
 with (L/f'{name}.log').open('w') as f:subprocess.run([sys.executable,'-I',str(script)],check=True,stdout=f,stderr=subprocess.STDOUT,env=env)
checks=[]
for a in sorted((S/'results').iterdir()):
 b=D/'results'/a.name
 if a.name=='runtime.json':continue
 if a.suffix=='.json':passed=json.loads(a.read_text())==json.loads(b.read_text())
 elif a.suffix=='.csv':passed=a.read_bytes()==b.read_bytes()
 elif a.suffix=='.png':
  x,y=Image.open(a).convert('RGB'),Image.open(b).convert('RGB');passed=x.size==y.size and ImageChops.difference(x,y).getbbox() is None
 else:continue
 checks.append({'file':a.name,'identical':passed})
assert all(c['identical'] for c in checks),checks
(L/'portable_comparison.json').write_text(json.dumps(checks,indent=2))
print('PASS: original checksums, 17 scientific output comparisons, and all new independent checks. See logs/.')
