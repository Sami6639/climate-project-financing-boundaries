from pathlib import Path
import hashlib,json,sys
R=Path(__file__).resolve().parent
bad=[]
for name,digest in json.loads((R/'manifest_sha256.json').read_text()).items():
 p=R/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:bad.append(name)
if bad: print('FAIL:',bad);sys.exit(1)
print('PASS: publication snapshot file hashes match')
