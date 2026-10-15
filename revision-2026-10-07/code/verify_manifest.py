"""Verify deposited package bytes before regenerating results."""
from pathlib import Path
import csv
import hashlib

root = Path(__file__).resolve().parents[1]
failures = []
with (root / "manifest_sha256.csv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
for row in rows:
    path = root / row["path"]
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        failures.append(row["path"])
if failures:
    raise SystemExit("Integrity check failed: " + ", ".join(failures))
print(f"Verified {len(rows)} files against manifest_sha256.csv")
