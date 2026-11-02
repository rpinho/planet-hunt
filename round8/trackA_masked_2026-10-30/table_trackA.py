"""Print the Track A table using only saved JSONs and the standard library."""
import json
from pathlib import Path

root = Path(__file__).resolve().parent
rows = []
for path in sorted((root / "results").glob("A_*.json")):
    if path.name == "A_status.json":
        continue
    data = json.loads(path.read_text())
    z = data["cases"]["A"]["z_null"]
    edge = data["cases"]["unmasked"]["edge_share"]["A"]
    assert data["complete"] and data["bar_stays"] == (z >= 2.5 and edge <= 0.5)
    rows.append((data["signal"]["tic"], data["signal"]["P"], z, data["bar_stays"]))
status = json.loads((root / "results" / "A_status.json").read_text())
assert status["complete"] and len(rows) == status["n_target"] == status["n_measured"] == 17
assert sum(row[3] for row in rows) == status["k"] == 14
print("| TIC | Period (days) | z_A | Kept |")
print("|---:|---:|---:|:---:|")
for tic, period, z, kept in rows:
    print(f"| {tic} | {period} | {z} | {'yes' if kept else 'no'} |")
print(f"\n{sum(row[3] for row in rows)}/{len(rows)} stay; {len(rows)}/{len(rows)} measured.")
