"""
Check that every gene of processor/builtin_markers.json is a PanglaoDB human marker of one of
the set's PanglaoDB cell types. Maintenance only (needs network; not part of the image or tests).

PanglaoDB's own download refuses scripted clients, so the marker table is read from OmniPath,
which redistributes it unchanged:

  curl -sSL -o panglao.tsv "https://omnipathdb.org/annotations?resources=PanglaoDB&format=tsv"
  python3 scripts/verify_builtin_markers.py panglao.tsv
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

BUILTIN = Path(__file__).resolve().parent.parent / "processor" / "builtin_markers.json"


def main(path: str) -> int:
    records: dict[str, dict[str, str]] = defaultdict(dict)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            records[(r["record_id"], r["genesymbol"])][r["label"]] = r["value"]
    human: dict[str, set[str]] = defaultdict(set)    # gene -> PanglaoDB cell types (human markers)
    for (_, gene), rec in records.items():
        if rec.get("human") == "True" and rec.get("cell_type"):
            human[gene].add(rec["cell_type"])
    sets = json.loads(BUILTIN.read_text())["sets"]
    bad = 0
    for s in sets:
        for g in s["genes"]:
            if not human[g] & set(s["panglaodb"]):
                bad += 1
                print(f"{s['cell_type']}: {g} is not a PanglaoDB human marker of {s['panglaodb']} "
                      f"(listed for {sorted(human[g]) or 'nothing'})")
    print(f"{sum(len(s['genes']) for s in sets)} genes in {len(sets)} sets checked, {bad} not found")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
