"""Run survey_structure.process() across all remaining CAK candidate
structures, one at a time, writing each result to <pdbid>_result.json as it
completes (so partial progress is never lost even if a later structure
fails or the batch is interrupted).
"""
import json
import traceback

from survey_structure import process

CANDIDATES = [
    ("8p6v", "I74"),
    ("8p6x", "NS9"),
    ("8p6z", "X4L"),
    ("8p70", "X4Q"),
    ("8p71", "X4F"),
    ("8p72", "X2Z"),
    ("8p74", "X3Z"),
    ("8p75", "X3Z"),
    ("8p76", "X1W"),
    ("8p77", "I73"),
    ("8p78", "1QK"),
    ("8p7l", "WD9"),
]

for pdbid, comp_id in CANDIDATES:
    print(f"\n=== {pdbid} ({comp_id}) ===", flush=True)
    try:
        result = process(pdbid, comp_id, out_dir=".")
    except Exception as e:
        result = {"pdbid": pdbid, "ligand_comp_id": comp_id, "error": f"{type(e).__name__}: {e}"}
        traceback.print_exc()
    json.dump(result, open(f"{pdbid}_result.json", "w"), indent=2)
    print(json.dumps(result, indent=2), flush=True)

print("\n=== BATCH DONE ===", flush=True)
