"""Batch driver for casf_validation.py over the first 20 PDB IDs (in the
order they appear in the externally-sourced CASF-2016 core-set CSV --
https://github.com/OptiMaL-PSE-Lab/DeepDock, DockingResults_CASF2016_CoreSet.csv)
-- an objective, non-cherry-picked selection rule fixed before any of these
structures was inspected individually.
"""
import json
import traceback

from casf_validation import process

PDBIDS = [
    "4k18", "4qac", "1o3f", "4ih7", "3dx1", "1syi", "2p4y", "3nq9", "3wtj",
    "4w9i", "3b65", "3pww", "4w9c", "3kwa", "3arq", "3d6q", "2w4x", "2j7h",
    "2xbv", "2cet",
]

for pdbid in PDBIDS:
    print(f"\n=== {pdbid} ===", flush=True)
    try:
        result = process(pdbid, out_dir=".")
    except Exception as e:
        result = {"pdbid": pdbid, "error": f"{type(e).__name__}: {e}"}
        traceback.print_exc()
    json.dump(result, open(f"{pdbid}_casf_result.json", "w"), indent=2)
    print(json.dumps(result, indent=2), flush=True)

print("\n=== CASF BATCH DONE ===", flush=True)
