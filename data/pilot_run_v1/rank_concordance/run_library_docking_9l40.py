"""Error-propagation demonstration (review round 4, JCAMD strategy item 3):
dock the SAME 40-compound pilot set (identical SMILES, identical box size,
identical exhaustiveness) against the structure the redocking validation
control excluded (9L40, failed at 2.53 A symmetry-corrected RMSD) instead
of the validated one (9L4B), to show what the scoring landscape would have
looked like had the validation control not been run and 9L40 been trusted
blindly for the pilot screen.

One seed per compound (not three) to keep this within budget -- this is a
demonstration of scoring-landscape distortion via rank correlation against
the already-3-seed-averaged 9L4B result, not a new primary campaign
requiring its own replicate-averaging.
"""
import json
import subprocess
import time

import pandas as pd
from vina import Vina

BOX_SIZE = 22.0
EXHAUSTIVENESS = 16
SEED = 2000
RECEPTOR_PDBQT = "9l40_protein.pdbqt"
CENTER = json.load(open("box_centers.json"))["9l40"]["center"]

docking_set = pd.read_csv("docking_candidate_set.csv")

records = []
t_start = time.time()
for idx, row in docking_set.iterrows():
    smi_path = f"lig_{row['compound_id']}.smi"
    pdbqt_path = f"lig_{row['compound_id']}.pdbqt"
    open(smi_path, "w").write(row["smiles"] + "\n")
    try:
        subprocess.run(["obabel", smi_path, "-O", pdbqt_path, "--gen3d", "-p", "7.4"],
                        check=True, capture_output=True, timeout=60)
    except subprocess.CalledProcessError:
        print(f"  [skip] {row['compound_id']}: 3D embedding failed")
        continue

    v = Vina(sf_name="vina", seed=SEED)
    v.set_receptor(RECEPTOR_PDBQT)
    try:
        v.set_ligand_from_file(pdbqt_path)
    except Exception as e:
        print(f"  [skip] {row['compound_id']}: {e}")
        continue
    v.compute_vina_maps(center=list(CENTER), box_size=[BOX_SIZE] * 3)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=5)
    energies = v.energies(n_poses=5)
    records.append({
        "compound_id": row["compound_id"], "source": row["source"],
        "smiles": row["smiles"], "seed": SEED,
        "best_affinity_kcal_mol_9l40": float(energies[0][0]),
    })
    elapsed = time.time() - t_start
    print(f"[{idx+1}/{len(docking_set)}] {row['compound_id']} done "
          f"(elapsed {elapsed/60:.1f} min)", flush=True)

results = pd.DataFrame(records)
results.to_csv("atr_docking_results_9l40_raw.csv", index=False)
print(f"\nTotal time: {(time.time()-t_start)/60:.1f} min")
print(results.sort_values('best_affinity_kcal_mol_9l40').to_string(index=False))
