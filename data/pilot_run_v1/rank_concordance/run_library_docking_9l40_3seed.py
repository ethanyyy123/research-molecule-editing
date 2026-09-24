"""Matched-design error-propagation re-run (review round 5, M4): the
original 9L40-vs-9L4B comparison used ONE seed for 9L40 but the existing
9L4B result averages THREE (2000, 2001, 2002) -- a one-seed ranking is
noisier by construction, so the weak rank correlation (Spearman rho=0.19)
was at least partly a design artifact, not purely a receptor-difference
effect. This re-runs 9L40 with the same three seeds as 9L4B for a matched
comparison.
"""
import json
import subprocess
import time

import pandas as pd
from vina import Vina

BOX_SIZE = 22.0
EXHAUSTIVENESS = 16
SEEDS = [2000, 2001, 2002]
RECEPTOR_PDBQT = "9l40_protein.pdbqt"
CENTER = json.load(open("box_centers.json"))["9l40"]["center"]

docking_set = pd.read_csv("docking_candidate_set.csv")

records = []
t_start = time.time()
for idx, row in docking_set.iterrows():
    smi_path = f"lig3_{row['compound_id']}.smi"
    pdbqt_path = f"lig3_{row['compound_id']}.pdbqt"
    open(smi_path, "w").write(row["smiles"] + "\n")
    try:
        subprocess.run(["obabel", smi_path, "-O", pdbqt_path, "--gen3d", "-p", "7.4"],
                        check=True, capture_output=True, timeout=60)
    except subprocess.CalledProcessError:
        print(f"  [skip] {row['compound_id']}: 3D embedding failed")
        continue

    for seed in SEEDS:
        v = Vina(sf_name="vina", seed=seed)
        v.set_receptor(RECEPTOR_PDBQT)
        try:
            v.set_ligand_from_file(pdbqt_path)
        except Exception as e:
            print(f"  [skip] {row['compound_id']} seed={seed}: {e}")
            continue
        v.compute_vina_maps(center=list(CENTER), box_size=[BOX_SIZE] * 3)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=5)
        energies = v.energies(n_poses=5)
        records.append({
            "compound_id": row["compound_id"], "source": row["source"],
            "smiles": row["smiles"], "seed": seed,
            "best_affinity_kcal_mol": float(energies[0][0]),
        })
    elapsed = time.time() - t_start
    print(f"[{idx+1}/{len(docking_set)}] {row['compound_id']} done "
          f"(elapsed {elapsed/60:.1f} min)", flush=True)

results = pd.DataFrame(records)
results.to_csv("atr_docking_results_9l40_3seed_raw.csv", index=False)

summary = (results.groupby(["compound_id", "source", "smiles"])["best_affinity_kcal_mol"]
           .agg(["mean", "std", "min"]).reset_index()
           .rename(columns={"mean": "affinity_mean", "std": "affinity_std", "min": "affinity_best"}))
summary = summary.sort_values("affinity_best")
summary.to_csv("atr_docking_results_9l40_3seed_summary.csv", index=False)

print(f"\nTotal time: {(time.time()-t_start)/60:.1f} min")
print(summary.to_string(index=False))
