"""Matched control for the flexible-sidechain redocking test (review round
5, M2): the flexible run (flexible_redocking_results.json) differs from the
original rigid redocking result (2.53 A, 22 A box, Open-Babel-prepared
ligand) in THREE ways at once -- side-chain flexibility, box size (22->30
A), and receptor composition (Meeko's automatic removal of ~400
unmatched-template residues). This script isolates flexibility alone by
running fully RIGID redocking with the same 30 A box and the same
Meeko-prepared receptor/ligand used in the flexible run (just without
declaring any residue flexible), so it can be compared directly against
both the original rigid result (2.53 A) and the flexible result (3.00 A).
"""
import json
import subprocess

import numpy as np
from spyrmsd import io, rmsd
from vina import Vina

BOX_CENTER = [140.15, 208.575, 172.671]
BOX_SIZE = 30.0
EXHAUSTIVENESS = 32
SEEDS = [1000, 1001, 1002]

# The Meeko-prepared receptor's RIGID part alone (no flex file passed) --
# same atom composition (post template-match filtering) as the flexible run.
RIGID_PDBQT = "9l40_flexdock_rigid.pdbqt"
LIGAND_PDBQT = "9l40_ligand_meeko.pdbqt"
NATIVE_SDF = "9l40_active_site_ligand_fixed.sdf"


def parse_pdbqt_coords(path, pose_index=0):
    coords, names = [], []
    current_pose = 0
    with open(path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                current_pose += 1
                continue
            if current_pose != pose_index:
                continue
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x, y, z))
                names.append(line[12:16].strip())
    return np.array(coords), names


native_mol = io.loadmol(NATIVE_SDF)
results = {"naive_rmsd": [], "symm_rmsd": [], "vina_score": []}

for seed in SEEDS:
    v = Vina(sf_name="vina", seed=seed)
    v.set_receptor(RIGID_PDBQT)  # NOTE: rigid only -- no flex file
    v.set_ligand_from_file(LIGAND_PDBQT)
    v.compute_vina_maps(center=BOX_CENTER, box_size=[BOX_SIZE] * 3)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)

    out_pdbqt = f"9l40_matched_rigid_seed{seed}_docked.pdbqt"
    v.write_poses(out_pdbqt, n_poses=1, overwrite=True)
    energies = v.energies(n_poses=1)

    native_coords, native_names = parse_pdbqt_coords(LIGAND_PDBQT, pose_index=0)
    docked_coords, docked_names = parse_pdbqt_coords(out_pdbqt, pose_index=0)

    naive = float("nan")
    if docked_names == native_names:
        diff = native_coords - docked_coords
        naive = float(np.sqrt((diff ** 2).sum(axis=1).mean()))

    ligand_pdb = f"9l40_matched_rigid_seed{seed}_ligand.pdb"
    with open(ligand_pdb, "w") as fh:
        for i, ((x, y, z), name) in enumerate(zip(docked_coords, docked_names), start=1):
            elem = "".join(c for c in name if c.isalpha())[:1] or "C"
            fh.write(f"HETATM{i:>5} {name:<4} LIG L   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2}\n")
        fh.write("END\n")
    docked_sdf = f"9l40_matched_rigid_seed{seed}_ligand.sdf"
    subprocess.run(["obabel", ligand_pdb, "-O", docked_sdf], check=True, capture_output=True)
    docked_mol = io.loadmol(docked_sdf)
    symm = float(rmsd.rmsdwrapper(native_mol, docked_mol, symmetry=True, minimize=False, strip=True)[0])

    results["naive_rmsd"].append(naive)
    results["symm_rmsd"].append(symm)
    results["vina_score"].append(float(energies[0][0]))
    print(f"seed={seed}: score={energies[0][0]:.3f}, naive={naive:.3f}, symm={symm:.3f}", flush=True)

summary = {
    "description": "Matched rigid control: 30 A box + Meeko-prepared receptor/ligand, NO flexible residues",
    "naive_rmsd_mean": float(np.mean(results["naive_rmsd"])),
    "symm_rmsd_mean": float(np.mean(results["symm_rmsd"])),
    "symm_rmsd_std": float(np.std(results["symm_rmsd"])),
    "comparison_original_rigid_22A_openbabel": 2.53,
    "comparison_flexible_30A_meeko": 3.00,
    "raw": results,
}
json.dump(summary, open("matched_rigid_control_results.json", "w"), indent=2)
print(json.dumps(summary, indent=2))
