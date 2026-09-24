"""Run ONE seed of the flexible-sidechain redocking test and print a JSON
result line.

The native crash (terminate called after throwing an instance of
'swig::stop_iteration') traced back NOT to reusing multiple Vina() objects
in one process, but to import order: `import vina` followed by
`from spyrmsd import ...` aborts at the `spyrmsd` import itself (confirmed
by a step-by-step trace -- no user code had run yet), while importing
spyrmsd first and vina second works every time. Both SWIG-wrapped native
extensions apparently register conflicting C++ exception-translation state
at load time. Fixed by reordering the imports below; per-seed subprocess
isolation is kept anyway as a defensive measure.
"""
import json
import subprocess
import sys

import numpy as np
from spyrmsd import io, rmsd
from vina import Vina

BOX_CENTER = [140.15, 208.575, 172.671]
BOX_SIZE = 30.0
EXHAUSTIVENESS = 32
RIGID_PDBQT = "9l40_flexdock_rigid.pdbqt"
FLEX_PDBQT = "9l40_flexdock_flex.pdbqt"
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


def main(seed):
    v = Vina(sf_name="vina", seed=seed)
    v.set_receptor(RIGID_PDBQT, FLEX_PDBQT)
    v.set_ligand_from_file(LIGAND_PDBQT)
    v.compute_vina_maps(center=BOX_CENTER, box_size=[BOX_SIZE] * 3)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)

    out_pdbqt = f"9l40_flexdock_seed{seed}_docked.pdbqt"
    v.write_poses(out_pdbqt, n_poses=1, overwrite=True)
    energies = v.energies(n_poses=1)
    score = float(energies[0][0])

    native_coords, native_names = parse_pdbqt_coords(LIGAND_PDBQT, pose_index=0)
    all_coords, all_names = parse_pdbqt_coords(out_pdbqt, pose_index=0)
    docked_coords = all_coords[: len(native_coords)]
    docked_names = all_names[: len(native_coords)]

    naive = float("nan")
    if docked_names == native_names:
        diff = native_coords - docked_coords
        naive = float(np.sqrt((diff ** 2).sum(axis=1).mean()))

    ligand_pdb = f"9l40_flexdock_seed{seed}_ligand.pdb"
    with open(ligand_pdb, "w") as fh:
        for i, ((x, y, z), name) in enumerate(zip(docked_coords, docked_names), start=1):
            elem = "".join(c for c in name if c.isalpha())[:1] or "C"
            fh.write(f"HETATM{i:>5} {name:<4} LIG L   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2}\n")
        fh.write("END\n")
    docked_sdf = f"9l40_flexdock_seed{seed}_ligand.sdf"
    subprocess.run(["obabel", ligand_pdb, "-O", docked_sdf], check=True, capture_output=True)

    native_mol = io.loadmol(NATIVE_SDF)
    docked_mol = io.loadmol(docked_sdf)
    symm = float(rmsd.rmsdwrapper(native_mol, docked_mol, symmetry=True, minimize=False, strip=True)[0])

    print(json.dumps({"seed": seed, "vina_score": score, "naive_rmsd": naive, "symm_rmsd": symm}))


if __name__ == "__main__":
    main(int(sys.argv[1]))
