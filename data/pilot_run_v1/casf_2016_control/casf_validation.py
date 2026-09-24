"""Positive-control validation demanded by review round 5: run the
UNMODIFIED redocking-validation pipeline (same protein/ligand extraction,
same Open Babel PDBQT prep, same 3-seed exhaustiveness-16 Vina redocking,
same symmetry-corrected RMSD) on a real, externally-sourced,
non-cherry-picked subset of the CASF-2016 / PDBbind core set (285
complexes) -- the exact benchmark the reviewer's own cited 58% Vina
self-docking baseline (Zhang et al. 2023, JCIM, 10.1021/acs.jcim.3c00054)
is measured on. If this pipeline reproduces something close to that
baseline here, the 82%/93% CAK failure rate is a real finding about that
structure series, not a pipeline artifact. If it reproduces a similarly
low pass rate, the pipeline itself is the problem.

Unlike survey_structure.py, the ligand is NOT specified by the caller --
CASF-2016 entries were not looked up individually to avoid any appearance
of picking easy targets. Instead, the largest non-polymer, non-solvent,
non-common-crystallization-additive HETATM group in the structure is
selected automatically, using a fixed, pre-specified exclusion list (ions,
water, and standard buffer/cryoprotectant components) -- the same
principle PDBbind/CASF preprocessing pipelines use to identify "the"
ligand of interest in a crystal structure solved for a specific complex.

Usage: python3 casf_validation.py <pdbid>
"""
import gzip
import json
import subprocess
import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from spyrmsd import io, rmsd
from vina import Vina

warnings.filterwarnings("ignore")

S3_BASE = "https://pdbsnapshots.s3.us-west-2.amazonaws.com/20260101/pub/pdb"
EXHAUSTIVENESS = 16
SEEDS = [1000, 1001, 1002]
BOX_SIZE = 22.0

# Water, common ions, and common crystallization buffers/cryoprotectants --
# never the ligand of pharmacological interest. Fixed before any structure
# in this batch was inspected.
EXCLUDE_RESNAMES = {
    "HOH", "WAT", "DOD",
    "ZN", "MG", "CA", "NA", "K", "CL", "MN", "FE", "FE2", "CU", "CU1",
    "CO", "NI", "CD", "HG", "BR", "IOD", "F", "AL", "BA", "CS", "LI",
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE", "MPD", "DMS", "ACT",
    "TRS", "EPE", "HEPES", "BME", "DTT", "IMD", "FMT", "ACY", "CIT",
    "TAM", "BEZ", "1PE", "PEO", "P6G", "MES", "MRD", "BOG", "OCT", "NH4",
    "CO3", "NO3", "AZI", "SCN", "UNX", "UNL", "PEG2", "PEG3", "PEG4",
    "PE4", "P33", "PLM", "MYR", "OLA", "STE", "12P", "15P", "1PG",
    # Glycosylation sugars -- covalently attached N-/O-linked glycans, not
    # a freely-bound small-molecule ligand of pharmacological interest.
    # Missing this category caused a real mispick on 3dx1 (NAG selected
    # over the actual bound inhibitor YHO) -- found and fixed this round.
    "NAG", "NDG", "MAN", "BMA", "FUC", "FUL", "GAL", "GLA", "GLC", "BGC",
    "XYL", "XYS", "SIA", "NAN", "A2G", "FCA",
}


def fetch(pdbid, category, ext):
    d = pdbid[1:3]
    if category == "mmCIF":
        url = f"{S3_BASE}/data/structures/divided/mmCIF/{d}/{pdbid}.cif.gz"
    elif category == "validation":
        url = f"{S3_BASE}/validation_reports/{d}/{pdbid}/{pdbid}_validation.xml.gz"
    else:
        raise ValueError(category)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read()
        return gzip.decompress(data).decode("utf-8", errors="replace")
    except Exception as e:
        print(f"  [fetch failed] {url}: {e}")
        return None


def write_pdb_atom(f, serial, name, resname, chain, resseq, x, y, z, element):
    name_field = name if len(name) <= 4 else name[:4]
    f.write(f"HETATM{serial:>5} {name_field:<4} {resname:<3} {chain}{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n")


def process(pdbid, out_dir="."):
    out = Path(out_dir) / pdbid
    out.mkdir(parents=True, exist_ok=True)
    result = {"pdbid": pdbid}

    cif_text = fetch(pdbid, "mmCIF", "cif")
    if cif_text is None:
        result["error"] = "mmCIF fetch failed"
        return result
    cif_path = out / f"{pdbid}.cif"
    cif_path.write_text(cif_text)

    for line in cif_text.splitlines():
        if line.strip().startswith("_reflns.d_resolution_high") or \
           line.strip().startswith("_em_3d_reconstruction.resolution "):
            parts = line.split()
            try:
                result["resolution"] = float(parts[-1])
                break
            except ValueError:
                pass

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdbid, str(cif_path))
    model = structure[0]

    # Auto-detect the ligand: largest non-excluded HETATM residue.
    candidates = {}
    for chain in model:
        for res in chain:
            if res.id[0] == " ":
                continue
            resname = res.resname
            if resname in EXCLUDE_RESNAMES:
                continue
            n_heavy = sum(1 for a in res if (a.element or "").strip() != "H")
            if n_heavy < 5:
                continue
            candidates[f"{chain.id}{res.id[1]}_{resname}"] = (chain.id, res, n_heavy)
    if not candidates:
        result["error"] = "no candidate ligand found after exclusion filter"
        return result

    best_key = max(candidates, key=lambda k: candidates[k][2])
    chain_id, ligand_res, n_heavy = candidates[best_key]
    ligand_comp_id = ligand_res.resname
    result["ligand_comp_id"] = ligand_comp_id
    result["ligand_copy_used"] = best_key
    result["ligand_n_heavy_atoms"] = n_heavy
    result["n_candidates_considered"] = len(candidates)

    ligand_pdb = out / "ligand.pdb"
    atoms = [a for a in ligand_res if not (a.is_disordered() and a.get_altloc() not in ("A", " "))]
    with open(ligand_pdb, "w") as f:
        for i, atom in enumerate(atoms, start=1):
            elem = (atom.element or atom.get_name()[0]).strip()
            x, y, z = atom.coord
            write_pdb_atom(f, i, atom.get_name(), "LIG", "L", 1, x, y, z, elem)
        f.write("END\n")

    class ProteinSelect(Select):
        def accept_residue(self, residue):
            return residue.id[0] == " "

    io_ = PDBIO()
    io_.set_structure(structure)
    protein_pdb = out / "protein.pdb"
    io_.save(str(protein_pdb), ProteinSelect())

    ligand_sdf = out / "ligand.sdf"
    ligand_pdbqt = out / "ligand.pdbqt"
    protein_pdbqt = out / "protein.pdbqt"
    try:
        subprocess.run(["obabel", str(ligand_pdb), "-O", str(ligand_sdf)], check=True, capture_output=True, timeout=60)
        subprocess.run(["obabel", str(ligand_pdb), "-O", str(ligand_pdbqt), "--partialcharge", "gasteiger"],
                        check=True, capture_output=True, timeout=60)
        subprocess.run(["obabel", str(protein_pdb), "-O", str(protein_pdbqt), "-xr", "--partialcharge", "gasteiger"],
                        check=True, capture_output=True, timeout=120)
    except Exception as e:
        result["error"] = f"obabel prep failed: {e}"
        return result

    coords = np.array([a.coord for a in atoms])
    center = coords.mean(axis=0).tolist()
    result["box_center"] = center

    native_mol = io.loadmol(str(ligand_sdf))
    symm_vals, scores = [], []
    for seed in SEEDS:
        v = Vina(sf_name="vina", seed=seed)
        v.set_receptor(str(protein_pdbqt))
        try:
            v.set_ligand_from_file(str(ligand_pdbqt))
        except Exception as e:
            result["error"] = f"ligand load failed: {e}"
            return result
        v.compute_vina_maps(center=center, box_size=[BOX_SIZE] * 3)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=5)
        out_pdbqt = out / f"docked_seed{seed}.pdbqt"
        v.write_poses(str(out_pdbqt), n_poses=1, overwrite=True)
        energies = v.energies(n_poses=1)
        scores.append(float(energies[0][0]))

        docked_sdf = out / f"docked_seed{seed}.sdf"
        subprocess.run(["obabel", str(out_pdbqt), "-O", str(docked_sdf), "-f", "1", "-l", "1"],
                        check=True, capture_output=True)
        try:
            docked_mol = io.loadmol(str(docked_sdf))
            symm = float(rmsd.rmsdwrapper(native_mol, docked_mol, symmetry=True,
                                           minimize=False, strip=True)[0])
        except Exception as e:
            symm = float("nan")
        symm_vals.append(symm)

    result["symm_rmsd_values"] = symm_vals
    result["symm_rmsd_mean"] = float(np.nanmean(symm_vals))
    result["symm_rmsd_std"] = float(np.nanstd(symm_vals))
    result["vina_scores"] = scores
    result["passes_2.0A"] = bool(np.nanmean(symm_vals) <= 2.0)
    return result


if __name__ == "__main__":
    pdbid = sys.argv[1].lower()
    result = process(pdbid)
    print(json.dumps(result, indent=2))
    json.dump(result, open(f"{pdbid}_casf_result.json", "w"), indent=2)
