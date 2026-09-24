"""Process one PDB ID through the redocking-validation pipeline for the
n>1 survey (review round 4): download mmCIF from the AWS PDB mirror,
auto-detect the bound-ligand chemical component, extract protein+ligand,
redock (3 seeded replicates), pull wwPDB Q-score/residue_inclusion, and
record resolution/deposition year. Reuses the exact protocol validated on
9L40/9L4B/2YM8 (contact-analysis disambiguation, placeholder residue name
for 5-char comp_ids, symmetry-corrected RMSD) so results are comparable.

Usage: python3 survey_structure.py <pdbid> <expected_ligand_comp_id>
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
EXHAUSTIVENESS = 16  # matches the confirmed-consistent setting from round 3
SEEDS = [1000, 1001, 1002]
BOX_SIZE = 22.0


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


def process(pdbid, ligand_comp_id, out_dir="."):
    out = Path(out_dir) / pdbid
    out.mkdir(parents=True, exist_ok=True)
    result = {"pdbid": pdbid, "ligand_comp_id": ligand_comp_id}

    cif_text = fetch(pdbid, "mmCIF", "cif")
    if cif_text is None:
        result["error"] = "mmCIF fetch failed"
        return result
    cif_path = out / f"{pdbid}.cif"
    cif_path.write_text(cif_text)

    for line in cif_text.splitlines():
        if line.strip().startswith("_em_3d_reconstruction.resolution ") or \
           line.strip().startswith("_reflns.d_resolution_high"):
            parts = line.split()
            try:
                result["resolution"] = float(parts[-1])
                break
            except ValueError:
                pass
    for line in cif_text.splitlines():
        if "_pdbx_database_status.recvd_initial_deposition_date" in line:
            result["deposition_date"] = line.split()[-1].strip("'\"")

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdbid, str(cif_path))
    model = structure[0]

    ligand_copies = {}
    for chain in model:
        for res in chain:
            if res.resname == ligand_comp_id:
                ligand_copies[f"{chain.id}{res.id[1]}"] = (chain.id, res)
    if not ligand_copies:
        result["error"] = f"ligand {ligand_comp_id} not found"
        return result

    # If multiple copies, pick the one with fewest cross-chain contacts
    # (most likely the canonical single-site pocket, not an interface/
    # secondary site) -- same principle as the 9L40 A2701 vs A2702 check.
    protein_atoms_by_chain = {}
    for chain in model:
        protein_atoms_by_chain[chain.id] = [a for r in chain if r.id[0] == " " for a in r]

    def cross_chain_contacts(chain_id, res):
        own = sum(1 for a in res if any((a - pa) < 4.5 for pa in protein_atoms_by_chain.get(chain_id, [])))
        other = 0
        for cid, atoms in protein_atoms_by_chain.items():
            if cid == chain_id:
                continue
            other += sum(1 for a in res if any((a - pa) < 4.5 for pa in atoms))
        return other

    best_key = min(ligand_copies, key=lambda k: cross_chain_contacts(*ligand_copies[k]))
    chain_id, ligand_res = ligand_copies[best_key]
    result["ligand_copy_used"] = best_key
    result["n_ligand_copies_total"] = len(ligand_copies)

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
    subprocess.run(["obabel", str(ligand_pdb), "-O", str(ligand_sdf)], check=True, capture_output=True)
    subprocess.run(["obabel", str(ligand_pdb), "-O", str(ligand_pdbqt), "--partialcharge", "gasteiger"],
                    check=True, capture_output=True)
    subprocess.run(["obabel", str(protein_pdb), "-O", str(protein_pdbqt), "-xr", "--partialcharge", "gasteiger"],
                    check=True, capture_output=True)

    coords = np.array([a.coord for a in atoms])
    center = coords.mean(axis=0).tolist()
    result["box_center"] = center

    native_mol = io.loadmol(str(ligand_sdf))
    symm_vals, naive_vals, scores = [], [], []
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

    # wwPDB validation metrics for the ligand copy used. Attribute order in
    # the real XML is NOT resname-then-resnum (it's actually
    # ...model/chain/resnum/said/ent/seq/icode/resname/altcode/...), so a
    # regex that requires resname to appear before resnum never matches --
    # found by inspecting the raw 8P73 validation XML directly. Fixed by
    # extracting each <ModelledSubgroup .../> tag whole and matching
    # attributes as an order-independent dict, filtered on chain+resnum+
    # resname together (resnum alone is not unique across chains).
    val_xml = fetch(pdbid, "validation", "xml")
    if val_xml:
        import re
        resnum = int("".join(c for c in best_key if c.isdigit()))
        for tag in re.findall(r'<ModelledSubgroup\b[^>]*/?>', val_xml):
            attrs = dict(re.findall(r'(\w+)="([^"]*)"', tag))
            if (attrs.get("resname") == ligand_comp_id
                    and attrs.get("resnum") == str(resnum)
                    and attrs.get("chain") == chain_id):
                result["Q_score"] = float(attrs.get("Q_score", "nan"))
                result["residue_inclusion"] = float(attrs.get("residue_inclusion", "nan"))
                break

    return result


if __name__ == "__main__":
    pdbid, ligand_comp_id = sys.argv[1], sys.argv[2]
    result = process(pdbid.lower(), ligand_comp_id)
    print(json.dumps(result, indent=2))
    json.dump(result, open(f"{pdbid.lower()}_result.json", "w"), indent=2)
