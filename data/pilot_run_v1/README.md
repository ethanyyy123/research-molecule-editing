# Pilot Run v1 — Executed Results

All files in this directory are actual outputs from the cognate redocking
validation pipeline. See the repository README for the full data mapping
to manuscript tables and sections.

## Subdirectories

- `ligands/` — Co-crystallized ligands extracted from PDB structures
- `validation_analysis/` — Symmetry-corrected RMSD, wwPDB reports, redocked poses
- `cak_survey/` — 14-structure CDK7/cyclin H/MAT1 survey (manuscript Table 2)
- `flexible_redocking/` — Flexible-sidechain redocking test (manuscript §3.4)
- `rank_concordance/` — 9L40 vs. 9L4B rank comparison (manuscript §3.5)
- `casf_2016_control/` — 20-complex CASF-2016 positive control (manuscript §3.3)

Receptor structures are not included (regenerable from RCSB PDB AWS mirror,
snapshot 2026-01-01). See repository README for PDB IDs and retrieval details.
