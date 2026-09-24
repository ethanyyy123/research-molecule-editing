# Review round 5 — external re-review of the round-4 survey

This directory holds the real data behind the response to a fifth
external review, which caught two real errors (an arithmetic mistake in
the CAK fail count, and a mischaracterization of the 8P74/8P75 structure
pair) and one confounded experiment (the flexible-sidechain redocking
comparison), and demanded a decisive positive control before the
survey's 82%/93% redocking-failure claim could be trusted: run the same,
unmodified pipeline on an external, unbiased benchmark set and compare
its pass rate to the literature.

- `casf_validation.py` — the redocking-validation pipeline generalized
  to auto-detect the ligand of interest (largest non-excluded heteroatom
  group, after excluding water/ions/common crystallization
  additives/glycosylation sugars) rather than requiring it as an
  argument, since CASF-2016 entries were not looked up individually
  (avoiding any appearance of picking easy targets). `run_casf_batch.py`
  runs it over the first 20 PDB IDs (in file order) of a real,
  externally hosted CASF-2016 core-set listing
  (`https://github.com/OptiMaL-PSE-Lab/DeepDock`,
  `Validation_Docking/DockingResults_CASF2016_CoreSet.csv`) -- an
  objective selection rule fixed before any of these 20 structures was
  inspected individually. `*_casf_result.json` are the per-structure
  results; `casf_combined_results.json` merges them.

  **Result: 8/20 pass (40%).** Statistically indistinguishable from
  AutoDock Vina's documented 58% self-docking success rate on the full
  285-complex PDBbind Core Set (one-sided binomial p=0.08) -- so this
  pipeline is not badly broken. Against this pipeline's *own* measured
  40% baseline, the CAK survey's 1/14 (7%) passing rate remains a
  significant outlier (p=0.008), while the original ATR pair's 1/2 (50%)
  is fully consistent with a 40% baseline (p=0.84). The CAK finding
  survives the review's own proposed decisive test.

  One real bug was found and fixed while reviewing this batch: 3DX1's
  auto-detected ligand was initially a glycosylation sugar (NAG) rather
  than the deposited inhibitor (YHO), because glycan chemical component
  codes were not on the original exclusion list. Fixed by adding a
  glycan-code category; the corrected result (3.29 Å against the real
  ligand) still fails, so this bug did not change the pass/fail count,
  but would have been a real error in the paper's Table 6 had it gone
  unnoticed.

- `run_rigid_control_matched.py` — the M2 matched control for the
  flexible-sidechain redocking test: fully rigid redocking with the same
  enlarged 30 Å box and the same Meeko-prepared receptor/ligand used in
  the flexible run, but no residues declared flexible. Result: mean
  symmetry-corrected RMSD 4.50±0.05 Å, *worse* than both the original
  rigid result (2.53 Å, 22 Å box, Open Babel prep) and the flexible
  result (3.00 Å) — meaning the original comparison (flexible vs.
  original-rigid) was confounded by box size and receptor preparation,
  not just flexibility. Against this correctly matched baseline,
  flexibility genuinely improves the pose (4.50→3.00 Å), the opposite
  conclusion from what the paper originally reported.

See `../error_propagation_9l40/run_library_docking_9l40_3seed.py` for the
M4 matched (3-seed) re-run of the error-propagation comparison, addressing
the review's concern that the original 1-seed-vs-3-seed comparison was
asymmetric.

See `docs/08_review_response_round5.md` for the full response, including
the two claims withdrawn (the CAK fail count arithmetic, and the
8P74/8P75 characterization) and the reasoning behind what was kept.
