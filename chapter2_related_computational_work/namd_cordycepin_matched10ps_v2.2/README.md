# NAMD Cordycepin MAPK/ADK — matched 10-ps analysis pipeline v2.2

This directory contains the **final post-production analysis and provenance code** used for the Chapter 2 MAPK/ADK matched-temporal-resolution workflow.

## Sampling policy
- MAPK raw trajectory: 10 ps/frame; analysis stride = 1; 10,000 production frames.
- ADK raw trajectory: 2 ps/frame; analysis stride = 5; 10,000 analysis frames.
- Common analysis grid: 0.01–100.00 ns at 10 ps spacing.
- Production only; NVT/NPT are excluded from reported structural observables.
- Raw trajectory files are read-only. No interpolation, smoothing, filtering, imputation, outlier removal, or frame balancing is applied.

## Code included
The complete final analysis chain in this directory is indexed in `CODE_MANIFEST.tsv`:
1. exact T1/T2/T3 state discovery;
2. topology lineage and exact snapshot generation;
3. receptor-only/apo topology-mode audit;
4. matched-10ps RMSD, backbone RMSD, RMSF and Rg;
5. exact PCA, RMSIP and non-forced cluster-quality diagnostics;
6. publication/thesis figure generation;
7. updated Methods/Results/Discussion generation;
8. environment and end-to-end launcher;
9. known-answer verification scripts for the provenance/sampling gates.

## Interpretation lock
The verified 100-ns production systems are receptor-only/apo trajectories. They support receptor conformational sampling and do **not** constitute ligand-bound MD. T1/T2/T3 docking states are exact 0/50/100 ns receptor boundary states and are docked independently.

## Clustering lock
For the final real data, Calinski–Harabasz and Davies–Bouldin criteria do not select the same k; therefore no unique cluster count or representative structure is forced.

## Reproducibility
The final real-data matched-10ps Stage 03 and Stage 04 runs passed their internal 10-PASS checks. Generated analysis outputs and multi-GB raw DCD trajectories are intentionally not committed in this code directory. Output ZIPs retain SHA-256 manifests.

See `CODE_MANIFEST.tsv` for the code inventory and the stage-specific `references/REFERENCES_VERIFIED.tsv` files for frozen bibliographic identities.
