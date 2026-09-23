# NAMD Cordycepin MAPK/ADK — matched 10-ps analysis pipeline v2.2

This directory contains the analysis code used for the final matched-temporal-resolution Chapter 2 MAPK/ADK molecular-dynamics workflow.

## Sampling policy
- MAPK raw trajectory: 10 ps/frame; analysis stride = 1; 10,000 production frames.
- ADK raw trajectory: 2 ps/frame; analysis stride = 5; 10,000 analysis frames.
- Common analysis grid: 0.01–100.00 ns at 10 ps spacing.
- Production only; NVT/NPT are excluded from reported structural observables.
- Raw trajectory files are read-only. No interpolation, smoothing, filtering, imputation, outlier removal, or frame balancing is applied.

## Interpretation lock
The verified production systems are receptor-only/apo trajectories. These trajectories support receptor conformational sampling and do not constitute ligand-bound MD. T1/T2/T3 docking states are exact 0/50/100 ns receptor boundary states and are docked independently.

## Reproducibility
All stages include verification code. The final real-data matched-10-ps Stage 03 and Stage 04 runs passed their internal 10-PASS checks. Generated analysis outputs and multi-GB raw DCD trajectories are not committed in this source-code directory.
