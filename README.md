# Doktora_Tezi — Thesis code and reproducibility materials

Repository for source code and reproducibility materials associated with the PhD thesis of **Kaan Adacan**.

Canonical repository URL: https://github.com/KaanAdacan/Doktora_Tezi

This repository separates recovered historical analysis code, later provenance/QC gates, verified upstream pipeline provenance, and related computational work so that different code roles are not conflated.

## Chapter 1 — RNA-seq and miRNA-seq

- `chapter1_omics/executed_analysis/r_scripts/`
  - recovered final downstream RNA-seq v15 R analysis
  - recovered final downstream miRNA-seq v9 R analysis
- `chapter1_omics/r_scripts/`
  - sample-provenance validation and gate scripts
- `chapter1_omics/upstream/recovered_exact/`
  - reserved only for byte-identical upstream historical scripts that can be recovered
- `chapter1_omics/upstream/verified_provenance/`
  - STAR/Bowtie/featureCounts/count-matrix provenance supported by retained analysis records

The historical `seq.r` working launcher is not presented as a self-contained final pipeline.

## Chapter 3 — Cyx-KA CohortMaster and Cyx-KA Vesseller

Source code is intentionally withheld until publication of the associated scientific work. See `chapter3_tools/README.md`.

## Molecular-dynamics and related computational work

MD utilities are kept outside the thesis chapter numbering to avoid incorrectly implying that they constitute thesis Chapter 2:

- `related_computational_work/molecular_dynamics/mapk/`
- `related_computational_work/molecular_dynamics/adk/`

These directories include NAMD configurations/runtime utilities and analysis scripts for RMSD, RMSF, radius of gyration (Rg), residue dynamics, and snapshots.

## Scope and data policy

This repository does not distribute raw FASTQ/BAM files, count matrices, MD trajectories, coordinate/topology data, or other large experimental/source datasets unless a later release explicitly states otherwise.

## Integrity

See:

- `QC_REPORT.txt`
- `SHA256SUMS.txt`
- `REPOSITORY_FILE_MANIFEST_SHA256.tsv`
- `chapter1_omics/CODE_MANIFEST.tsv`
- `chapter1_omics/executed_analysis/PROVENANCE_SHA256.txt`

For thesis citation, use a frozen GitHub release/tag and commit SHA rather than the moving default branch.

## License

No open-source license has been selected in this release. Source-code reuse rights should be set deliberately before publication.
