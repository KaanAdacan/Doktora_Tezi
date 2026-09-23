#!/usr/bin/env bash
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV="${HOME}/NAMD/envs/cordy_md_analysis"
S3="/home/kaan/NAMD/Publication_Analysis/03_APO_PROTEIN_METRICS_MATCHED10PS"
S4="/home/kaan/NAMD/Publication_Analysis/04_PCA_CLUSTERING_MATCHED10PS"
S5="/home/kaan/NAMD/Publication_Analysis/05_THESIS_FIGURES_MATCHED10PS"
[[ "${WSL_DISTRO_NAME:-}" == "NAMD-Ubuntu" ]] || { echo "[FAIL] Run inside NAMD-Ubuntu"; exit 2; }
[[ -x "$ENV/bin/python" ]] || { echo "[FAIL] Missing env: $ENV"; exit 2; }

echo "=== 0. Known-answer sampling verification ==="
"$ENV/bin/python" "$SRC/03_APO_PROTEIN_METRICS_MATCHED10PS/verification/verify_matched10ps.py"

echo "=== 1. Stage 03 matched 10 ps ==="
rm -rf "$S3"
"$ENV/bin/python" "$SRC/03_APO_PROTEIN_METRICS_MATCHED10PS/code/apo_protein_metrics_matched10ps.py" \
 --root /home/kaan/NAMD/systems \
 --qc /home/kaan/NAMD/Publication_Analysis/00_QC \
 --stage02 /home/kaan/NAMD/Publication_Analysis/02_TOPOLOGY_LIGAND_SNAPSHOTS \
 --stage02a /home/kaan/NAMD/Publication_Analysis/02A_TOPOLOGY_MODE_AUDIT \
 --out "$S3" \
 --references "$SRC/03_APO_PROTEIN_METRICS_MATCHED10PS/references/REFERENCES_VERIFIED.tsv"
cat "$S3/10PASS_REPORT.txt"

echo "=== 2. Stage 04 matched 10 ps ==="
rm -rf "$S4"
"$ENV/bin/python" "$SRC/04_PCA_CLUSTERING_MATCHED10PS/code/pca_clustering_matched10ps.py" \
 --root /home/kaan/NAMD/systems \
 --qc /home/kaan/NAMD/Publication_Analysis/00_QC \
 --stage02 /home/kaan/NAMD/Publication_Analysis/02_TOPOLOGY_LIGAND_SNAPSHOTS \
 --stage02a /home/kaan/NAMD/Publication_Analysis/02A_TOPOLOGY_MODE_AUDIT \
 --stage03 "$S3" --out "$S4" \
 --references "$SRC/04_PCA_CLUSTERING_MATCHED10PS/references/REFERENCES_VERIFIED.tsv"
cat "$S4/10PASS_REPORT.txt"

echo "=== 3. Thesis figures ==="
rm -rf "$S5"; mkdir -p "$S5"
"$ENV/bin/python" "$SRC/05_THESIS_FIGURES_MATCHED10PS/code/build_thesis_figures_matched10ps.py" --stage03 "$S3" --stage04 "$S4" --out "$S5"

echo "=== 4. Updated method/results/discussion ==="
"$ENV/bin/python" "$SRC/05_THESIS_FIGURES_MATCHED10PS/code/generate_updated_text.py" --stage03 "$S3" --stage04 "$S4" --out "$S5/UPDATED_METHODS_RESULTS_DISCUSSION_TR.md"

echo "=== 5. Final real-results zip ==="
DEST="/mnt/c/Users/kaana/Downloads/NAMD_Cordycepin_MATCHED10PS_REAL_RESULTS"
rm -rf "$DEST" "$DEST.zip"; mkdir -p "$DEST"
cp -a "$S3" "$DEST/"; cp -a "$S4" "$DEST/"; cp -a "$S5" "$DEST/"; cp -a "$SRC/METHODS_RESULTS_DISCUSSION_UPDATE_TR.md" "$DEST/"
python3 - <<'PY'
from pathlib import Path
import hashlib,csv,zipfile,os
root=Path('/mnt/c/Users/kaana/Downloads/NAMD_Cordycepin_MATCHED10PS_REAL_RESULTS')
rows=[]
for p in sorted(root.rglob('*')):
    if p.is_file(): rows.append([str(p.relative_to(root)),hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_size])
with (root/'SHA256_MANIFEST.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['relative_path','sha256','bytes']);w.writerows(rows)
z=Path(str(root)+'.zip')
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as zz:
    for p in root.rglob('*'):
        if p.is_file():zz.write(p,arcname=str(root.name+'/'+str(p.relative_to(root))))
print('[PASS]',z)
PY

echo "FINAL ZIP: /mnt/c/Users/kaana/Downloads/NAMD_Cordycepin_MATCHED10PS_REAL_RESULTS.zip"
