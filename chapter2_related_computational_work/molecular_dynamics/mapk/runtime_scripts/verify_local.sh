#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

REQ=(
  inputs/mapk_ion.pdb
  inputs/mapk_ion.psf
  inputs/par_water_ions_ale.prm
  inputs/par_all36_prot.prm
  inputs/par_all36_cgenff.prm
  01_NVT.conf 02_NPT.conf
  03_PROD_000_025ns.conf 04_PROD_025_050ns.conf
  05_PROD_050_075ns.conf 06_PROD_075_100ns.conf
)
for f in "${REQ[@]}"; do
  [[ -s "$f" ]] || { echo "MISSING=$f"; exit 1; }
done

NAMD_BIN="${NAMD_BIN:-$(command -v namd3 || true)}"
[[ -n "$NAMD_BIN" && -x "$NAMD_BIN" ]] || {
  echo "FATAL: namd3 bulunamadi. ADK'de kullandigin NAMD 3 CUDA binary PATH icinde olmali." >&2
  exit 2
}
command -v nvidia-smi >/dev/null 2>&1 || { echo "FATAL: nvidia-smi bulunamadi" >&2; exit 3; }

for f in run_chain_local.sh start_mapk.sh status_local.sh watch_mapk.sh stop_mapk.sh verify_local.sh; do
  bash -n "$f"
done

atoms=$(awk 'BEGIN{n=0} /^ATOM  |^HETATM/{n++} END{print n}' inputs/mapk_ion.pdb)
echo "PACKAGE_STATIC_CHECK=PASS"
echo "PWD=$PWD"
echo "NAMD_BIN=$NAMD_BIN"
echo "PDB_ATOMS=$atoms"
echo "GPU:"
nvidia-smi -L
echo "CONFIG_SUMMARY:"
for f in 0*.conf; do
  printf '%-28s ' "$f"
  grep -E '^(run|minimize|dcdfreq|restartfreq|GPUresident)' "$f" | tr '\n' '; '
  echo
done
