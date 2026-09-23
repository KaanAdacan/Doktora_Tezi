#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

bash ./verify_package.sh
mkdir -p runs/01_NVT runs/02_NPT runs/03_PROD_000_025ns runs/04_PROD_025_050ns runs/05_PROD_050_075ns runs/06_PROD_075_100ns

if [[ -e .ADK_CHAIN_SUBMITTED && "${FORCE_SUBMIT:-0}" != "1" ]]; then
  echo "REFUSING_DUPLICATE_SUBMISSION: .ADK_CHAIN_SUBMITTED exists"
  echo "Inspect squeue -u $USER. If you intentionally need a second chain, run FORCE_SUBMIT=1 bash submit_chain.sh"
  exit 2
fi

J1=$(sbatch --parsable 01_NVT.batch)
J2=$(sbatch --parsable --dependency=afterok:$J1 02_NPT.batch)
J3=$(sbatch --parsable --dependency=afterok:$J2 03_PROD_000_025ns.batch)
J4=$(sbatch --parsable --dependency=afterok:$J3 04_PROD_025_050ns.batch)
J5=$(sbatch --parsable --dependency=afterok:$J4 05_PROD_050_075ns.batch)
J6=$(sbatch --parsable --dependency=afterok:$J5 06_PROD_075_100ns.batch)

STAMP=$(date +%Y%m%dT%H%M%S)
REPORT="ADK_CHAIN_${STAMP}.txt"
{
  echo "ADK_NVT=$J1"
  echo "ADK_NPT=$J2"
  echo "ADK_000_025ns=$J3"
  echo "ADK_025_050ns=$J4"
  echo "ADK_050_075ns=$J5"
  echo "ADK_075_100ns=$J6"
} | tee "$REPORT"
cp "$REPORT" .ADK_CHAIN_SUBMITTED

echo
echo "QUEUE_STATUS:"
squeue -u "$USER"
