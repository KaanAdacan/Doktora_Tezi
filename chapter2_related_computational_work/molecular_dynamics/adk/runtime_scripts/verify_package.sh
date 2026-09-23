#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

REQ=(
  inputs/ADK_ion.pdb
  inputs/ADK_ion.psf
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

for f in *.batch submit_chain.sh status.sh verify_package.sh; do
  bash -n "$f"
done

sha256sum -c SOURCE_INPUT_SHA256SUMS.txt >/dev/null

grep -Eq '^[[:space:]]*cellBasisVector1[[:space:]]+126\.16[[:space:]]+0[[:space:]]+0([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*cellBasisVector2[[:space:]]+0[[:space:]]+118\.25[[:space:]]+0([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*cellBasisVector3[[:space:]]+0[[:space:]]+0[[:space:]]+117\.03([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*cellOrigin[[:space:]]+73\.74[[:space:]]+76\.31[[:space:]]+78\.04([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*PMEGridSizeX[[:space:]]+128([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*PMEGridSizeY[[:space:]]+128([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*PMEGridSizeZ[[:space:]]+128([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*minimize[[:space:]]+50000([[:space:]]|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*run[[:space:]]+5000000([[:space:]]|;|$)' 01_NVT.conf
grep -Eq '^[[:space:]]*run[[:space:]]+5000000([[:space:]]|;|$)' 02_NPT.conf
for f in 03_PROD_000_025ns.conf 04_PROD_025_050ns.conf 05_PROD_050_075ns.conf 06_PROD_075_100ns.conf; do
  grep -Eq '^[[:space:]]*run[[:space:]]+12500000([[:space:]]|;|$)' "$f"
done

for f in 02_NPT.conf 03_PROD_000_025ns.conf 04_PROD_025_050ns.conf 05_PROD_050_075ns.conf 06_PROD_075_100ns.conf; do
  grep -q 'binCoordinates[[:space:]]\+\$inputname\.coor' "$f"
  grep -q 'binVelocities[[:space:]]\+\$inputname\.vel' "$f"
  grep -q 'extendedSystem[[:space:]]\+\$inputname\.xsc' "$f"
  grep -q 'firsttimestep[[:space:]]\+\[get_first_ts "\$inputname\.xsc"\]' "$f"
done

for f in *.batch; do
  grep -q '^#SBATCH --account=kadacan$' "$f"
  grep -q '^#SBATCH --partition=akya-cuda$' "$f"
  grep -q '^#SBATCH --cpus-per-task=10$' "$f"
  grep -q '^#SBATCH --gres=gpu:1$' "$f"
  grep -q '^#SBATCH --time=3-00:00:00$' "$f"
  grep -q 'module load apps/namd/3.0.1-multicore-CUDA' "$f"
  grep -q 'namd3 +p${SLURM_CPUS_PER_TASK} +setcpuaffinity +devices 0' "$f"
done

if grep -Eni 'mapk|partition=protein|account=compute|namd2|namd/2\.14' \
  01_NVT.conf 02_NPT.conf 03_PROD_000_025ns.conf 04_PROD_025_050ns.conf \
  05_PROD_050_075ns.conf 06_PROD_075_100ns.conf *.batch submit_chain.sh status.sh; then
  echo "FATAL: legacy/MAPK token found in executable package surface" >&2
  exit 30
fi

echo "PACKAGE_STATIC_CHECK=PASS"
echo "ADK_VECTOR_ORIGIN_CONTRACT=PASS"
echo "ADK_STEP_COUNT_CONTRACT=PASS"
echo "TRUBA_EXECUTION_CONTRACT=PASS"
echo "INPUT_SHA256_CONTRACT=PASS"
echo "PWD=$PWD"
echo "NOTE=This does not execute NAMD; compute-node CUDA validation occurs when the first job starts."
