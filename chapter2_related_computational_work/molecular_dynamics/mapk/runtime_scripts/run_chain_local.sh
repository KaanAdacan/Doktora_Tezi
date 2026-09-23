#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

PE="${PE:-8}"
GPU_DEVICE="${GPU_DEVICE:-0}"
NAMD_BIN="${NAMD_BIN:-$(command -v namd3 || true)}"

[[ -n "$NAMD_BIN" && -x "$NAMD_BIN" ]] || { echo "FATAL: namd3 bulunamadi" >&2; exit 10; }
command -v nvidia-smi >/dev/null 2>&1 || { echo "FATAL: nvidia-smi bulunamadi" >&2; exit 11; }

mkdir -p queue_logs failed_runs \
  runs/01_NVT runs/02_NPT runs/03_PROD_000_025ns \
  runs/04_PROD_025_050ns runs/05_PROD_050_075ns runs/06_PROD_075_100ns

exec 9>.MAPK_LOCAL_CHAIN.lock
if ! flock -n 9; then
  echo "FATAL: MAPK local chain zaten calisiyor (lock aktif)." >&2
  exit 12
fi

MASTER="queue_logs/resume_master.log"
log_master() { echo "[$(date '+%F %T')] $*" | tee -a "$MASTER"; }

STAGES=(
  "01_NVT|01_NVT.conf|runs/01_NVT|MAPK_NVT"
  "02_NPT|02_NPT.conf|runs/02_NPT|MAPK_NPT"
  "03_PROD_000_025ns|03_PROD_000_025ns.conf|runs/03_PROD_000_025ns|MAPK_PROD_000_025ns"
  "04_PROD_025_050ns|04_PROD_025_050ns.conf|runs/04_PROD_025_050ns|MAPK_PROD_025_050ns"
  "05_PROD_050_075ns|05_PROD_050_075ns.conf|runs/05_PROD_050_075ns|MAPK_PROD_050_075ns"
  "06_PROD_075_100ns|06_PROD_075_100ns.conf|runs/06_PROD_075_100ns|MAPK_PROD_075_100ns"
)

stage_complete() {
  local stage="$1" dir="$2" prefix="$3"
  local log="queue_logs/${stage}.log"
  [[ -s "$dir/$prefix.coor" && -s "$dir/$prefix.vel" && -s "$dir/$prefix.xsc" ]] || return 1
  [[ -s "$log" ]] || return 1
  grep -q '\[Partition 0\]\[Node 0\] End of program' "$log" || return 1
}

archive_incomplete_stage() {
  local stage="$1" dir="$2"
  if find "$dir" -mindepth 1 -maxdepth 1 -type f -print -quit | grep -q .; then
    local stamp dest
    stamp=$(date '+%Y%m%dT%H%M%S')
    dest="failed_runs/${stage}_${stamp}"
    mkdir -p "$dest"
    find "$dir" -mindepth 1 -maxdepth 1 -type f -exec mv -t "$dest" {} +
    log_master "ARCHIVED incomplete $stage -> $dest"
  fi
}

for f in inputs/mapk_ion.psf inputs/mapk_ion.pdb \
         inputs/par_water_ions_ale.prm inputs/par_all36_prot.prm inputs/par_all36_cgenff.prm; do
  [[ -s "$f" ]] || { log_master "FATAL missing input: $f"; exit 20; }
done

log_master "=========================================="
log_master "MAPK LOCAL WSL GPU CHAIN START"
log_master "NAMD=$NAMD_BIN PE=+p$PE GPU_DEVICE=$GPU_DEVICE"
nvidia-smi -L | tee -a "$MASTER"

for rec in "${STAGES[@]}"; do
  IFS='|' read -r stage conf dir prefix <<< "$rec"
  log="queue_logs/${stage}.log"

  if stage_complete "$stage" "$dir" "$prefix"; then
    log_master "SKIP/PASS $stage (already complete)"
    continue
  fi

  archive_incomplete_stage "$stage" "$dir"
  rm -f "$log"

  log_master "START $stage"
  log_master "COMMAND: $NAMD_BIN +p$PE +setcpuaffinity +devices $GPU_DEVICE $conf"

  set +e
  "$NAMD_BIN" "+p${PE}" +setcpuaffinity +devices "$GPU_DEVICE" "$conf" > "$log" 2>&1
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    log_master "FAIL $stage exit=$rc"
    tail -40 "$log" | tee -a "$MASTER" || true
    exit "$rc"
  fi

  if ! stage_complete "$stage" "$dir" "$prefix"; then
    log_master "FAIL $stage: NAMD exited 0 but final state/completion marker missing"
    tail -40 "$log" | tee -a "$MASTER" || true
    exit 90
  fi

  perf=$(grep '^PERFORMANCE:' "$log" | tail -1 || true)
  timing=$(grep '^TIMING:' "$log" | tail -1 || true)
  [[ -n "$timing" ]] && log_master "$timing"
  [[ -n "$perf" ]] && log_master "$perf"
  log_master "DONE $stage"
done

log_master "=========================================="
log_master "ALL PASS: NVT + NPT + 100 ns PRODUCTION TAMAMLANDI"
log_master "PE: +p$PE"
last_perf=$(grep '^PERFORMANCE:' queue_logs/06_PROD_075_100ns.log | tail -1 || true)
[[ -n "$last_perf" ]] && log_master "FINAL $last_perf"
log_master "=========================================="
