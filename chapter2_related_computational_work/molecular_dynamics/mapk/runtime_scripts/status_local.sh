#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"

echo "===== MAPK LOCAL CHAIN ====="
if [[ -s .MAPK_CHAIN.pid ]]; then
  pid=$(cat .MAPK_CHAIN.pid 2>/dev/null || true)
  if [[ -n "${pid:-}" ]] && kill -0 "$pid" 2>/dev/null; then
    echo "Launcher: RUNNING pid=$pid"
  else
    echo "Launcher: not running (stale pid=${pid:-unknown})"
  fi
else
  echo "Launcher: no pid file"
fi

echo
echo "===== ACTIVE NAMD ====="
pgrep -af '[n]amd3.*(01_NVT|02_NPT|03_PROD_000_025ns|04_PROD_025_050ns|05_PROD_050_075ns|06_PROD_075_100ns)\.conf' || echo "Aktif MAPK NAMD bulunamadi."

CONF=$(pgrep -af '[n]amd3.*(01_NVT|02_NPT|03_PROD_000_025ns|04_PROD_025_050ns|05_PROD_050_075ns|06_PROD_075_100ns)\.conf' \
  | grep -oE '(01_NVT|02_NPT|03_PROD_000_025ns|04_PROD_025_050ns|05_PROD_050_075ns|06_PROD_075_100ns)\.conf' | head -1 || true)
if [[ -n "$CONF" ]]; then
  stage="${CONF%.conf}"
  log="queue_logs/${stage}.log"
  echo; echo "===== ACTIVE CONFIG ====="; echo "$CONF"
  echo; echo "===== TIMING ====="; grep '^TIMING:' "$log" | tail -3 2>/dev/null || true
  echo; echo "===== ENERGY ====="; grep '^ENERGY:' "$log" | tail -1 2>/dev/null || true
fi

echo; echo "===== QUEUE ====="
tail -12 queue_logs/resume_master.log 2>/dev/null || echo "resume_master.log henuz yok."

echo; echo "===== GPU ====="
nvidia-smi --query-gpu=temperature.gpu,utilization.gpu,memory.used,memory.total,power.draw --format=csv,noheader 2>/dev/null || true

echo; echo "===== DISK ====="
du -sh runs 2>/dev/null || true
df -h . | tail -1
