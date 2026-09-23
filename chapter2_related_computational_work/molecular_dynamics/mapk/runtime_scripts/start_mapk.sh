#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p queue_logs
bash ./verify_local.sh

if [[ -s .MAPK_CHAIN.pid ]]; then
  oldpid=$(cat .MAPK_CHAIN.pid 2>/dev/null || true)
  if [[ -n "${oldpid:-}" ]] && kill -0 "$oldpid" 2>/dev/null; then
    echo "MAPK chain zaten calisiyor. PID=$oldpid"
    exit 2
  fi
fi

PE="${PE:-8}" GPU_DEVICE="${GPU_DEVICE:-0}" \
  nohup bash ./run_chain_local.sh >> queue_logs/launcher.log 2>&1 &
pid=$!
echo "$pid" > .MAPK_CHAIN.pid
echo "MAPK_LOCAL_CHAIN_STARTED PID=$pid PE=${PE:-8} GPU_DEVICE=${GPU_DEVICE:-0}"
echo "Durum: bash status_local.sh"
echo "Canli izleme: bash watch_mapk.sh"
