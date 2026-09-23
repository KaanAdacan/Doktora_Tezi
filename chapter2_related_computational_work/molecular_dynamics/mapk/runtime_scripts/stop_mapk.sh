#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
echo "Bu komut MAPK NAMD prosesini sonlandirir; mevcut stage daha sonra temizden tekrar kosar."

pids=$(pgrep -f '[n]amd3.*(01_NVT|02_NPT|03_PROD_000_025ns|04_PROD_025_050ns|05_PROD_050_075ns|06_PROD_075_100ns)\.conf' || true)
if [[ -n "$pids" ]]; then
  echo "Stopping NAMD PID(s): $pids"
  kill $pids || true
fi
if [[ -s .MAPK_CHAIN.pid ]]; then
  pid=$(cat .MAPK_CHAIN.pid || true)
  [[ -n "${pid:-}" ]] && kill "$pid" 2>/dev/null || true
fi
echo "STOP signal sent."
