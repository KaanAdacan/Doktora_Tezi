#!/bin/bash
set -u
cd "$(dirname "$0")"
echo "===== SQUEUE ====="
squeue -u "$USER"
echo
echo "===== RECENT ADK SACCT ====="
sacct -S now-7days -u "$USER" --name=ADK_01_NVT,ADK_02_NPT,ADK_03_P000_025,ADK_04_P025_050,ADK_05_P050_075,ADK_06_P075_100 \
  --format=JobID,JobName%18,Partition%12,State%18,ExitCode,Elapsed,NodeList -X 2>/dev/null || true
echo
echo "===== FINAL STATES ====="
find runs -maxdepth 2 -type f \( -name '*.coor' -o -name '*.vel' -o -name '*.xsc' \) -printf '%p %s bytes\n' 2>/dev/null | sort || true
