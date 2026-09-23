#!/usr/bin/env bash
set -euo pipefail

ENV="${HOME}/NAMD/envs/cordy_md_analysis"
PYTHON="${PYTHON:-python3}"

echo "======================================================================"
echo "CORDYCEPIN MD ANALYSIS ENVIRONMENT SETUP"
echo "======================================================================"
echo "Python: $($PYTHON --version 2>&1)"
echo "Target: $ENV"

if ! "$PYTHON" -m venv --help >/dev/null 2>&1; then
  echo "[FAIL] Python venv module is unavailable."
  echo "On Ubuntu, install python3-venv for this Python version, then rerun."
  exit 2
fi

rm -rf "$ENV"
"$PYTHON" -m venv "$ENV"
"$ENV/bin/python" -m pip install --upgrade pip setuptools wheel
"$ENV/bin/python" -m pip install "MDAnalysis[analysis]==2.10.0" matplotlib

echo
echo "=== IMPORT VERIFICATION ==="
"$ENV/bin/python" - <<'PY'
import MDAnalysis, numpy, matplotlib
print("[PASS] MDAnalysis", MDAnalysis.__version__)
print("[PASS] numpy", numpy.__version__)
print("[PASS] matplotlib", matplotlib.__version__)
assert MDAnalysis.__version__ == "2.10.0"
PY

"$ENV/bin/python" -m pip freeze | sort > "${ENV}/ENV_FREEZE.txt"
sha256sum "${ENV}/ENV_FREEZE.txt" > "${ENV}/ENV_FREEZE.sha256"

echo
echo "[PASS] Environment ready: $ENV"
echo "[PASS] Freeze: $ENV/ENV_FREEZE.txt"
