set -euo pipefail
source "$HOME/venvs/md-analysis/bin/activate"
python rmsd.py
python rmsf.py
python rg.py
python residue_dynamics.py
python snapshots.py
