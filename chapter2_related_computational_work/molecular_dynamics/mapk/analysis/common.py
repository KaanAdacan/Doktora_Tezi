from pathlib import Path
import math
import csv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd as calc_rmsd

ROOT = Path('/home/kaan/NAMD/systems/MAPK/MAPK_NVT_NPT_GPU')
PSF = ROOT / 'inputs/mapk_ion.psf'
NPT_COOR = ROOT / 'runs/02_NPT/MAPK_NPT.coor'
NPT_XSC = ROOT / 'runs/02_NPT/MAPK_NPT.xsc'
DCDS = [ROOT / 'runs/03_PROD_000_025ns/MAPK_PROD_000_025ns.dcd',
    ROOT / 'runs/04_PROD_025_050ns/MAPK_PROD_025_050ns.dcd',
    ROOT / 'runs/05_PROD_050_075ns/MAPK_PROD_050_075ns.dcd',
    ROOT / 'runs/06_PROD_075_100ns/MAPK_PROD_075_100ns.dcd',
    ROOT / 'runs/06_PROD_075_100ns_RESUME_20260920_102439/MAPK_PROD_075_100ns_RESUME.dcd']
EXPECTED_ATOMS = 181249
EXPECTED_SEGMENT_FRAMES = [2500, 2500, 2500, 465, 2035]
EXPECTED_TOTAL_FRAMES = 10000
RAW_FRAME_PS = 10.0
SAMPLE_INDICES = range(0, 10000, 1)
SYSTEM = 'MAPK'
OUTPUT = ROOT / "analysis_0_100ns"


def die(message):
    raise SystemExit("FATAL: " + str(message))


def xsc_dimensions(path):
    if not path.exists():
        return None
    last = None
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s[:1] == chr(35):
            continue
        parts = s.split()
        if len(parts) >= 13:
            try:
                vals = [float(x) for x in parts[:13]]
                last = vals
            except ValueError:
                pass
    if last is None:
        return None
    a = np.array(last[1:4], dtype=float)
    b = np.array(last[4:7], dtype=float)
    c = np.array(last[7:10], dtype=float)
    def angle(u, v):
        z = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return math.degrees(math.acos(float(np.clip(z, -1.0, 1.0))))
    return np.array([np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c), angle(b,c), angle(a,c), angle(a,b)], dtype=np.float32)


def validate():
    required = [PSF, NPT_COOR, *DCDS]
    for path in required:
        if not path.exists():
            die("Missing file: " + str(path))
    rows=[]
    total=0
    for i,(dcd,expected) in enumerate(zip(DCDS,EXPECTED_SEGMENT_FRAMES),start=1):
        u=mda.Universe(str(PSF),str(dcd))
        atoms=u.atoms.n_atoms
        frames=len(u.trajectory)
        if atoms != EXPECTED_ATOMS:
            die(f"Unexpected atom count in {dcd}: {atoms}")
        if frames != expected:
            die(f"Unexpected frame count in {dcd}: {frames}; expected {expected}")
        total += frames
        rows.append([i,str(dcd),atoms,frames,"PASS"])
    if total != EXPECTED_TOTAL_FRAMES:
        die(f"Total frames {total}; expected {EXPECTED_TOTAL_FRAMES}")
    OUTPUT.mkdir(parents=True,exist_ok=True)
    with (OUTPUT/'qc_segments.csv').open('w',newline='') as fh:
        w=csv.writer(fh)
        w.writerow(['segment','dcd','n_atoms','n_frames','status'])
        w.writerows(rows)
    return total


def reference():
    ref=mda.Universe(str(PSF),str(NPT_COOR))
    dims=xsc_dimensions(NPT_XSC)
    if dims is not None:
        ref.dimensions=dims
    protein=ref.select_atoms('protein')
    backbone=ref.select_atoms('protein and backbone')
    ca=ref.select_atoms('protein and name CA')
    if not protein.n_atoms or not backbone.n_atoms or not ca.n_atoms:
        die('Protein selection is empty')
    if dims is not None:
        protein.unwrap(compound='fragments')
    return ref,protein,backbone,ca


def trajectory():
    return mda.Universe(str(PSF),[str(x) for x in DCDS])


def sample_time_ns(frame_index):
    return (frame_index + 1) * RAW_FRAME_PS / 1000.0


def whole(protein):
    protein.unwrap(compound='fragments')


def fit(backbone, reference_backbone, protein):
    align.alignto(backbone,reference_backbone,select=None,weights='mass',subselection=protein,strict=True)


def save_csv(path,header,rows):
    with path.open('w',newline='') as fh:
        w=csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
