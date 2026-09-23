import warnings
import MDAnalysis as mda
from common import validate, reference, trajectory, OUTPUT, SYSTEM, EXPECTED_TOTAL_FRAMES, RAW_FRAME_PS, whole, fit

validate()
ref,ref_protein,ref_bb,ref_ca=reference()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ref_protein.write(str(OUTPUT/f"{SYSTEM}_000ns_protein.pdb"))
u=trajectory()
protein=u.select_atoms("protein")
bb=u.select_atoms("protein and backbone")
for ns in (50,100):
    idx=int(round(ns*1000.0/RAW_FRAME_PS))-1
    if idx < 0 or idx >= EXPECTED_TOTAL_FRAMES:
        raise SystemExit(f"FATAL: target frame outside trajectory for {ns} ns")
    u.trajectory[idx]
    whole(protein)
    fit(bb,ref_bb,protein)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        protein.write(str(OUTPUT/f"{SYSTEM}_{ns:03d}ns_protein.pdb"))
    print(f"Snapshot {ns} ns PASS",flush=True)
print("Snapshots PASS")
