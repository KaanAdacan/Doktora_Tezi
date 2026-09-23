import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import validate, reference, trajectory, SAMPLE_INDICES, OUTPUT, sample_time_ns, whole, fit, calc_rmsd

validate()
ref,ref_protein,ref_bb,ref_ca=reference()
ref_bb_pos=ref_bb.positions.copy()
ref_ca_pos=ref_ca.positions.copy()
u=trajectory()
protein=u.select_atoms("protein")
bb=u.select_atoms("protein and backbone")
ca=u.select_atoms("protein and name CA")
rows=[]
for n,idx in enumerate(SAMPLE_INDICES,start=1):
    u.trajectory[idx]
    whole(protein)
    fit(bb,ref_bb,protein)
    rows.append([idx,sample_time_ns(idx),float(calc_rmsd(ca.positions,ref_ca_pos,center=False,superposition=False)),float(calc_rmsd(bb.positions,ref_bb_pos,center=False,superposition=False))])
    if n % 250 == 0:
        print(f"RMSD {n} frames {sample_time_ns(idx):.3f} ns",flush=True)
from common import save_csv
save_csv(OUTPUT/"rmsd_timeseries.csv",["frame_index_0based","time_ns","ca_rmsd_A","backbone_rmsd_A"],rows)
a=np.asarray(rows,dtype=float)
plt.figure(figsize=(8,4.8))
plt.plot(a[:,1],a[:,2],label="Cα")
plt.plot(a[:,1],a[:,3],label="Backbone")
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (Å)")
plt.legend()
plt.tight_layout()
plt.savefig(OUTPUT/"rmsd.png",dpi=300)
plt.close()
print("RMSD PASS")
