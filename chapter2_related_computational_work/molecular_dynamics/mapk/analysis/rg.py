import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import validate, trajectory, SAMPLE_INDICES, OUTPUT, sample_time_ns, whole, save_csv

validate()
u=trajectory()
protein=u.select_atoms("protein")
rows=[]
for n,idx in enumerate(SAMPLE_INDICES,start=1):
    u.trajectory[idx]
    whole(protein)
    rows.append([idx,sample_time_ns(idx),float(protein.radius_of_gyration())])
    if n % 250 == 0:
        print(f"Rg {n} frames {sample_time_ns(idx):.3f} ns",flush=True)
save_csv(OUTPUT/"rg_timeseries.csv",["frame_index_0based","time_ns","protein_rg_A"],rows)
a=np.asarray(rows,dtype=float)
plt.figure(figsize=(8,4.8))
plt.plot(a[:,1],a[:,2])
plt.xlabel("Time (ns)")
plt.ylabel("Protein Rg (Å)")
plt.tight_layout()
plt.savefig(OUTPUT/"rg.png",dpi=300)
plt.close()
print("Rg PASS")
