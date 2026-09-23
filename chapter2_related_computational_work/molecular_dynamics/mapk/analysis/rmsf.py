import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import validate, reference, trajectory, SAMPLE_INDICES, OUTPUT, whole, fit, save_csv

validate()
ref,ref_protein,ref_bb,ref_ca=reference()
u=trajectory()
protein=u.select_atoms("protein")
bb=u.select_atoms("protein and backbone")
ca=u.select_atoms("protein and name CA")
ca_mean=np.zeros((ca.n_atoms,3),dtype=float)
ca_m2=np.zeros_like(ca_mean)
bb_mean=np.zeros((bb.n_atoms,3),dtype=float)
bb_m2=np.zeros_like(bb_mean)
n=0
for k,idx in enumerate(SAMPLE_INDICES,start=1):
    u.trajectory[idx]
    whole(protein)
    fit(bb,ref_bb,protein)
    n += 1
    x=ca.positions.astype(float)
    d=x-ca_mean
    ca_mean += d/n
    ca_m2 += d*(x-ca_mean)
    xb=bb.positions.astype(float)
    db=xb-bb_mean
    bb_mean += db/n
    bb_m2 += db*(xb-bb_mean)
    if k % 250 == 0:
        print(f"RMSF {k} frames",flush=True)
ca_rmsf=np.sqrt(np.sum(ca_m2/n,axis=1))
bb_atom_rmsf=np.sqrt(np.sum(bb_m2/n,axis=1))
by_res={}
for value,ridx in zip(bb_atom_rmsf,bb.resindices):
    by_res.setdefault(int(ridx),[]).append(float(value))
rows=[]
for i,a in enumerate(ca):
    rows.append([a.segid,a.resid,a.resname,int(a.resindex),float(ca_rmsf[i]),float(np.mean(by_res.get(int(a.resindex),[np.nan])))])
save_csv(OUTPUT/"rmsf_residue.csv",["segid","resid","resname","resindex","ca_rmsf_A","mean_backbone_atom_rmsf_A"],rows)
plt.figure(figsize=(9,4.8))
plt.plot(np.arange(ca.n_atoms),ca_rmsf)
plt.xlabel("Protein residue order")
plt.ylabel("Cα RMSF (Å)")
plt.tight_layout()
plt.savefig(OUTPUT/"rmsf.png",dpi=300)
plt.close()
print("RMSF PASS")
