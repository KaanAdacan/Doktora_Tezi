import numpy as np
from common import validate, reference, trajectory, SAMPLE_INDICES, OUTPUT, whole, fit, save_csv

validate()
ref,ref_protein,ref_bb,ref_ca=reference()
ref_ca_pos=ref_ca.positions.copy().astype(float)
u=trajectory()
protein=u.select_atoms("protein")
bb=u.select_atoms("protein and backbone")
ca=u.select_atoms("protein and name CA")
sum_d=np.zeros(ca.n_atoms,dtype=float)
sum_d2=np.zeros(ca.n_atoms,dtype=float)
max_d=np.zeros(ca.n_atoms,dtype=float)
n=0
for k,idx in enumerate(SAMPLE_INDICES,start=1):
    u.trajectory[idx]
    whole(protein)
    fit(bb,ref_bb,protein)
    d=np.linalg.norm(ca.positions.astype(float)-ref_ca_pos,axis=1)
    sum_d += d
    sum_d2 += d*d
    max_d=np.maximum(max_d,d)
    n += 1
    if k % 250 == 0:
        print(f"Residue dynamics {k} frames",flush=True)
rows=[]
for i,a in enumerate(ca):
    rows.append([a.segid,a.resid,a.resname,int(a.resindex),float(sum_d[i]/n),float(np.sqrt(sum_d2[i]/n)),float(max_d[i])])
save_csv(OUTPUT/"residue_dynamics.csv",["segid","resid","resname","resindex","mean_ca_displacement_from_0ns_A","rms_ca_displacement_from_0ns_A","max_ca_displacement_from_0ns_A"],rows)
print("Residue dynamics PASS")
