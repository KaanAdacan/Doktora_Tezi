#!/usr/bin/env python3
from pathlib import Path
import tempfile,struct,subprocess,sys,json,csv
HERE=Path(__file__).resolve(); SCRIPT=HERE.parents[1]/"code"/"topology_ligand_snapshots.py"; REFS=HERE.parents[1]/"references"/"REFERENCES_VERIFIED.tsv"
def pdb_line(serial,name,res,chain,resid,x,y,z,segid,element="C"):
    return f"ATOM  {serial:5d} {name:>4s} {res[:3]:>3s} {chain:1s}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      {segid:<4s}{element:>2s}"
def write_pair(base):
    atoms=[("PROA","1","ALA","CA","CT1"),("PROA","2","ALA","CA","CT1"),("PROA","3","GLY","CA","CT2"),("PROA","4","SER","CA","CT1"),
           ("LIGA","900","COR","C1","CG3R"),("LIGA","900","COR","C2","CG3R"),("LIGA","900","COR","C3","CG3R"),("LIGA","900","COR","C4","CG3R"),("LIGA","900","COR","C5","CG3R"),("LIGA","900","COR","C6","CG3R"),
           ("WT1","2","TIP3","OH2","OT"),("WT1","2","TIP3","H1","HT"),("WT1","2","TIP3","H2","HT")]
    base.mkdir(parents=True,exist_ok=True); psf=base/"system.psf"; pdb=base/"system.pdb"
    lines=["PSF","",f"{len(atoms):8d} !NATOM"]; pl=[]
    for i,(seg,resid,res,an,typ) in enumerate(atoms,1):
        lines.append(f"{i:8d} {seg:<4s} {resid:>4s} {res:<4s} {an:<4s} {typ:<6s} 0.000000 12.0110 0")
        pl.append(pdb_line(i,an,res,"A",int(resid),float(i),float(i+1),float(i+2),seg))
    psf.write_text("\n".join(lines)+"\n"); pdb.write_text("\n".join(pl)+"\nEND\n"); return psf,pdb,len(atoms)
def coor(p,n,offset):
    vals=[]
    for i in range(n): vals.extend([offset+i,offset+i+0.1,offset+i+0.2])
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("wb") as f: f.write(struct.pack("<i",n)); f.write(struct.pack("<"+f"{n*3}d",*vals))
with tempfile.TemporaryDirectory() as td:
    td=Path(td); root=td/"systems"; s1=td/"stage01"; out=td/"out"; s1.mkdir(); state_rows=[]
    for target in ("ADK","MAPK"):
        active=root/target/"NVT_NPT_GPU"; _,_,n=write_pair(active/"inputs")
        (active/"04_PROD_025_050ns.conf").write_text("structure inputs/system.psf\ncoordinates inputs/system.pdb\n")
        arch=root/target/"ORIGINAL"/target; write_pair(arch/"inputs"); (arch/"04_PROD_025_050ns.conf").write_text("structure inputs/system.psf\ncoordinates inputs/system.pdb\n")
        cfgs={"T1":f"{target}/NVT_NPT_GPU/03_PROD_000_025ns.conf","T2":f"{target}/NVT_NPT_GPU/05_PROD_050_075ns.conf","T3":f"{target}/NVT_NPT_GPU/06_PROD_075_100ns.conf"}
        for rel in cfgs.values(): (root/rel).write_text("# active production config\n")
        for lab,ns,off in [("T1",0,0.0),("T2",50,50.0),("T3",100,100.0)]:
            p=active/"runs"/f"{lab}.coor"; coor(p,n,off)
            state_rows.append({"target":target,"label":lab,"time_ns":float(ns),"status":"RESOLVED","config":cfgs[lab],"source_kind":"","state_path":str(p.relative_to(root)),"state_exists":"True","state_natom":n,"expected_natom":n,"sha256":"","candidate_count":1})
    fields=["target","label","time_ns","status","config","source_kind","state_path","state_exists","state_natom","expected_natom","sha256","candidate_count"]
    with (s1/"T1_T2_T3_RESOLUTION.tsv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(state_rows)
    (s1/"STATE_DISCOVERY_SUMMARY.json").write_text(json.dumps({"final_status":"PASS","snapshot_extraction_eligible":True}))
    cmd=[sys.executable,str(SCRIPT),"--root",str(root),"--stage01",str(s1),"--out",str(out),"--references",str(REFS)]
    r1=subprocess.run(cmd,capture_output=True,text=True)
    if not (out/"STAGE02_SUMMARY.json").exists(): print(r1.stdout); print(r1.stderr); raise SystemExit(2)
    summ=json.loads((out/"STAGE02_SUMMARY.json").read_text()); snaps=list(csv.DictReader((out/"SNAPSHOT_MANIFEST.tsv").open(),delimiter="\t")); lig=list(csv.DictReader((out/"LIGAND_CANDIDATES.tsv").open(),delimiter="\t")); audit=list(csv.DictReader((out/"TOPOLOGY_PAIR_AUDIT.tsv").open(),delimiter="\t"))
    out2=td/"out2"; r2=subprocess.run([sys.executable,str(SCRIPT),"--root",str(root),"--stage01",str(s1),"--out",str(out2),"--references",str(REFS)],capture_output=True,text=True); snaps2=list(csv.DictReader((out2/"SNAPSHOT_MANIFEST.tsv").open(),delimiter="\t"))
    key=lambda rr:[(x["target"],x["label"],x["time_ns"],x["natom"],x["endianness"],x["sha256"]) for x in rr]
    checks={"V1_run1_exit_zero":r1.returncode==0,"V2_run2_exit_zero":r2.returncode==0,"V3_all_10_pass":all(x["status"]=="PASS" for x in summ["passes"]),"V4_TIP3_to_TIP_signature_normalized":all(x["signature_match"]=="True" for x in audit),"V5_active_lineage_selected":sum(x["selected_active_pair"]=="True" for x in audit)==2 and all("NVT_NPT_GPU" in x["config"] for x in audit if x["selected_active_pair"]=="True"),"V6_archive_not_selected":all(x["selected_active_pair"]=="False" for x in audit if "ORIGINAL" in x["config"]),"V7_unique_ligand_each_target":len(lig)==2 and all(x["resname"]=="COR" for x in lig),"V8_six_snapshots_valid":len(snaps)==6 and all(x["topology_signature_match"]=="True" for x in snaps),"V9_deterministic_within_target_hashes":key(snaps)==key(snaps2) and all(len({x["sha256"] for x in snaps if x["target"]==t})==3 for t in ("ADK","MAPK")),"V10_final_eligible":summ["rmsd_stage_eligible"] is True}
    rep=HERE.parent/"KNOWN_ANSWER_VERIFICATION.txt"
    with rep.open("w") as f:
        f.write("02_TOPOLOGY_LIGAND_SNAPSHOTS KNOWN-ANSWER VERIFICATION — LINEAGE/TIP3 AWARE\n"+"="*92+"\n")
        for k,v in checks.items(): f.write(f"{k:<58} {'PASS' if v else 'FAIL'}\n")
        f.write("-"*92+"\nFINAL: "+("10/10 PASS" if all(checks.values()) else "FAIL")+"\n")
    print(rep.read_text())
    if not all(checks.values()): print(json.dumps(summ,indent=2)); print(audit); print(r1.stdout,r1.stderr); raise SystemExit(2)
