#!/usr/bin/env python3
from pathlib import Path
import tempfile, subprocess, sys, json, csv, struct, hashlib, os

HERE=Path(__file__).resolve()
SCRIPT=HERE.parents[1]/"code"/"state_discovery.py"
REFS=HERE.parents[1]/"references"/"REFERENCES_VERIFIED.tsv"

def namdbin(p,n):
    with p.open("wb") as f:
        f.write(struct.pack("<i",n))
        for i in range(n*3): f.write(struct.pack("<d",float(i)))

def conf(p, bincoord, output, first, run, dt=2.0, dcdfreq=5000):
    p.write_text(f"""binCoordinates {bincoord}
outputName {output}
firsttimestep {first}
timestep {dt}
dcdFreq {dcdfreq}
run {run}
""")

with tempfile.TemporaryDirectory() as td:
    td=Path(td)
    root=td/"systems"; qc=td/"qc"; out=td/"out"; qc.mkdir()
    expected={"ADK":10,"MAPK":12}
    dcd_headers=[]; chain=[]
    for target in ("ADK","MAPK"):
        base=root/target/"RUN"; base.mkdir(parents=True); n=expected[target]
        namdbin(base/"state_000.coor",n);namdbin(base/"state_050.coor",n);namdbin(base/"state_100.coor",n)
        steps0=1000000 if target=="ADK" else 2000000
        for idx,(a,b) in enumerate([(0,25),(25,50),(50,75),(75,100)], start=3):
            d=base/f"{idx:02d}_PROD_{a:03d}_{b:03d}ns"; d.mkdir()
            if a==0: inp="../state_000.coor"
            elif a==50: inp="../state_050.coor"
            else:
                mid=base/f"state_{a:03d}.coor"
                if not mid.exists(): namdbin(mid,n)
                inp=f"../state_{a:03d}.coor"
            outprefix="../state_100" if b==100 else f"prod_{a:03d}_{b:03d}"
            conf(d/"run.conf",inp,outprefix,steps0+int(a*500000),int((b-a)*500000))
        for k,(a,b) in enumerate([(0,25),(25,50),(50,75),(75,100)]):
            rel=f"{target}/dummy_{a}_{b}.dcd";dcd_headers.append({"target":target,"relative_path":rel,"natom":n});chain.append({"target":target,"relative_path":rel})
    with (qc/"DCD_HEADERS.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["target","relative_path","natom"],delimiter="	"); w.writeheader(); w.writerows(dcd_headers)
    with (qc/"PRODUCTION_CHAIN.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["target","relative_path"],delimiter="	"); w.writeheader(); w.writerows(chain)
    (qc/"QC_SUMMARY.json").write_text(json.dumps({"final_status":"PASS","publication_eligible_for_next_stage":True}))
    cmd=[sys.executable,str(SCRIPT),"--root",str(root),"--qc",str(qc),"--out",str(out),"--references",str(REFS)]
    r1=subprocess.run(cmd,capture_output=True,text=True);s1=json.loads((out/"STATE_DISCOVERY_SUMMARY.json").read_text());res=list(csv.DictReader((out/"T1_T2_T3_RESOLUTION.tsv").open(),delimiter="	"))
    out2=td/"out2";cmd2=[sys.executable,str(SCRIPT),"--root",str(root),"--qc",str(qc),"--out",str(out2),"--references",str(REFS)];r2=subprocess.run(cmd2,capture_output=True,text=True);res2=list(csv.DictReader((out2/"T1_T2_T3_RESOLUTION.tsv").open(),delimiter="	"))
    key=lambda rows:[(x["target"],x["label"],x["time_ns"],x["status"],Path(x["state_path"]).name,x["state_natom"],x["sha256"]) for x in rows]
    checks={"V1_run1_exit_zero":r1.returncode==0,"V2_run2_exit_zero":r2.returncode==0,"V3_all_10_pass":all(x["status"]=="PASS" for x in s1["passes"]),"V4_six_states_resolved":len(res)==6 and all(x["status"]=="RESOLVED" for x in res),"V5_T1_is_exact_input_state":all(Path(x["state_path"]).name=="state_000.coor" for x in res if x["label"]=="T1"),"V6_T2_is_exact_input_state":all(Path(x["state_path"]).name=="state_050.coor" for x in res if x["label"]=="T2"),"V7_T3_is_final_output_state":all(Path(x["state_path"]).name=="state_100.coor" for x in res if x["label"]=="T3"),"V8_atomcounts_match":all(x["state_natom"]==x["expected_natom"] for x in res),"V9_hashes_and_determinism":key(res)==key(res2) and all(len(x["sha256"])==64 for x in res),"V10_final_eligible":s1["snapshot_extraction_eligible"] is True}
    report=HERE.parent/"KNOWN_ANSWER_VERIFICATION.txt"
    with report.open("w") as f:
        f.write("01_STATE_DISCOVERY KNOWN-ANSWER VERIFICATION\n"+"="*72+"\n")
        for k,v in checks.items():f.write(f"{k:<45} {'PASS' if v else 'FAIL'}\n")
        f.write("-"*72+"\nFINAL: "+("10/10 PASS" if all(checks.values()) else "FAIL")+"\n")
    print(report.read_text())
    if not all(checks.values()):print("RUN1\n",r1.stdout,r1.stderr);raise SystemExit(2)
