#!/usr/bin/env python3
"""
01_STATE_DISCOVERY
Resolve exact NAMD state sources for T1=0 ns, T2=50 ns, T3=100 ns.

Key rule:
NAMD DCD initial coordinates are not included, so exact boundary states are
resolved from NAMD binCoordinates/input states and final outputName.coor files,
not by relabeling nearest DCD frames.

No scientific observable is calculated here.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, re, struct, sys
from pathlib import Path
from typing import Dict, List, Optional

VERSION="1.0.0"

def sha256_file(p: Path) -> str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()

def read_tsv(p: Path):
    with p.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="	"))

def write_tsv(p: Path, rows, fields):
    with p.open("w", encoding="utf-8", newline="") as f:
        w=csv.DictWriter(f, fieldnames=fields, delimiter="	", extrasaction="ignore")
        w.writeheader()
        for r in rows: w.writerow({k:r.get(k,"") for k in fields})

def strip_comment(line: str) -> str:
    return line.split("#",1)[0].strip()

VAR1=re.compile(r"${([A-Za-z_][A-Za-z0-9_]*)}")
VAR2=re.compile(r"$([A-Za-z_][A-Za-z0-9_]*)")

def substitute(s: str, vars: Dict[str,str]) -> str:
    prev=None
    for _ in range(20):
        if s==prev: break
        prev=s
        s=VAR1.sub(lambda m: vars.get(m.group(1), m.group(0)), s)
        s=VAR2.sub(lambda m: vars.get(m.group(1), m.group(0)), s)
    return s.strip().strip('"').strip("'")

def parse_conf(p: Path) -> Dict[str,object]:
    vars={};directives={};runs=[]
    try: lines=p.read_text(errors="replace").splitlines()
    except Exception as e: return {"path":str(p),"error":str(e)}
    for raw in lines:
        line=strip_comment(raw)
        if not line: continue
        m=re.match(r"(?i)^sets+([A-Za-z_][A-Za-z0-9_]*)s+(.+?)s*$", line)
        if m: vars[m.group(1)] = substitute(m.group(2), vars)
    keys=("coordinates","bincoordinates","binvelocities","extendedsystem","outputname","restartname","cwd","firsttimestep","timestep","dcdfreq","restartfreq","binaryoutput")
    for raw in lines:
        line=strip_comment(raw)
        if not line: continue
        m=re.match(r"^(S+)s+(.+?)s*$", line)
        if not m: continue
        k=m.group(1).lower();v=substitute(m.group(2), vars)
        if k in keys: directives[k]=v
        elif k=="run":
            try:runs.append(int(v.split()[0]))
            except:pass
    directives["run_steps"]=runs;directives["path"]=str(p);directives["variables"]=vars
    return directives

def resolve_path(raw: Optional[str], conf: Path, cwd: Optional[str]) -> Optional[Path]:
    if not raw:return None
    if "$" in raw or "[" in raw or "]" in raw or "{" in raw or "}" in raw:return None
    p=Path(raw)
    if p.is_absolute():return p.resolve()
    if cwd:
        c=Path(cwd)
        if c.is_absolute():return (c/p).resolve()
        return (conf.parent/c/p).resolve()
    return (conf.parent/p).resolve()

def namdbin_natom(p: Path) -> Optional[int]:
    try:
        size=p.stat().st_size
        with p.open("rb") as f:b=f.read(4)
        if len(b)!=4:return None
        vals=[struct.unpack("<i",b)[0],struct.unpack(">i",b)[0]]
        plausible=[n for n in vals if n>0 and size in (4+n*3*8,4+n*3*4)]
        return plausible[0] if plausible else None
    except Exception:return None

def conf_prod_interval(p: Path):
    s=str(p).replace("\\","/");m=re.search(r"PROD_(d{3})_(d{3})ns",s,re.I)
    if not m:return None
    a,b=float(m.group(1)),float(m.group(2));return (a,b) if b>a else None

def load_expected_atoms(qcdir: Path):
    rows=read_tsv(qcdir/"PRODUCTION_CHAIN.tsv");dcd={r["relative_path"]:r for r in read_tsv(qcdir/"DCD_HEADERS.tsv")};by={}
    for r in rows:by.setdefault(r["target"],set()).add(int(dcd[r["relative_path"]]["natom"]))
    return {t:(next(iter(s)) if len(s)==1 else None) for t,s in by.items()}

def candidate_confs(root: Path,target: str):
    arr=[]
    for p in (root/target).rglob("*"):
        if p.is_file() and p.suffix.lower() in (".conf",".namd") and "PROD_" in str(p).upper():
            iv=conf_prod_interval(p)
            if iv:arr.append((p,iv))
    return sorted(arr,key=lambda x:(x[1][0],x[1][1],str(x[0])))

def choose_boundary_conf(confs,boundary: float,role: str):
    x=[(p,iv) for p,iv in confs if abs((iv[0] if role=="start" else iv[1])-boundary)<1e-9]
    if role=="end" and x:
        resume=[z for z in x if "RESUME" in str(z[0]).upper()]
        if resume:x=resume
    return x

def resolve_input_state(conf: Path,parsed):
    raw=parsed.get("bincoordinates") or parsed.get("coordinates");return raw,resolve_path(raw,conf,parsed.get("cwd"))
def resolve_output_state(conf: Path,parsed):
    raw=parsed.get("outputname");prefix=resolve_path(raw,conf,parsed.get("cwd"));return raw,(Path(str(prefix)+".coor") if prefix else None)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--root",default="/home/kaan/NAMD/systems");ap.add_argument("--qc",default="/home/kaan/NAMD/Publication_Analysis/00_QC");ap.add_argument("--out",default="/home/kaan/NAMD/Publication_Analysis/01_STATE_DISCOVERY");ap.add_argument("--references",default=None);args=ap.parse_args()
    root=Path(args.root).resolve();qc=Path(args.qc).resolve();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True);refs=Path(args.references).resolve() if args.references else Path(__file__).resolve().parents[1]/"references"/"REFERENCES_VERIFIED.tsv"
    qc_summary=json.loads((qc/"QC_SUMMARY.json").read_text()) if (qc/"QC_SUMMARY.json").exists() else {};expected_atoms=load_expected_atoms(qc) if (qc/"PRODUCTION_CHAIN.tsv").exists() and (qc/"DCD_HEADERS.tsv").exists() else {}
    config_rows=[];state_rows=[];resolution_rows=[]
    for target in ("ADK","MAPK"):
        confs=candidate_confs(root,target);parsed_map={}
        for p,iv in confs:
            q=parse_conf(p);parsed_map[p]=q
            config_rows.append({"target":target,"relative_path":str(p.relative_to(root)),"nominal_start_ns":iv[0],"nominal_end_ns":iv[1],"firsttimestep":q.get("firsttimestep",""),"timestep_fs":q.get("timestep",""),"dcdfreq":q.get("dcdfreq",""),"run_steps":",".join(map(str,q.get("run_steps",[]))),"bincoordinates_raw":q.get("bincoordinates",""),"coordinates_raw":q.get("coordinates",""),"outputname_raw":q.get("outputname",""),"cwd":q.get("cwd","")})
        for label,ns,role in [("T1",0.0,"start"),("T2",50.0,"start"),("T3",100.0,"end")]:
            choices=choose_boundary_conf(confs,ns,role);candidates=[];exp=expected_atoms.get(target)
            for p,iv in choices:
                q=parsed_map[p]
                if role=="start":raw,sp=resolve_input_state(p,q);kind="binCoordinates" if q.get("bincoordinates") else "coordinates"
                else:raw,sp=resolve_output_state(p,q);kind="outputName.coor"
                exists=bool(sp and sp.exists());nat=namdbin_natom(sp) if exists and sp.suffix.lower()==".coor" else None;candidates.append((p,iv,q,raw,sp,kind,exists,nat))
            good=[c for c in candidates if c[6] and (c[7] in (None,exp))];chosen=good[0] if len(good)==1 else None
            if chosen:
                p,iv,q,raw,sp,kind,exists,nat=chosen;status="RESOLVED";sha=sha256_file(sp);relstate=str(sp.relative_to(root)) if root in sp.parents else str(sp)
            else:
                p=choices[0][0] if len(choices)==1 else None;status="AMBIGUOUS_OR_MISSING";kind=nat=sha=relstate=""
            resolution_rows.append({"target":target,"label":label,"time_ns":ns,"status":status,"config":str(p.relative_to(root)) if p else "","source_kind":kind if chosen else "","state_path":relstate if chosen else "","state_exists":True if chosen else False,"state_natom":nat if chosen else "","expected_natom":exp or "","sha256":sha if chosen else "","candidate_count":len(candidates)})
            for c in candidates:
                cp,civ,cq,craw,csp,ckind,cexists,cnat=c
                state_rows.append({"target":target,"label":label,"time_ns":ns,"config":str(cp.relative_to(root)),"source_kind":ckind,"raw_value":craw or "","resolved_path":str(csp) if csp else "","exists":cexists,"natom":cnat or "","expected_natom":exp or "","atomcount_match":(cnat==exp) if cnat is not None and exp is not None else ""})
    p1=root.exists() and all((root/t).is_dir() for t in ("ADK","MAPK"));p2=qc_summary.get("final_status")=="PASS" and qc_summary.get("publication_eligible_for_next_stage") is True;p3=all(any(r["target"]==t for r in config_rows) for t in ("ADK","MAPK"));p4=all(sum(1 for r in resolution_rows if r["target"]==t)==3 for t in ("ADK","MAPK"));p5=all(r["status"]=="RESOLVED" for r in resolution_rows if r["label"] in ("T1","T2"));p6=all(r["status"]=="RESOLVED" for r in resolution_rows if r["label"]=="T3");p7=all(str(r["state_natom"])==str(r["expected_natom"]) for r in resolution_rows);p8=(all(float(r["time_ns"]) in (0.0,50.0,100.0) for r in resolution_rows) and [r["label"] for r in resolution_rows if r["target"]=="ADK"]==["T1","T2","T3"] and [r["label"] for r in resolution_rows if r["target"]=="MAPK"]==["T1","T2","T3"]);p9=True
    for t in ("ADK","MAPK"):
        rr=[r for r in resolution_rows if r["target"]==t]
        if len(set(r["state_path"] for r in rr))!=3 or not all(re.fullmatch(r"[0-9a-f]{64}",r["sha256"]) for r in rr):p9=False
    refrows=read_tsv(refs) if refs.exists() else [];p10=bool(refrows) and all(r.get("title") and r.get("doi") and r.get("pmid") and r.get("verified_status")=="VERIFIED" for r in refrows)
    passes=[("P1_DATA_ROOT_AND_TARGETS",p1),("P2_PRIOR_QC_10PASS",p2),("P3_PRODUCTION_CONFIG_DISCOVERY",p3),("P4_THREE_STATE_SLOTS_PER_TARGET",p4),("P5_T1_T2_INPUT_STATE_RESOLUTION",p5),("P6_T3_FINAL_STATE_RESOLUTION",p6),("P7_STATE_ATOMCOUNT_MATCH",p7),("P8_T1_T2_T3_TIME_IDENTITY",p8),("P9_STATE_HASH_AND_UNIQUENESS",p9),("P10_REFERENCE_IDENTITY_MANIFEST",p10)]
    write_tsv(out/"PRODUCTION_CONFIGS.tsv",config_rows,["target","relative_path","nominal_start_ns","nominal_end_ns","firsttimestep","timestep_fs","dcdfreq","run_steps","bincoordinates_raw","coordinates_raw","outputname_raw","cwd"]);write_tsv(out/"STATE_CANDIDATES.tsv",state_rows,["target","label","time_ns","config","source_kind","raw_value","resolved_path","exists","natom","expected_natom","atomcount_match"]);write_tsv(out/"T1_T2_T3_RESOLUTION.tsv",resolution_rows,["target","label","time_ns","status","config","source_kind","state_path","state_exists","state_natom","expected_natom","sha256","candidate_count"])
    summary={"version":VERSION,"passes":[{"id":k,"status":"PASS" if v else "FAIL"} for k,v in passes],"final_status":"PASS" if all(v for _,v in passes) else "FAIL","snapshot_extraction_eligible":bool(all(v for _,v in passes)),"note":"Exact boundary states are resolved from NAMD input/final coordinate states; DCD nearest frames are not relabeled as exact 0/50/100 ns states."};(out/"STATE_DISCOVERY_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    with (out/"10PASS_REPORT.txt").open("w") as f:
        f.write("NAMD CORDYCEPIN MAPK/ADK — 01_STATE_DISCOVERY 10-PASS REPORT\n"+"="*74+"\n")
        for k,v in passes:f.write(f"{k:<44} {'PASS' if v else 'FAIL'}\n")
        f.write("-"*74+f"\nFINAL STATUS: {summary['final_status']}\nSNAPSHOT EXTRACTION ELIGIBLE: {'YES' if summary['snapshot_extraction_eligible'] else 'NO'}\n")
    print(json.dumps(summary,indent=2));return 0 if all(v for _,v in passes) else 2
if __name__=="__main__":raise SystemExit(main())
