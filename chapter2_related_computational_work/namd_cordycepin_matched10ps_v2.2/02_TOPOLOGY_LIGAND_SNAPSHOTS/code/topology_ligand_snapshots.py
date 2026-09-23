#!/usr/bin/env python3
"""
02_TOPOLOGY_LIGAND_SNAPSHOTS

Purpose:
1) Resolve a simulation PSF + template PDB for each target.
2) Verify atom ordering/signature at topology level.
3) Discover the simulated non-protein ligand without assuming its residue name.
4) Convert exact T1/T2/T3 NAMD binary coordinate states to PDB by replacing
   coordinates in the verified template PDB atom order.
5) Fail closed on ambiguity.

This module does not calculate RMSD or any other scientific observable.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, re, struct
from pathlib import Path
from collections import Counter, defaultdict

VERSION="1.0.0"

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HSD","HSE","HSP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","CYX","ASH","GLH","LYN"}
WATER = {"TIP3","TIP3P","HOH","WAT","SOL"}
IONS = {"SOD","CLA","NA","CL","K","POT","CAL","MG","ZN","ZN2","CA","MG2"}
CAPS = {"ACE","NME","CT3","NH2"}

def sha256_file(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def read_tsv(p):
    with Path(p).open(encoding="utf-8") as f:return list(csv.DictReader(f,delimiter="	"))
def write_tsv(p,rows,fields):
    with Path(p).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="	",extrasaction="ignore");w.writeheader()
        for r in rows:w.writerow({k:r.get(k,"") for k in fields})
def strip_comment(s):return s.split("#",1)[0].strip()
VAR1=re.compile(r"${([A-Za-z_][A-Za-z0-9_]*)}");VAR2=re.compile(r"$([A-Za-z_][A-Za-z0-9_]*)")
def subst(s,vars):
    old=None
    for _ in range(20):
        if s==old:break
        old=s;s=VAR1.sub(lambda m:vars.get(m.group(1),m.group(0)),s);s=VAR2.sub(lambda m:vars.get(m.group(1),m.group(0)),s)
    return s.strip().strip('"').strip("'")
def parse_conf(p):
    vars={};d={}
    try:lines=Path(p).read_text(errors="replace").splitlines()
    except:return {}
    for raw in lines:
        line=strip_comment(raw);m=re.match(r"(?i)^sets+([A-Za-z_][A-Za-z0-9_]*)s+(.+)$",line)
        if m:vars[m.group(1)]=subst(m.group(2),vars)
    for raw in lines:
        line=strip_comment(raw);m=re.match(r"^(S+)s+(.+?)s*$",line)
        if not m:continue
        k=m.group(1).lower();v=subst(m.group(2),vars)
        if k in ("structure","coordinates","cwd","bincoordinates","outputname"):d[k]=v
    d["vars"]=vars;return d
def resolve(raw,conf,cwd=None):
    if not raw or any(x in raw for x in ("$","[","]","{","}")):return None
    p=Path(raw)
    if p.is_absolute():return p.resolve()
    if cwd:
        c=Path(cwd)
        if c.is_absolute():return (c/p).resolve()
        return (Path(conf).parent/c/p).resolve()
    return (Path(conf).parent/p).resolve()
def parse_psf_atoms(p):
    lines=Path(p).read_text(errors="replace").splitlines();start=None;n=None
    for i,line in enumerate(lines):
        if "!NATOM" in line:n=int(line.split()[0]);start=i+1;break
    if n is None:raise ValueError("PSF !NATOM not found")
    atoms=[]
    for line in lines[start:start+n]:
        sp=line.split()
        if len(sp)<6:raise ValueError("Malformed PSF atom line")
        atoms.append({"index":int(sp[0]),"segid":sp[1],"resid":sp[2],"resname":sp[3],"atomname":sp[4],"atomtype":sp[5]})
    if len(atoms)!=n:raise ValueError("PSF atom count mismatch")
    return atoms
def parse_pdb_atoms(p):
    rows=[]
    for line in Path(p).read_text(errors="replace").splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            rows.append({"line":line,"atomname":line[12:16].strip(),"resname":line[17:20].strip(),"chain":line[21:22].strip(),"resid":line[22:26].strip(),"segid":line[72:76].strip() if len(line)>=76 else "","element":line[76:78].strip() if len(line)>=78 else ""})
    return rows
def canonical_pdb_resname(psf_resname):return str(psf_resname)[:3]
def signature_match(psf,pdb):
    if len(psf)!=len(pdb):return False,0,[]
    mism=[];nmatch=0
    for i,(a,b) in enumerate(zip(psf,pdb),start=1):
        ok=a["resid"]==b["resid"] and canonical_pdb_resname(a["resname"])==b["resname"] and a["atomname"]==b["atomname"]
        if ok:nmatch+=1
        elif len(mism)<20:mism.append((i,a["resid"],a["resname"],a["atomname"],b["resid"],b["resname"],b["atomname"]))
    return nmatch==len(psf),nmatch,mism
def common_prefix_depth(a,b):
    n=0
    for x,y in zip(a.parts,b.parts):
        if x!=y:break
        n+=1
    return n
def select_active_topology_pair(exact_pairs,root,active_config_relpaths):
    by_pair={}
    for conf,psf,pdb in exact_pairs:by_pair.setdefault((str(psf.resolve()),str(pdb.resolve())),[]).append((conf,psf,pdb))
    active=[Path(x) for x in active_config_relpaths if x];scored=[]
    for key,items in by_pair.items():
        best_item=None;best_score=-1
        for item in items:
            conf_rel=item[0].relative_to(root);score=max((common_prefix_depth(conf_rel,a) for a in active),default=0)
            if score>best_score:best_score=score;best_item=item
        scored.append((best_score,best_item,key))
    if not scored:return None,[]
    mx=max(x[0] for x in scored);top=[x for x in scored if x[0]==mx]
    return (top[0][1] if len(top)==1 else None),scored
def read_namdbin(p,expected_n):
    p=Path(p)
    with p.open("rb") as f:b=f.read()
    candidates=[]
    for endian in ("<",">"):
        n=struct.unpack(endian+"i",b[:4])[0]
        if n==expected_n and len(b)==4+n*3*8:
            vals=struct.unpack_from(endian+f"{n*3}d",b,4);coords=[(vals[i],vals[i+1],vals[i+2]) for i in range(0,len(vals),3)];candidates.append((endian,coords))
    if len(candidates)!=1:raise ValueError(f"NAMDBIN endianness/size resolution failed for {p}; candidates={len(candidates)} size={p.stat().st_size}")
    return candidates[0]
def write_pdb_from_template(template,coords,out):
    lines=Path(template).read_text(errors="replace").splitlines();j=0;result=[]
    for line in lines:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            x,y,z=coords[j];j+=1;line=line.ljust(80) if len(line)<80 else line;line=line[:30]+f"{x:8.3f}{y:8.3f}{z:8.3f}"+line[54:]
        result.append(line)
    if j!=len(coords):raise ValueError("Template coordinate replacement atom count mismatch")
    Path(out).write_text("\n".join(result)+"\n",encoding="utf-8")
def find_topology_pairs(root,target,expected_n):
    confs=[p for p in (root/target).rglob("*") if p.is_file() and p.suffix.lower() in (".conf",".namd")];explicit=[]
    for c in confs:
        q=parse_conf(c);psf=resolve(q.get("structure"),c,q.get("cwd"));pdb=resolve(q.get("coordinates"),c,q.get("cwd"))
        if psf and pdb and psf.exists() and pdb.exists() and psf.suffix.lower()==".psf" and pdb.suffix.lower()==".pdb":
            try:pa=parse_psf_atoms(psf);pb=parse_pdb_atoms(pdb);explicit.append((c,psf,pdb,len(pa),len(pb)))
            except:pass
    seen=set();uniq=[]
    for x in explicit:
        key=(str(x[1]),str(x[2]))
        if key not in seen:seen.add(key);uniq.append(x)
    return uniq,[x for x in uniq if x[3]==expected_n and x[4]==expected_n]
def residue_instances(psf):
    d=defaultdict(list)
    for a in psf:d[(a["segid"],a["resid"],a["resname"])].append(a)
    return d
def ligand_candidates(psf):
    rows=[]
    for (segid,resid,resname),atoms in residue_instances(psf).items():
        cls="OTHER"
        if resname in AA:cls="PROTEIN"
        elif resname in WATER:cls="WATER"
        elif resname in IONS:cls="ION"
        elif resname in CAPS:cls="CAP"
        if cls=="OTHER":rows.append({"segid":segid,"resid":resid,"resname":resname,"atom_count":len(atoms),"atom_names":",".join(a["atomname"] for a in atoms)})
    return rows
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--root",default="/home/kaan/NAMD/systems");ap.add_argument("--stage01",default="/home/kaan/NAMD/Publication_Analysis/01_STATE_DISCOVERY");ap.add_argument("--out",default="/home/kaan/NAMD/Publication_Analysis/02_TOPOLOGY_LIGAND_SNAPSHOTS");ap.add_argument("--references",default=None);args=ap.parse_args()
    root=Path(args.root).resolve();s1=Path(args.stage01).resolve();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True);refs=Path(args.references).resolve() if args.references else Path(__file__).resolve().parents[1]/"references"/"REFERENCES_VERIFIED.tsv"
    s1sum=json.loads((s1/"STATE_DISCOVERY_SUMMARY.json").read_text()) if (s1/"STATE_DISCOVERY_SUMMARY.json").exists() else {};states=read_tsv(s1/"T1_T2_T3_RESOLUTION.tsv") if (s1/"T1_T2_T3_RESOLUTION.tsv").exists() else [];expected={}
    for r in states:expected.setdefault(r["target"],set()).add(int(r["expected_natom"]))
    expected={t:(next(iter(v)) if len(v)==1 else None) for t,v in expected.items()};pair_rows=[];chosen={};ligand_rows=[];snapshot_rows=[];signature_fail=False;pair_ambig=False
    for target in ("ADK","MAPK"):
        allpairs,matching=find_topology_pairs(root,target,expected.get(target));exact=[]
        for conf,psf,pdb,np,nq in matching:
            try:pa=parse_psf_atoms(psf);pb=parse_pdb_atoms(pdb);ok,_,mism=signature_match(pa,pb)
            except Exception:ok=False;mism=[("parse_error",)]
            pair_rows.append({"target":target,"config":str(conf.relative_to(root)),"psf":str(psf.relative_to(root)) if root in psf.parents else str(psf),"pdb":str(pdb.relative_to(root)) if root in pdb.parents else str(pdb),"psf_natom":np,"pdb_natom":nq,"signature_match":ok,"mismatch_preview":"; ".join(map(str,mism[:5])),"lineage_score":"","selected_active_pair":False})
            if ok:exact.append((conf,psf,pdb))
        active_configs=[r["config"] for r in states if r["target"]==target and r.get("config")];selected,scored=select_active_topology_pair(exact,root,active_configs);score_by_pair={(str(item[1].resolve()),str(item[2].resolve())):score for score,item,key in scored}
        for rr in pair_rows:
            if rr["target"]==target:
                psf_abs=(root/rr["psf"]).resolve() if not Path(rr["psf"]).is_absolute() else Path(rr["psf"]).resolve();pdb_abs=(root/rr["pdb"]).resolve() if not Path(rr["pdb"]).is_absolute() else Path(rr["pdb"]).resolve();rr["lineage_score"]=score_by_pair.get((str(psf_abs),str(pdb_abs)),"");rr["selected_active_pair"]=bool(selected and psf_abs==selected[1].resolve() and pdb_abs==selected[2].resolve())
        if selected is None:pair_ambig=True;continue
        chosen[target]=selected;pa=parse_psf_atoms(selected[1]);cands=ligand_candidates(pa)
        for c in cands:ligand_rows.append({"target":target,**c})
        template=exact[0][2];snaps=out/"snapshots"/target;snaps.mkdir(parents=True,exist_ok=True)
        for r in [x for x in states if x["target"]==target]:
            state=(root/r["state_path"]).resolve() if not Path(r["state_path"]).is_absolute() else Path(r["state_path"]);endian,coords=read_namdbin(state,expected[target]);op=snaps/f"{target}_{r['label']}_{float(r['time_ns']):.0f}ns.pdb";write_pdb_from_template(template,coords,op);op_atoms=parse_pdb_atoms(op);ok2,_,mism2=signature_match(pa,op_atoms);finite=all(math.isfinite(v) for xyz in coords for v in xyz);snapshot_rows.append({"target":target,"label":r["label"],"time_ns":r["time_ns"],"state_path":r["state_path"],"template_pdb":str(template.relative_to(root)),"output_pdb":str(op),"natom":len(coords),"endianness":"little" if endian=="<" else "big","finite_coordinates":finite,"topology_signature_match":ok2,"sha256":sha256_file(op)})
            if not ok2:signature_fail=True
    refrows=read_tsv(refs) if refs.exists() else [];p1=root.exists() and all((root/t).is_dir() for t in ("ADK","MAPK"));p2=s1sum.get("final_status")=="PASS" and s1sum.get("snapshot_extraction_eligible") is True;p3=all(t in chosen for t in ("ADK","MAPK")) and not pair_ambig;p4=all(r["signature_match"] is True for r in pair_rows if r["target"] in chosen);p5=all(len([r for r in ligand_rows if r["target"]==t and 5<=int(r["atom_count"])<=100])==1 for t in ("ADK","MAPK"));p6=len(snapshot_rows)==6 and all(r["finite_coordinates"] for r in snapshot_rows);p7=all(int(r["natom"])==expected[r["target"]] for r in snapshot_rows);p8=all(r["topology_signature_match"] for r in snapshot_rows) and not signature_fail;p9=all(re.fullmatch(r"[0-9a-f]{64}",r["sha256"]) for r in snapshot_rows)
    for t in ("ADK","MAPK"):
        rr=[r for r in snapshot_rows if r["target"]==t]
        if len(rr)!=3 or len({r["sha256"] for r in rr})!=3:p9=False
    p10=bool(refrows) and all(r.get("title") and r.get("doi") and r.get("pmid") and r.get("verified_status")=="VERIFIED" for r in refrows);passes=[("P1_DATA_ROOT_AND_TARGETS",p1),("P2_STAGE01_10PASS",p2),("P3_UNIQUE_SIMULATION_PSF_PDB_PAIR",p3),("P4_PSF_PDB_ATOM_ORDER_SIGNATURE",p4),("P5_UNIQUE_PLAUSIBLE_LIGAND_RESIDUE",p5),("P6_SIX_EXACT_SNAPSHOTS_WRITTEN",p6),("P7_SNAPSHOT_ATOMCOUNT_MATCH",p7),("P8_SNAPSHOT_TOPOLOGY_SIGNATURE",p8),("P9_WITHIN_TARGET_SNAPSHOT_HASH_UNIQUENESS",p9),("P10_REFERENCE_IDENTITY_MANIFEST",p10)]
    write_tsv(out/"TOPOLOGY_PAIR_AUDIT.tsv",pair_rows,["target","config","psf","pdb","psf_natom","pdb_natom","signature_match","lineage_score","selected_active_pair","mismatch_preview"]);write_tsv(out/"LIGAND_CANDIDATES.tsv",ligand_rows,["target","segid","resid","resname","atom_count","atom_names"]);write_tsv(out/"SNAPSHOT_MANIFEST.tsv",snapshot_rows,["target","label","time_ns","state_path","template_pdb","output_pdb","natom","endianness","finite_coordinates","topology_signature_match","sha256"])
    summary={"version":VERSION,"chosen_topology":{t:{"config":str(chosen[t][0].relative_to(root)),"psf":str(chosen[t][1].relative_to(root)),"pdb":str(chosen[t][2].relative_to(root))} for t in chosen},"ligand_candidates":{t:[r for r in ligand_rows if r["target"]==t] for t in ("ADK","MAPK")},"passes":[{"id":k,"status":"PASS" if v else "FAIL"} for k,v in passes],"final_status":"PASS" if all(v for _,v in passes) else "FAIL","rmsd_stage_eligible":bool(all(v for _,v in passes))};(out/"STAGE02_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    with (out/"10PASS_REPORT.txt").open("w") as f:
        f.write("NAMD CORDYCEPIN MAPK/ADK — 02_TOPOLOGY_LIGAND_SNAPSHOTS 10-PASS REPORT\n"+"="*82+"\n")
        for k,v in passes:f.write(f"{k:<50} {'PASS' if v else 'FAIL'}\n")
        f.write("-"*82+f"\nFINAL STATUS: {summary['final_status']}\nRMSD STAGE ELIGIBLE: {'YES' if summary['rmsd_stage_eligible'] else 'NO'}\n")
    print(json.dumps(summary,indent=2));return 0 if all(v for _,v in passes) else 2
if __name__=="__main__":raise SystemExit(main())
