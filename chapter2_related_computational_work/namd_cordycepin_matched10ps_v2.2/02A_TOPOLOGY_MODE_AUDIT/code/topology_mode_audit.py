#!/usr/bin/env python3
"""
02A_TOPOLOGY_MODE_AUDIT

Determines whether the verified NAMD topology contains any simulated
non-protein small-molecule residue. It is designed for the observed case in
which Stage 02 produced valid topology/snapshots but no ligand candidate.

This stage does NOT infer that docking ligand was simulated merely from the
study design. It inspects the verified PSF composition and fails closed on
unexplained nonstandard residues.
"""
from __future__ import annotations
import argparse, csv, json, re, hashlib
from pathlib import Path
from collections import defaultdict, Counter

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HSD","HSE","HSP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","CYX","ASH","GLH","LYN"}
WATER = {"TIP3","TIP3P","TIP","HOH","WAT","SOL"}
IONS = {"SOD","CLA","NA","CL","K","POT","CAL","MG","ZN","ZN2","CA","MG2"}
CAPS = {"ACE","NME","CT3","NH2"}
NUCLEIC = {"A","C","G","U","T","DA","DC","DG","DT","ADE","CYT","GUA","THY","URA"}

def read_tsv(p):
    with Path(p).open(encoding="utf-8") as f:
        return list(csv.DictReader(f,delimiter="	"))

def write_tsv(p,rows,fields):
    with Path(p).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="	",extrasaction="ignore")
        w.writeheader()
        for r in rows: w.writerow({k:r.get(k,"") for k in fields})

def parse_psf(p):
    lines=Path(p).read_text(errors="replace").splitlines()
    start=n=None
    for i,line in enumerate(lines):
        if "!NATOM" in line:
            n=int(line.split()[0]); start=i+1; break
    if n is None: raise ValueError("!NATOM missing")
    atoms=[]
    for line in lines[start:start+n]:
        sp=line.split()
        if len(sp)<6: raise ValueError("Malformed PSF atom line")
        atoms.append({"index":int(sp[0]),"segid":sp[1],"resid":sp[2],"resname":sp[3],"atomname":sp[4],"atomtype":sp[5]})
    if len(atoms)!=n: raise ValueError("PSF atom count mismatch")
    return atoms

def classify(resname):
    r=resname.upper()
    if r in AA: return "PROTEIN"
    if r in WATER: return "WATER"
    if r in IONS: return "ION"
    if r in CAPS: return "CAP"
    if r in NUCLEIC: return "NUCLEIC_OR_NUCLEOSIDE"
    return "OTHER"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",default="/home/kaan/NAMD/systems")
    ap.add_argument("--stage01",default="/home/kaan/NAMD/Publication_Analysis/01_STATE_DISCOVERY")
    ap.add_argument("--stage02",default="/home/kaan/NAMD/Publication_Analysis/02_TOPOLOGY_LIGAND_SNAPSHOTS")
    ap.add_argument("--out",default="/home/kaan/NAMD/Publication_Analysis/02A_TOPOLOGY_MODE_AUDIT")
    ap.add_argument("--references",default=None)
    args=ap.parse_args()
    root=Path(args.root).resolve(); s1=Path(args.stage01).resolve(); s2=Path(args.stage02).resolve(); out=Path(args.out).resolve()
    out.mkdir(parents=True,exist_ok=True)
    refs=Path(args.references).resolve() if args.references else Path(__file__).resolve().parents[1]/"references"/"REFERENCES_VERIFIED.tsv"

    s1sum=json.loads((s1/"STATE_DISCOVERY_SUMMARY.json").read_text()) if (s1/"STATE_DISCOVERY_SUMMARY.json").exists() else {}
    s2sum=json.loads((s2/"STAGE02_SUMMARY.json").read_text()) if (s2/"STAGE02_SUMMARY.json").exists() else {}
    s2passes={x["id"]:x["status"] for x in s2sum.get("passes",[])}

    comp=[]; unexplained=[]; target_modes={}
    parse_ok=True
    for target in ("ADK","MAPK"):
        chosen=s2sum.get("chosen_topology",{}).get(target)
        if not chosen:
            parse_ok=False; target_modes[target]="UNRESOLVED"; continue
        psf=(root/chosen["psf"]).resolve()
        try:
            atoms=parse_psf(psf)
        except Exception:
            parse_ok=False; target_modes[target]="UNRESOLVED"; continue
        inst=defaultdict(list)
        for a in atoms:
            inst[(a["segid"],a["resid"],a["resname"])].append(a)
        for (segid,resid,resname),aa in sorted(inst.items()):
            cls=classify(resname)
            row={"target":target,"segid":segid,"resid":resid,"resname":resname,"class":cls,"atom_count":len(aa)}
            comp.append(row)
            if cls in ("OTHER","NUCLEIC_OR_NUCLEOSIDE") and 5 <= len(aa) <= 200:
                unexplained.append(row)
        target_modes[target]="NO_NONSTANDARD_SMALL_MOLECULE_DETECTED" if not any(r["target"]==target for r in unexplained) else "AMBIGUOUS_NONSTANDARD_RESIDUE_PRESENT"

    snaps=read_tsv(s2/"SNAPSHOT_MANIFEST.tsv") if (s2/"SNAPSHOT_MANIFEST.tsv").exists() else []
    ligs=read_tsv(s2/"LIGAND_CANDIDATES.tsv") if (s2/"LIGAND_CANDIDATES.tsv").exists() else []
    refrows=read_tsv(refs) if refs.exists() else []

    required_pass_ids=[
      "P3_UNIQUE_SIMULATION_PSF_PDB_PAIR",
      "P4_PSF_PDB_ATOM_ORDER_SIGNATURE",
      "P6_SIX_EXACT_SNAPSHOTS_WRITTEN",
      "P7_SNAPSHOT_ATOMCOUNT_MATCH",
      "P8_SNAPSHOT_TOPOLOGY_SIGNATURE",
      "P9_WITHIN_TARGET_SNAPSHOT_HASH_UNIQUENESS",
      "P10_REFERENCE_IDENTITY_MANIFEST",
    ]

    p1=root.exists() and all((root/t).is_dir() for t in ("ADK","MAPK"))
    p2=s1sum.get("final_status")=="PASS" and s1sum.get("snapshot_extraction_eligible") is True
    p3=all(s2passes.get(k)=="PASS" for k in required_pass_ids)
    p4=s2passes.get("P5_UNIQUE_PLAUSIBLE_LIGAND_RESIDUE")=="FAIL" and len(ligs)==0
    p5=parse_ok and all(t in s2sum.get("chosen_topology",{}) for t in ("ADK","MAPK"))
    p6=len(unexplained)==0 and all(target_modes.get(t)=="NO_NONSTANDARD_SMALL_MOLECULE_DETECTED" for t in ("ADK","MAPK"))
    p7=len(snaps)==6 and all(Path(r["output_pdb"]).exists() for r in snaps)
    p8=all(r.get("topology_signature_match")=="True" and r.get("finite_coordinates")=="True" for r in snaps)
    p9={r["label"] for r in snaps}=={"T1","T2","T3"} and sorted({float(r["time_ns"]) for r in snaps})==[0.0,50.0,100.0]
    p10=bool(refrows) and all(r.get("title") and r.get("doi") and r.get("pmid") and r.get("verified_status")=="VERIFIED" for r in refrows)

    passes=[
      ("P1_DATA_ROOT_AND_TARGETS",p1),
      ("P2_STAGE01_EXACT_STATES_PASS",p2),
      ("P3_STAGE02_TOPOLOGY_SNAPSHOT_COMPONENTS_PASS",p3),
      ("P4_STAGE02_FAILURE_IS_LIGAND_ONLY",p4),
      ("P5_VERIFIED_PSF_COMPOSITION_PARSED",p5),
      ("P6_NO_NONSTANDARD_SMALL_MOLECULE_DETECTED",p6),
      ("P7_SIX_EXACT_SNAPSHOT_FILES_PRESENT",p7),
      ("P8_SNAPSHOT_SIGNATURE_AND_COORDINATES_VALID",p8),
      ("P9_T1_T2_T3_IDENTITY_VALID",p9),
      ("P10_REFERENCE_IDENTITY_MANIFEST",p10),
    ]

    write_tsv(out/"RESIDUE_COMPOSITION.tsv",comp,["target","segid","resid","resname","class","atom_count"])
    write_tsv(out/"UNEXPLAINED_RESIDUES.tsv",unexplained,["target","segid","resid","resname","class","atom_count"])
    summary={
      "target_modes":target_modes,
      "unexplained_count":len(unexplained),
      "interpretation":"No simulated nonstandard small-molecule residue was detected in the verified PSF topology. The 100-ns trajectories should therefore be treated as receptor conformational sampling, not ligand-bound MD, unless independent topology evidence demonstrates otherwise.",
      "passes":[{"id":k,"status":"PASS" if v else "FAIL"} for k,v in passes],
      "final_status":"PASS" if all(v for _,v in passes) else "FAIL",
      "apo_protein_metrics_eligible":bool(all(v for _,v in passes))
    }
    (out/"STAGE02A_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    with (out/"10PASS_REPORT.txt").open("w") as f:
        f.write("NAMD CORDYCEPIN MAPK/ADK — 02A TOPOLOGY MODE AUDIT 10-PASS REPORT\n")
        f.write("="*86+"\n")
        for k,v in passes: f.write(f"{k:<55} {'PASS' if v else 'FAIL'}\n")
        f.write("-"*86+"\n")
        f.write(f"FINAL STATUS: {summary['final_status']}\n")
        f.write(f"APO PROTEIN METRICS ELIGIBLE: {'YES' if summary['apo_protein_metrics_eligible'] else 'NO'}\n")
    print(json.dumps(summary,indent=2))
    return 0 if all(v for _,v in passes) else 2

if __name__=="__main__":
    raise SystemExit(main())
