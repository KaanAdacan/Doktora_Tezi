#!/usr/bin/env python3
from pathlib import Path
import tempfile,json,csv,subprocess,sys
HERE=Path(__file__).resolve();SCRIPT=HERE.parents[1]/"code"/"topology_mode_audit.py";REFS=HERE.parents[1]/"references"/"REFERENCES_VERIFIED.tsv"
with tempfile.TemporaryDirectory() as td:
    td=Path(td);root=td/"systems";s1=td/"s1";s2=td/"s2";out=td/"out";s1.mkdir();s2.mkdir()
    for t in ("ADK","MAPK"):
        d=root/t/"active";d.mkdir(parents=True)
        atoms=[("PROA","1","ALA","N","NH1"),("PROA","1","ALA","CA","CT1"),("PROA","1","ALA","C","C"),("PROA","1","ALA","O","O"),("WAT","1","TIP3","OH2","OT"),("WAT","1","TIP3","H1","HT"),("WAT","1","TIP3","H2","HT"),("ION","1","SOD","SOD","SOD")]
        lines=["PSF","",f"{len(atoms):8d} !NATOM"]
        for i,a in enumerate(atoms,1):
            sg,ri,rn,an,at=a;lines.append(f"{i:8d} {sg} {ri} {rn} {an} {at} 0.0 1.0 0")
        (d/"system.psf").write_text("\n".join(lines)+"\n")
    (s1/"STATE_DISCOVERY_SUMMARY.json").write_text(json.dumps({"final_status":"PASS","snapshot_extraction_eligible":True}))
    chosen={t:{"psf":f"{t}/active/system.psf","pdb":f"{t}/active/system.pdb","config":"x"} for t in ("ADK","MAPK")}
    passes=[{"id":"P3_UNIQUE_SIMULATION_PSF_PDB_PAIR","status":"PASS"},{"id":"P4_PSF_PDB_ATOM_ORDER_SIGNATURE","status":"PASS"},{"id":"P5_UNIQUE_PLAUSIBLE_LIGAND_RESIDUE","status":"FAIL"},{"id":"P6_SIX_EXACT_SNAPSHOTS_WRITTEN","status":"PASS"},{"id":"P7_SNAPSHOT_ATOMCOUNT_MATCH","status":"PASS"},{"id":"P8_SNAPSHOT_TOPOLOGY_SIGNATURE","status":"PASS"},{"id":"P9_WITHIN_TARGET_SNAPSHOT_HASH_UNIQUENESS","status":"PASS"},{"id":"P10_REFERENCE_IDENTITY_MANIFEST","status":"PASS"}]
    (s2/"STAGE02_SUMMARY.json").write_text(json.dumps({"chosen_topology":chosen,"passes":passes,"final_status":"FAIL","rmsd_stage_eligible":False}))
    (s2/"LIGAND_CANDIDATES.tsv").write_text("target\tsegid\tresid\tresname\tatom_count\tatom_names\n")
    rows=[]
    for t in ("ADK","MAPK"):
      for lab,ns in [("T1",0),("T2",50),("T3",100)]:
        p=s2/f"{t}_{lab}.pdb";p.write_text("END\n");rows.append({"target":t,"label":lab,"time_ns":float(ns),"output_pdb":str(p),"topology_signature_match":"True","finite_coordinates":"True"})
    with (s2/"SNAPSHOT_MANIFEST.tsv").open("w",newline="") as f:
      w=csv.DictWriter(f,fieldnames=["target","label","time_ns","output_pdb","topology_signature_match","finite_coordinates"],delimiter="\t");w.writeheader();w.writerows(rows)
    cmd=[sys.executable,str(SCRIPT),"--root",str(root),"--stage01",str(s1),"--stage02",str(s2),"--out",str(out),"--references",str(REFS)]
    r1=subprocess.run(cmd,capture_output=True,text=True);summ=json.loads((out/"STAGE02A_SUMMARY.json").read_text());residue_rows=list(csv.DictReader((out/"RESIDUE_COMPOSITION.tsv").open(),delimiter="\t"));unexplained_rows=list(csv.DictReader((out/"UNEXPLAINED_RESIDUES.tsv").open(),delimiter="\t"))
    checks={"V1_run_exit_zero":r1.returncode==0,"V2_all_internal_10pass":all(x["status"]=="PASS" for x in summ["passes"]),"V3_both_targets_resolved":set(summ["target_modes"])=={"ADK","MAPK"},"V4_both_targets_no_small_molecule":all(v=="NO_NONSTANDARD_SMALL_MOLECULE_DETECTED" for v in summ["target_modes"].values()),"V5_unexplained_count_zero":summ["unexplained_count"]==0,"V6_unexplained_table_empty":len(unexplained_rows)==0,"V7_protein_rows_present":all(any(x["target"]==t and x["class"]=="PROTEIN" for x in residue_rows) for t in ("ADK","MAPK")),"V8_water_rows_recognized":all(any(x["target"]==t and x["class"]=="WATER" for x in residue_rows) for t in ("ADK","MAPK")),"V9_ion_rows_recognized":all(any(x["target"]==t and x["class"]=="ION" for x in residue_rows) for t in ("ADK","MAPK")),"V10_final_eligible":summ["apo_protein_metrics_eligible"] is True}
    rep=HERE.parent/"KNOWN_ANSWER_VERIFICATION.txt";rep.write_text("\n".join(["02A TOPOLOGY MODE AUDIT VERIFICATION","="*70]+[f"{k:<48} {'PASS' if v else 'FAIL'}" for k,v in checks.items()]+["-"*70,"FINAL: 10/10 PASS" if all(checks.values()) else "FINAL: FAIL"])+"\n");print(rep.read_text())
    if not all(checks.values()):print(r1.stdout,r1.stderr);raise SystemExit(2)
