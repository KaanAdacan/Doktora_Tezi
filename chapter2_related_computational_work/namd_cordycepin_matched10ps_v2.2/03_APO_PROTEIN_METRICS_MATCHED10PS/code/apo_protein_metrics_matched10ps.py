#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,re,sys
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd as mda_rmsd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from sampling_policy import TARGET_SAMPLING, matched_indices


def read_tsv(p):
    with Path(p).open(encoding="utf-8") as f:return list(csv.DictReader(f,delimiter="\t"))

def write_tsv(p,rows,fields):
    with Path(p).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore");w.writeheader()
        for r in rows:w.writerow({k:r.get(k,"") for k in fields})

def _production_identity(rel):
    p=Path(rel);local=(p.parent.name+"/"+p.name).upper()
    return ("PROD_" in local) and not bool(re.search(r"(^|[/_.-])(NVT|NPT|BENCH)([/_.-]|$)",local))

def build_axis(chain,target,root):
    rr=sorted([r for r in chain if r["target"]==target],key=lambda r:float(r["actual_start_ns"]))
    paths=[];times=[];prev=0.0
    for i,r in enumerate(rr):
        rel=r["relative_path"]
        if not _production_identity(rel):raise ValueError(f"Non-production DCD: {rel}")
        a=float(r["actual_start_ns"]);b=float(r["actual_end_ns"]);n=int(r["nset"])
        if i==0 and abs(a)>1e-9:raise ValueError("Production does not start at 0")
        if i and abs(a-prev)>1e-9:raise ValueError("Production chain discontinuity")
        step=(b-a)/n
        times.extend(a+(j+1)*step for j in range(n))
        p=(root/rel).resolve()
        if not p.exists():raise FileNotFoundError(p)
        paths.append(str(p));prev=b
    if abs(prev-100.0)>1e-8:raise ValueError("Production does not end at 100 ns")
    return paths,np.asarray(times,dtype=np.float64),rr

def plot(x,y,xlab,ylab,title,png,pdf):
    import matplotlib;matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.2,4.6));ax.plot(x,y,linewidth=.8)
    ax.set_xlabel(xlab);ax.set_ylabel(ylab);ax.set_title(title);ax.grid(True,alpha=.15)
    fig.tight_layout();fig.savefig(png,dpi=600,bbox_inches="tight");fig.savefig(pdf,bbox_inches="tight");plt.close(fig)

def target_run(target,root,qc,s2,out):
    chain=read_tsv(qc/"PRODUCTION_CHAIN.tsv");snaps=read_tsv(s2/"SNAPSHOT_MANIFEST.tsv")
    s2sum=json.loads((s2/"STAGE02_SUMMARY.json").read_text());chosen=s2sum["chosen_topology"][target]
    top=(root/chosen["psf"]).resolve();
    if not top.exists():raise FileNotFoundError(top)
    t1=[r for r in snaps if r["target"]==target and r["label"]=="T1"]
    if len(t1)!=1:raise ValueError("T1 resolution error")
    refp=Path(t1[0]["output_pdb"]).resolve()
    dcds,native_times,chain_rows=build_axis(chain,target,root)
    sample_idx,times=matched_indices(native_times,target)
    sample_set=set(int(i) for i in sample_idx)

    u=mda.Universe(str(top),dcds,format="DCD");ref=mda.Universe(str(top),str(refp))
    ca=u.select_atoms("protein and name CA");bb=u.select_atoms("protein and backbone");prot=u.select_atoms("protein")
    rca=ref.select_atoms("protein and name CA");rbb=ref.select_atoms("protein and backbone")
    if len(u.trajectory)!=len(native_times):raise ValueError("Native frame/time mismatch")
    if min(len(ca),len(bb),len(prot))<=0:raise ValueError("Empty protein selection")
    pm=np.asarray(prot.masses,dtype=float)
    if len(pm)!=len(prot) or not np.all(np.isfinite(pm)) or np.any(pm<=0):raise ValueError("Invalid masses")
    refca=rca.positions.copy();refbb=rbb.positions.copy()
    mean=np.zeros((len(ca),3));M2=np.zeros((len(ca),3));n=0;rows=[];out_i=0
    for raw_i,ts in enumerate(u.trajectory):
        if raw_i not in sample_set:continue
        rg=float(prot.radius_of_gyration())
        car=float(mda_rmsd(ca.positions,refca,center=True,superposition=True))
        align.alignto(u,ref,select="protein and backbone",weights=None,match_atoms=False)
        bbr=float(mda_rmsd(bb.positions,refbb,center=False,superposition=False))
        n+=1;x=ca.positions.astype(float);delta=x-mean;mean+=delta/n;M2+=delta*(x-mean)
        rows.append({"target":target,"analysis_frame_index":out_i+1,"raw_frame_index":raw_i+1,
                     "time_ns":f"{times[out_i]:.8f}","native_stride":TARGET_SAMPLING[target]["stride"],
                     "ca_rmsd_A":f"{car:.8f}","backbone_rmsd_A":f"{bbr:.8f}","protein_rg_A":f"{rg:.8f}"})
        out_i+=1
    if out_i!=10000:raise ValueError(f"{target}: analyzed {out_i}, expected 10000")
    rmsf=np.sqrt(np.sum(M2/n,axis=1))
    rrows=[{"target":target,"segid":a.segid,"resid":a.resid,"resname":a.resname,"atomname":a.name,"ca_rmsf_A":f"{float(v):.8f}"} for a,v in zip(ca,rmsf)]
    td=out/target;(td/"figures").mkdir(parents=True,exist_ok=True)
    write_tsv(td/"PROTEIN_STRUCTURAL_METRICS_MATCHED10PS.tsv",rows,
              ["target","analysis_frame_index","raw_frame_index","time_ns","native_stride","ca_rmsd_A","backbone_rmsd_A","protein_rg_A"])
    write_tsv(td/"CA_RMSF_MATCHED10PS.tsv",rrows,["target","segid","resid","resname","atomname","ca_rmsf_A"])
    x=np.array([float(r["time_ns"]) for r in rows])
    for fld,ylab,stem in [("ca_rmsd_A","Cα RMSD (Å)","CA_RMSD"),("backbone_rmsd_A","Backbone RMSD (Å)","BACKBONE_RMSD"),("protein_rg_A","Protein radius of gyration (Å)","RG")]:
        y=np.array([float(r[fld]) for r in rows]);plot(x,y,"Production time (ns)",ylab,f"{target} — {stem.replace('_',' ')}, matched 10-ps sampling",td/"figures"/f"{stem}.png",td/"figures"/f"{stem}.pdf")
    plot(np.arange(1,len(rrows)+1),rmsf,"Cα residue index","Cα RMSF (Å)",f"{target} — Cα RMSF, matched 10-ps sampling",td/"figures"/"CA_RMSF.png",td/"figures"/"CA_RMSF.pdf")
    blocks=[]
    for name,rr in [("0-50",[r for r in rows if float(r["time_ns"])<=50+1e-9]),("50-100",[r for r in rows if float(r["time_ns"])>50+1e-9])]:
        rec={"target":target,"block_ns":name,"n_frames":len(rr)}
        for fld in ("ca_rmsd_A","backbone_rmsd_A","protein_rg_A"):
            a=np.array([float(r[fld]) for r in rr]);rec[fld+"_mean"]=f"{a.mean():.8f}";rec[fld+"_sd"]=f"{a.std(ddof=1):.8f}";rec[fld+"_median"]=f"{np.median(a):.8f}"
        blocks.append(rec)
    fields=["target","block_ns","n_frames"]
    for fld in ("ca_rmsd_A","backbone_rmsd_A","protein_rg_A"):fields += [fld+"_mean",fld+"_sd",fld+"_median"]
    write_tsv(td/"BLOCK_0_50_vs_50_100_MATCHED10PS.tsv",blocks,fields)
    manifest={"target":target,"sampling_policy":"DETERMINISTIC_MATCHED_10PS","native_interval_ps":TARGET_SAMPLING[target]["native_interval_ns"]*1000,
              "native_stride":TARGET_SAMPLING[target]["stride"],"analysis_interval_ps":10.0,"native_frames":len(native_times),"analysis_frames":len(rows),
              "first_analysis_time_ns":float(times[0]),"last_analysis_time_ns":float(times[-1]),"raw_trajectory_modified":False,
              "excluded_operations":["interpolation","smoothing","filtering","imputation","outlier removal"],"dcds":dcds}
    (td/"ANALYSIS_MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    return rows,rrows,blocks,dcds,native_times,times,manifest

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--root",default="/home/kaan/NAMD/systems");ap.add_argument("--qc",default="/home/kaan/NAMD/Publication_Analysis/00_QC")
    ap.add_argument("--stage02",default="/home/kaan/NAMD/Publication_Analysis/02_TOPOLOGY_LIGAND_SNAPSHOTS");ap.add_argument("--stage02a",default="/home/kaan/NAMD/Publication_Analysis/02A_TOPOLOGY_MODE_AUDIT")
    ap.add_argument("--out",default="/home/kaan/NAMD/Publication_Analysis/03_APO_PROTEIN_METRICS_MATCHED10PS");ap.add_argument("--references",default=None);args=ap.parse_args()
    root=Path(args.root).resolve();qc=Path(args.qc).resolve();s2=Path(args.stage02).resolve();s2a=Path(args.stage02a).resolve();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    refs=Path(args.references).resolve() if args.references else Path(__file__).resolve().parents[1]/"references"/"REFERENCES_VERIFIED.tsv"
    qsum=json.loads((qc/"QC_SUMMARY.json").read_text());asum=json.loads((s2a/"STAGE02A_SUMMARY.json").read_text());error=None;res={}
    try:
        for t in ("MAPK","ADK"):res[t]=target_run(t,root,qc,s2,out)
    except Exception as e:error=f"{type(e).__name__}: {e}"
    p1=qsum.get("final_status")=="PASS";p2=asum.get("final_status")=="PASS" and asum.get("apo_protein_metrics_eligible") is True
    p3=error is None and set(res)=={"MAPK","ADK"}
    p4=p3 and res["MAPK"][6]["native_stride"]==1 and res["ADK"][6]["native_stride"]==5
    p5=p3 and all(len(res[t][0])==10000 for t in res)
    p6=p3 and all(abs(res[t][5][0]-.01)<1e-8 and abs(res[t][5][-1]-100)<1e-8 and np.allclose(np.diff(res[t][5]),.01,atol=1e-8,rtol=0) for t in res)
    p7=p3 and all(all(np.isfinite(float(r[k])) and float(r[k])>=0 for r in res[t][0] for k in ("ca_rmsd_A","backbone_rmsd_A","protein_rg_A")) for t in res)
    p8=p3 and all(all(np.isfinite(float(r["ca_rmsf_A"])) and float(r["ca_rmsf_A"])>=0 for r in res[t][1]) for t in res)
    p9=p3 and all([res[t][2][0]["n_frames"]==5000 and res[t][2][1]["n_frames"]==5000 for t in res])
    refrows=read_tsv(refs) if refs.exists() else [];p10=bool(refrows) and all(r.get("title") and r.get("doi") and r.get("pmid") and r.get("verified_status")=="VERIFIED" for r in refrows)
    passes=[("P1_STAGE00_QC_PASS",p1),("P2_TOPOLOGY_MODE_AUDIT_PASS",p2),("P3_BOTH_TARGET_ANALYSES_COMPLETE",p3),("P4_TARGET_STRIDES_MAPK1_ADK5",p4),("P5_EXACTLY_10000_FRAMES_EACH",p5),("P6_COMMON_10PS_TIME_GRID_0p01_TO_100NS",p6),("P7_RMSD_RG_NUMERIC_SANITY",p7),("P8_RMSF_NUMERIC_SANITY",p8),("P9_TWO_5000_FRAME_TEMPORAL_BLOCKS",p9),("P10_REFERENCE_IDENTITY_MANIFEST",p10)]
    summary={"error":error,"passes":[{"id":k,"status":"PASS" if v else "FAIL"} for k,v in passes],"final_status":"PASS" if all(v for _,v in passes) else "FAIL","next_stage_eligible":bool(all(v for _,v in passes)),
             "sampling_policy":"Matched 10-ps production sampling: MAPK native 10 ps stride 1; ADK native 2 ps stride 5. Raw DCD files are read-only. No interpolation, smoothing, filtering, imputation or outlier removal."}
    (out/"STAGE03_MATCHED10PS_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    with (out/"10PASS_REPORT.txt").open("w") as f:
        f.write("NAMD CORDYCEPIN MAPK/ADK — 03 MATCHED 10-PS APO METRICS 10-PASS REPORT\n"+"="*94+"\n")
        for k,v in passes:f.write(f"{k:<66} {'PASS' if v else 'FAIL'}\n")
        f.write("-"*94+f"\nFINAL STATUS: {summary['final_status']}\nNEXT-STAGE ELIGIBLE: {'YES' if summary['next_stage_eligible'] else 'NO'}\n")
        if error:f.write("ERROR: "+error+"\n")
    print(json.dumps(summary,indent=2));return 0 if all(v for _,v in passes) else 2
if __name__=="__main__":raise SystemExit(main())
