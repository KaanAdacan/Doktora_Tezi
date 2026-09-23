#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,re,hashlib,sys
from pathlib import Path
from collections import Counter
import numpy as np
import MDAnalysis as mda
from sklearn.cluster import KMeans
from sklearn.metrics import calinski_harabasz_score,davies_bouldin_score
sys.path.insert(0,str(Path(__file__).resolve().parent))
from sampling_policy import TARGET_SAMPLING, matched_indices
RANDOM_STATE=20260923;MAX_PCS=50;CLUSTER_KS=tuple(range(2,7))

def read_tsv(p):
    with Path(p).open(encoding="utf-8") as f:return list(csv.DictReader(f,delimiter="\t"))
def write_tsv(p,rows,fields):
    with Path(p).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore");w.writeheader();[w.writerow({k:r.get(k,"") for k in fields}) for r in rows]
def sha256_file(p):
    h=hashlib.sha256();
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def _production_identity(rel):
    p=Path(rel);local=(p.parent.name+"/"+p.name).upper();return ("PROD_" in local) and not bool(re.search(r"(^|[/_.-])(NVT|NPT|BENCH)([/_.-]|$)",local))
def build_axis(chain,target,root):
    rr=sorted([r for r in chain if r["target"]==target],key=lambda r:float(r["actual_start_ns"]));paths=[];times=[];prev=0.0
    for i,r in enumerate(rr):
        rel=r["relative_path"]
        if not _production_identity(rel):raise ValueError(f"Non-production DCD: {rel}")
        a=float(r["actual_start_ns"]);b=float(r["actual_end_ns"]);n=int(r["nset"])
        if i==0 and abs(a)>1e-9:raise ValueError("Production does not start at 0")
        if i and abs(a-prev)>1e-9:raise ValueError("Production chain discontinuity")
        dt=(b-a)/n;times.extend(a+(j+1)*dt for j in range(n));p=(root/rel).resolve();
        if not p.exists():raise FileNotFoundError(p)
        paths.append(str(p));prev=b
    if abs(prev-100)>1e-8:raise ValueError("Production does not end at 100 ns")
    return paths,np.asarray(times,dtype=np.float64),rr
def kabsch(mobile,reference):
    mob=np.asarray(mobile,dtype=np.float64);ref=np.asarray(reference,dtype=np.float64);mc=mob.mean(0);rc=ref.mean(0);P=mob-mc;Q=ref-rc;C=P.T@Q;V,S,Wt=np.linalg.svd(C,full_matrices=False);D=np.eye(3);D[2,2]=-1.0 if np.linalg.det(V@Wt)<0 else 1.0;R=V@D@Wt;return P@R+rc,R,mc,rc
def exact_pca(X,max_pcs=MAX_PCS):
    X=np.asarray(X,dtype=np.float64);mean=X.mean(0);Y=X-mean;cov=(Y.T@Y)/(len(Y)-1);vals,vecs=np.linalg.eigh(cov);order=np.argsort(vals)[::-1];vals=np.clip(vals[order],0,None);vecs=vecs[:,order];total=float(vals.sum());
    if total<=0:raise ValueError("PCA total variance is zero")
    ncomp=min(max_pcs,len(vals),len(X)-1);vals=vals[:ncomp];comps=vecs[:,:ncomp].T;ev=vals/total;scores=Y@comps.T;return mean,comps,ev,scores
def exact_pca_components(X,n=3):
    X=np.asarray(X,dtype=np.float64);Y=X-X.mean(0);cov=(Y.T@Y)/(len(Y)-1);vals,vecs=np.linalg.eigh(cov);order=np.argsort(vals)[::-1];return vecs[:,order[:n]].T
def rmsip(A,B,n=3):return float(np.sqrt(np.sum((np.asarray(A[:n])@np.asarray(B[:n]).T)**2)/n))
def plot_scree(ev,png,pdf,target):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt;x=np.arange(1,len(ev)+1);fig,ax=plt.subplots(figsize=(7.2,4.6));ax.plot(x,100*np.asarray(ev),marker="o",markersize=2.5,linewidth=1);ax.set_xlabel("Principal component");ax.set_ylabel("Explained variance (%)");ax.set_title(f"{target} — PCA, matched 10-ps sampling");fig.tight_layout();fig.savefig(png,dpi=600,bbox_inches="tight");fig.savefig(pdf,bbox_inches="tight");plt.close(fig)
def plot_pc(scores,times,png,pdf,target,labels=None):
    import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt;fig,ax=plt.subplots(figsize=(6.4,5.2));c=times if labels is None else labels;sc=ax.scatter(scores[:,0],scores[:,1],c=c,s=3,alpha=.55,rasterized=True);ax.set_xlabel("PC1");ax.set_ylabel("PC2");ax.set_title(f"{target} — PC1 vs PC2, matched 10-ps sampling");fig.colorbar(sc,ax=ax,label="Production time (ns)" if labels is None else "Cluster");fig.tight_layout();fig.savefig(png,dpi=600,bbox_inches="tight");fig.savefig(pdf,bbox_inches="tight");plt.close(fig)
def write_aligned_protein_frame(u,raw_frame_index,ca,ref_ca,protein,outp):
    u.trajectory[int(raw_frame_index)];_,R,mc,rc=kabsch(ca.positions,ref_ca);original=protein.positions.copy();protein.positions=((np.asarray(original,dtype=np.float64)-mc)@R+rc).astype(np.float32);protein.write(str(outp));protein.positions=original

def analyze_target(target,root,qc,s2,s2a,out):
    chain=read_tsv(qc/"PRODUCTION_CHAIN.tsv");snaps=read_tsv(s2/"SNAPSHOT_MANIFEST.tsv");s2sum=json.loads((s2/"STAGE02_SUMMARY.json").read_text());a2=json.loads((s2a/"STAGE02A_SUMMARY.json").read_text())
    if a2.get("final_status")!="PASS" or a2.get("apo_protein_metrics_eligible") is not True:raise RuntimeError("Stage02A not eligible")
    chosen=s2sum["chosen_topology"][target];psf=(root/chosen["psf"]).resolve();t1=[r for r in snaps if r["target"]==target and r["label"]=="T1"]
    if len(t1)!=1:raise ValueError("T1 identity error")
    refp=Path(t1[0]["output_pdb"]).resolve();dcds,native_times,chain_rows=build_axis(chain,target,root);sample_idx,times=matched_indices(native_times,target);sample_set=set(int(i) for i in sample_idx)
    u=mda.Universe(str(psf),dcds,format="DCD");ref=mda.Universe(str(psf),str(refp));ca=u.select_atoms("protein and name CA");protein=u.select_atoms("protein");rca=ref.select_atoms("protein and name CA")
    if len(ca)==0 or len(ca)!=len(rca):raise ValueError("C-alpha selection mismatch")
    ref_ca=np.asarray(rca.positions.copy(),dtype=np.float64);X=np.empty((10000,len(ca)*3),dtype=np.float32);raw_indices=np.empty(10000,dtype=int);j=0
    for raw_i,ts in enumerate(u.trajectory):
        if raw_i not in sample_set:continue
        aligned,_,_,_=kabsch(ca.positions,ref_ca);X[j]=aligned.reshape(-1).astype(np.float32);raw_indices[j]=raw_i;j+=1
    if j!=10000:raise ValueError(f"{target}: PCA sampled {j}, expected 10000")
    pca_mean,components,ev,scores=exact_pca(X);cum=np.cumsum(ev);hit=np.flatnonzero(cum>=.80);n_cluster_pc=int(hit[0]+1) if len(hit) else min(10,len(ev));n_cluster_pc=max(2,min(10,n_cluster_pc));Z=np.asarray(scores[:,:n_cluster_pc],dtype=np.float64)
    cluster_tests=[];models={}
    for k in CLUSTER_KS:
        km=KMeans(n_clusters=k,n_init=20,random_state=RANDOM_STATE,algorithm="lloyd");lab=km.fit_predict(Z);ch=float(calinski_harabasz_score(Z,lab));db=float(davies_bouldin_score(Z,lab));cluster_tests.append({"target":target,"k":k,"calinski_harabasz":f"{ch:.8f}","davies_bouldin":f"{db:.8f}","inertia":f"{float(km.inertia_):.8f}"});models[k]=(km,lab,ch,db)
    k_ch=max(models,key=lambda k:(models[k][2],-k));k_db=min(models,key=lambda k:(models[k][3],k));selected_k=k_ch if k_ch==k_db else None;labels=np.full(10000,-1,dtype=int);counts=Counter();dominant=None;rep_pdb=None;rep_time=None;occ=[]
    if selected_k is not None:
        km,labels,_,_=models[selected_k];counts=Counter(labels);dominant=max(counts,key=lambda c:(counts[c],-c));dom=np.flatnonzero(labels==dominant);mean_xyz=np.asarray(X[dom],dtype=np.float64).mean(0);dist2=np.sum((np.asarray(X[dom],dtype=np.float64)-mean_xyz)**2,axis=1);rep_local=int(dom[int(np.argmin(dist2))]);rep_time=float(times[rep_local])
        for c in sorted(counts):
            alln=int(np.sum(labels==c));em=times<=50+1e-9;lm=times>50+1e-9;e=int(np.sum(labels[em]==c));l=int(np.sum(labels[lm]==c));occ.append({"target":target,"cluster":int(c),"overall_n":alln,"overall_pct":100*alln/10000,"early_0_50_n":e,"early_0_50_pct":100*e/5000,"late_50_100_n":l,"late_50_100_pct":100*l/5000,"dominant":bool(c==dominant)})
    early=X[times<=50+1e-9];late=X[times>50+1e-9];half_rmsip=rmsip(exact_pca_components(early,3),exact_pca_components(late,3),3)
    td=out/target;(td/"figures").mkdir(parents=True,exist_ok=True);(td/"representatives").mkdir(parents=True,exist_ok=True)
    srows=[]
    for i,t in enumerate(times):
        row={"target":target,"analysis_frame_index":i+1,"raw_frame_index":int(raw_indices[i])+1,"time_ns":f"{t:.8f}","cluster":int(labels[i])};
        for pc in range(min(10,scores.shape[1])):row[f"PC{pc+1}"]=f"{float(scores[i,pc]):.8f}"
        srows.append(row)
    write_tsv(td/"PCA_SCORES_MATCHED10PS.tsv",srows,["target","analysis_frame_index","raw_frame_index","time_ns","cluster"]+[f"PC{i}" for i in range(1,min(10,scores.shape[1])+1)])
    evrows=[{"target":target,"pc":i+1,"explained_variance_ratio":f"{float(ev[i]):.10f}","cumulative_variance":f"{float(cum[i]):.10f}"} for i in range(len(ev))];write_tsv(td/"PCA_VARIANCE_MATCHED10PS.tsv",evrows,["target","pc","explained_variance_ratio","cumulative_variance"])
    write_tsv(td/"CLUSTER_K_SELECTION_MATCHED10PS.tsv",cluster_tests,["target","k","calinski_harabasz","davies_bouldin","inertia"]);write_tsv(td/"CLUSTER_OCCUPANCY_MATCHED10PS.tsv",occ,["target","cluster","overall_n","overall_pct","early_0_50_n","early_0_50_pct","late_50_100_n","late_50_100_pct","dominant"])
    if selected_k is not None:
        rep_local=int(np.argmin([np.sum((np.asarray(X[i],dtype=np.float64)-np.asarray(X[np.flatnonzero(labels==dominant)],dtype=np.float64).mean(0))**2) if labels[i]==dominant else np.inf for i in range(10000)]));rep_time=float(times[rep_local]);rep_pdb=td/"representatives"/f"{target}_DOMINANT_CLUSTER_REP_{rep_time:.3f}ns_PROTEIN.pdb";write_aligned_protein_frame(u,raw_indices[rep_local],ca,ref_ca,protein,rep_pdb)
    plot_scree(ev[:min(20,len(ev))],td/"figures"/"PCA_SCREE.png",td/"figures"/"PCA_SCREE.pdf",target);plot_pc(scores,times,td/"figures"/"PC1_PC2_TIME.png",td/"figures"/"PC1_PC2_TIME.pdf",target)
    if selected_k is not None:plot_pc(scores,times,td/"figures"/"PC1_PC2_CLUSTER.png",td/"figures"/"PC1_PC2_CLUSTER.pdf",target,labels=labels)
    manifest={"target":target,"sampling_policy":"DETERMINISTIC_MATCHED_10PS","native_interval_ps":TARGET_SAMPLING[target]["native_interval_ns"]*1000,"native_stride":TARGET_SAMPLING[target]["stride"],"analysis_interval_ps":10.0,"native_frames":len(native_times),"analysis_frames":10000,"first_analysis_time_ns":float(times[0]),"last_analysis_time_ns":float(times[-1]),"pca_method":"exact covariance eigendecomposition on matched 10-ps frames","pca_components_saved":int(len(ev)),"clustering_pcs":int(n_cluster_pc),"clustering_cumulative_variance":float(cum[n_cluster_pc-1]),"k_by_max_calinski_harabasz":int(k_ch),"k_by_min_davies_bouldin":int(k_db),"selected_k":int(selected_k) if selected_k is not None else None,"selection_rule":"selected only when max-CH and min-DB agree; otherwise no k is forced","dominant_cluster":int(dominant) if dominant is not None else None,"dominant_cluster_occupancy_pct":float(100*counts[dominant]/10000) if dominant is not None else None,"dominant_cluster_representative_time_ns":rep_time,"half_trajectory_RMSIP_top3":float(half_rmsip),"representative_pdb":str(rep_pdb) if rep_pdb is not None else None,"representative_sha256":sha256_file(rep_pdb) if rep_pdb is not None else None,"raw_file_policy":"read-only","excluded_operations":["interpolation","smoothing","filtering","imputation","outlier removal"]}
    (td/"ANALYSIS_MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8");return manifest

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--root",default="/home/kaan/NAMD/systems");ap.add_argument("--qc",default="/home/kaan/NAMD/Publication_Analysis/00_QC");ap.add_argument("--stage02",default="/home/kaan/NAMD/Publication_Analysis/02_TOPOLOGY_LIGAND_SNAPSHOTS");ap.add_argument("--stage02a",default="/home/kaan/NAMD/Publication_Analysis/02A_TOPOLOGY_MODE_AUDIT");ap.add_argument("--stage03",default="/home/kaan/NAMD/Publication_Analysis/03_APO_PROTEIN_METRICS_MATCHED10PS");ap.add_argument("--out",default="/home/kaan/NAMD/Publication_Analysis/04_PCA_CLUSTERING_MATCHED10PS");ap.add_argument("--references",default=None);args=ap.parse_args()
    root=Path(args.root).resolve();qc=Path(args.qc).resolve();s2=Path(args.stage02).resolve();s2a=Path(args.stage02a).resolve();s3=Path(args.stage03).resolve();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True);refs=Path(args.references).resolve() if args.references else Path(__file__).resolve().parents[1]/"references"/"REFERENCES_VERIFIED.tsv";s3sum=json.loads((s3/"STAGE03_MATCHED10PS_SUMMARY.json").read_text()) if (s3/"STAGE03_MATCHED10PS_SUMMARY.json").exists() else {};error=None;results={}
    try:
        if s3sum.get("final_status")!="PASS" or s3sum.get("next_stage_eligible") is not True:raise RuntimeError("Matched-10ps Stage03 not eligible")
        for t in ("MAPK","ADK"):results[t]=analyze_target(t,root,qc,s2,s2a,out)
    except Exception as e:error=f"{type(e).__name__}: {e}"
    p1=s3sum.get("final_status")=="PASS" and s3sum.get("next_stage_eligible") is True;p2=error is None and set(results)=={"MAPK","ADK"};p3=p2 and results["MAPK"]["native_stride"]==1 and results["ADK"]["native_stride"]==5;p4=p2 and all(results[t]["analysis_frames"]==10000 and abs(results[t]["first_analysis_time_ns"]-.01)<1e-8 and abs(results[t]["last_analysis_time_ns"]-100)<1e-8 for t in results);p5=p2 and all(results[t]["pca_components_saved"]>=3 and np.isfinite(results[t]["clustering_cumulative_variance"]) for t in results);p6=p2 and all(2<=results[t]["k_by_max_calinski_harabasz"]<=6 and 2<=results[t]["k_by_min_davies_bouldin"]<=6 for t in results);p7=p2 and all(0<=results[t]["half_trajectory_RMSIP_top3"]<=1 for t in results);p8=p2 and all(results[t]["selected_k"] is None or (results[t]["representative_pdb"] and Path(results[t]["representative_pdb"]).exists()) for t in results);p9=p2 and all((out/t/"PCA_SCORES_MATCHED10PS.tsv").exists() and (out/t/"PCA_VARIANCE_MATCHED10PS.tsv").exists() for t in results);refrows=read_tsv(refs) if refs.exists() else [];p10=bool(refrows) and all(r.get("title") and r.get("doi") and r.get("pmid") and r.get("verified_status")=="VERIFIED" for r in refrows)
    passes=[("P1_MATCHED10PS_STAGE03_PASS",p1),("P2_BOTH_TARGET_ANALYSES_COMPLETE",p2),("P3_TARGET_STRIDES_MAPK1_ADK5",p3),("P4_10000_FRAMES_EACH_COMMON_10PS_GRID",p4),("P5_EXACT_PCA_VALID",p5),("P6_ALL_FRAME_CLUSTER_METRICS_VALID",p6),("P7_HALF_TRAJECTORY_RMSIP_VALID",p7),("P8_NO_FORCED_REPRESENTATIVE_POLICY",p8),("P9_OUTPUT_INTEGRITY",p9),("P10_REFERENCE_IDENTITY_MANIFEST",p10)];summary={"error":error,"passes":[{"id":k,"status":"PASS" if v else "FAIL"} for k,v in passes],"final_status":"PASS" if all(v for _,v in passes) else "FAIL","next_stage_eligible":bool(all(v for _,v in passes)),"sampling_policy":"Matched 10-ps production sampling: MAPK stride 1, ADK stride 5; 10,000 frames each; no interpolation or smoothing.","targets":results};(out/"STAGE04_MATCHED10PS_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    with (out/"10PASS_REPORT.txt").open("w") as f:
        f.write("NAMD CORDYCEPIN MAPK/ADK — 04 MATCHED 10-PS PCA/CLUSTERING 10-PASS REPORT\n"+"="*96+"\n");[f.write(f"{k:<68} {'PASS' if v else 'FAIL'}\n") for k,v in passes];f.write("-"*96+f"\nFINAL STATUS: {summary['final_status']}\nNEXT-STAGE ELIGIBLE: {'YES' if summary['next_stage_eligible'] else 'NO'}\n");
        if error:f.write("ERROR: "+error+"\n")
    print(json.dumps(summary,indent=2));return 0 if all(v for _,v in passes) else 2
if __name__=="__main__":raise SystemExit(main())
