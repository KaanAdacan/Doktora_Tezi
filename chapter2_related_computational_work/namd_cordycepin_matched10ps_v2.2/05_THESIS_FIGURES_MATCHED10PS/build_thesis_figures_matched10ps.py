#!/usr/bin/env python3
from pathlib import Path
import argparse, pandas as pd, numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage03',default='/home/kaan/NAMD/Publication_Analysis/03_APO_PROTEIN_METRICS_MATCHED10PS');ap.add_argument('--stage04',default='/home/kaan/NAMD/Publication_Analysis/04_PCA_CLUSTERING_MATCHED10PS');ap.add_argument('--out',default='/home/kaan/NAMD/Publication_Analysis/05_THESIS_FIGURES_MATCHED10PS');args=ap.parse_args();s3=Path(args.stage03);s4=Path(args.stage04);out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    # FIG02: MAPK first, ADK second, no smoothing.
    fig,axs=plt.subplots(2,3,figsize=(15,8.5))
    for r,target in enumerate(['MAPK','ADK']):
        m=pd.read_csv(s3/target/'PROTEIN_STRUCTURAL_METRICS_MATCHED10PS.tsv',sep='\t');rf=pd.read_csv(s3/target/'CA_RMSF_MATCHED10PS.tsv',sep='\t')
        axs[r,0].plot(m.time_ns,m.ca_rmsd_A,linewidth=.7);axs[r,0].set_xlabel('Time (ns)');axs[r,0].set_ylabel('Cα RMSD (Å)');axs[r,0].set_title(f'{target} — Cα RMSD')
        axs[r,1].plot(rf.resid,rf.ca_rmsf_A,linewidth=.7);axs[r,1].set_xlabel('Residue number');axs[r,1].set_ylabel('Cα RMSF (Å)');axs[r,1].set_title(f'{target} — Cα RMSF')
        axs[r,2].plot(m.time_ns,m.protein_rg_A,linewidth=.7);axs[r,2].set_xlabel('Time (ns)');axs[r,2].set_ylabel('Rg (Å)');axs[r,2].set_title(f'{target} — Rg')
    for i,ax in enumerate(axs.ravel()):ax.text(-.13,1.04,'ABCDEF'[i],transform=ax.transAxes,fontweight='bold',fontsize=13);ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(out/'FIG02_ApoMD_MATCHED10PS_MAPK_then_ADK.png',dpi=600,bbox_inches='tight');fig.savefig(out/'FIG02_ApoMD_MATCHED10PS_MAPK_then_ADK.pdf',bbox_inches='tight');plt.close(fig)
    # FIG03 matched PCA.
    fig,axs=plt.subplots(2,2,figsize=(11,9))
    for r,target in enumerate(['MAPK','ADK']):
        v=pd.read_csv(s4/target/'PCA_VARIANCE_MATCHED10PS.tsv',sep='\t').iloc[:20];p=pd.read_csv(s4/target/'PCA_SCORES_MATCHED10PS.tsv',sep='\t')
        axs[r,0].plot(v.pc,100*v.explained_variance_ratio,marker='o',markersize=2.5,linewidth=.8);axs[r,0].set_xlabel('Principal component');axs[r,0].set_ylabel('Explained variance (%)');axs[r,0].set_title(f'{target} — PCA variance')
        sc=axs[r,1].scatter(p.PC1,p.PC2,c=p.time_ns,s=3,alpha=.55,rasterized=True);axs[r,1].set_xlabel('PC1');axs[r,1].set_ylabel('PC2');axs[r,1].set_title(f'{target} — PC1/PC2');fig.colorbar(sc,ax=axs[r,1],label='Time (ns)')
    for i,ax in enumerate(axs.ravel()):ax.text(-.13,1.04,'ABCD'[i],transform=ax.transAxes,fontweight='bold',fontsize=13);ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(out/'FIG03_PCA_MATCHED10PS_MAPK_then_ADK.png',dpi=600,bbox_inches='tight');fig.savefig(out/'FIG03_PCA_MATCHED10PS_MAPK_then_ADK.pdf',bbox_inches='tight');plt.close(fig)
    print(out)
if __name__=='__main__':main()
