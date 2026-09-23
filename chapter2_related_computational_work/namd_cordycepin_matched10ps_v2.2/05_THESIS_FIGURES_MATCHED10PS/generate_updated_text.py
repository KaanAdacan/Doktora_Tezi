#!/usr/bin/env python3
from pathlib import Path
import argparse,json,pandas as pd,numpy as np

def f(x,n=3):return f"{float(x):.{n}f}"
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage03',default='/home/kaan/NAMD/Publication_Analysis/03_APO_PROTEIN_METRICS_MATCHED10PS');ap.add_argument('--stage04',default='/home/kaan/NAMD/Publication_Analysis/04_PCA_CLUSTERING_MATCHED10PS');ap.add_argument('--out',default='/home/kaan/NAMD/Publication_Analysis/05_THESIS_FIGURES_MATCHED10PS/UPDATED_METHODS_RESULTS_DISCUSSION_TR.md');args=ap.parse_args();s3=Path(args.stage03);s4=Path(args.stage04)
    d={}
    for t in ['MAPK','ADK']:
        m=pd.read_csv(s3/t/'PROTEIN_STRUCTURAL_METRICS_MATCHED10PS.tsv',sep='\t');r=pd.read_csv(s3/t/'CA_RMSF_MATCHED10PS.tsv',sep='\t');b=pd.read_csv(s3/t/'BLOCK_0_50_vs_50_100_MATCHED10PS.tsv',sep='\t');v=pd.read_csv(s4/t/'PCA_VARIANCE_MATCHED10PS.tsv',sep='\t');sm=json.loads((s4/t/'ANALYSIS_MANIFEST.json').read_text())
        d[t]={'rmsd_mean':m.ca_rmsd_A.mean(),'rmsd_sd':m.ca_rmsd_A.std(ddof=1),'rg_mean':m.protein_rg_A.mean(),'rg_sd':m.protein_rg_A.std(ddof=1),'rmsf_mean':r.ca_rmsf_A.mean(),'rmsf_median':r.ca_rmsf_A.median(),'rmsf_max':r.ca_rmsf_A.max(),'early_rmsd':float(b.loc[b.block_ns=='0-50','ca_rmsd_A_mean'].iloc[0]),'late_rmsd':float(b.loc[b.block_ns=='50-100','ca_rmsd_A_mean'].iloc[0]),'early_rg':float(b.loc[b.block_ns=='0-50','protein_rg_A_mean'].iloc[0]),'late_rg':float(b.loc[b.block_ns=='50-100','protein_rg_A_mean'].iloc[0]),'pc1':100*v.loc[0,'explained_variance_ratio'],'pc2':100*v.loc[1,'explained_variance_ratio'],'top5':100*v.loc[4,'cumulative_variance'],'top10':100*v.loc[9,'cumulative_variance'],'rmsip':sm['half_trajectory_RMSIP_top3'],'k_ch':sm['k_by_max_calinski_harabasz'],'k_db':sm['k_by_min_davies_bouldin'],'selected_k':sm['selected_k']}
    text=f'''# GÜNCELLENMİŞ YÖNTEM, SONUÇ VE TARTIŞMA — MATCHED 10-PS ANALYSIS

## Yöntem güncellemesi
MAPK ve ADK üretim trajelerinin farklı kayıt aralıklarına sahip olması nedeniyle analizler ortak 10 ps zaman çözünürlüğünde harmonize edilmiştir. MAPK trajesi doğal olarak 10 ps/frame kaydedildiğinden stride 1 ile bütün kayıtlı üretim kareleri kullanılmıştır. ADK trajesi 2 ps/frame kaydedildiğinden her beşinci kayıtlı üretim karesi seçilmiş (stride 5) ve böylece her iki hedef için 0.01–100.00 ns aralığında 10,000 analiz karesi elde edilmiştir. Seçim deterministik olup ham DCD dosyaları değiştirilmemiştir; interpolasyon, smoothing, filtering, imputation veya outlier çıkarma uygulanmamıştır. RMSD, RMSF, Rg, PCA ve cluster-quality değerlendirmeleri bu aynı matched 10-ps frame setleri üzerinde yeniden hesaplanmıştır. T1/T2/T3 exact boundary states ve docking analizleri bu örnekleme değişikliğinden etkilenmemiştir.

## Sonuçlar
MAPK için matched-10ps Cα RMSD {f(d['MAPK']['rmsd_mean'])} ± {f(d['MAPK']['rmsd_sd'])} Å olup 0–50 ns ve 50–100 ns ortalamaları sırasıyla {f(d['MAPK']['early_rmsd'])} ve {f(d['MAPK']['late_rmsd'])} Å'dır. Ortalama Rg {f(d['MAPK']['rg_mean'])} ± {f(d['MAPK']['rg_sd'])} Å; Cα RMSF ortalaması {f(d['MAPK']['rmsf_mean'])} Å ve medyanı {f(d['MAPK']['rmsf_median'])} Å'dır. PCA'da PC1 ve PC2 sırasıyla %{f(d['MAPK']['pc1'],2)} ve %{f(d['MAPK']['pc2'],2)} varyans açıklamış; ilk 10 PC kümülatif %{f(d['MAPK']['top10'],2)} varyans açıklamıştır. Half-trajectory top-3 RMSIP {f(d['MAPK']['rmsip'])}'dir.

ADK için matched-10ps Cα RMSD {f(d['ADK']['rmsd_mean'])} ± {f(d['ADK']['rmsd_sd'])} Å olup 0–50 ns ve 50–100 ns ortalamaları sırasıyla {f(d['ADK']['early_rmsd'])} ve {f(d['ADK']['late_rmsd'])} Å'dır. Ortalama Rg {f(d['ADK']['rg_mean'])} ± {f(d['ADK']['rg_sd'])} Å; Cα RMSF ortalaması {f(d['ADK']['rmsf_mean'])} Å ve medyanı {f(d['ADK']['rmsf_median'])} Å'dır. PCA'da PC1 ve PC2 sırasıyla %{f(d['ADK']['pc1'],2)} ve %{f(d['ADK']['pc2'],2)} varyans açıklamış; ilk 10 PC kümülatif %{f(d['ADK']['top10'],2)} varyans açıklamıştır. Half-trajectory top-3 RMSIP {f(d['ADK']['rmsip'])}'dir.

## Tartışma güncellemesi
Matched 10-ps yeniden analiz, iki hedef arasındaki görsel örnekleme yoğunluğunu eşitlemekte ve karşılaştırılabilir bir zaman çözünürlüğü sağlamaktadır. Bu işlem veri normalizasyonu veya smoothing değildir; ADK için önceden belirlenmiş sabit aralıklı deterministik temporal subsampling'dir. MAPK ve ADK için bütün raporlanan MD metrikleri aynı 10 ps analiz gridinden türetilmiştir. Sonuçların yorumu yine tek bağımsız 100 ns apo trajectory sınırı içerisinde tutulmalı; temporal bloklar bağımsız replicate olarak değerlendirilmemelidir.

Cluster-quality sonuçları: MAPK CH k={d['MAPK']['k_ch']}, DB k={d['MAPK']['k_db']}, selected_k={d['MAPK']['selected_k']}; ADK CH k={d['ADK']['k_ch']}, DB k={d['ADK']['k_db']}, selected_k={d['ADK']['selected_k']}. İki kriter uyuşmadığında tek bir k zorlanmamalıdır.
'''
    p=Path(args.out);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text,encoding='utf-8');print(p)
if __name__=='__main__':main()
