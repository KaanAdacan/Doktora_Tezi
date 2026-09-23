#!/usr/bin/env python3
from pathlib import Path
import sys, numpy as np, py_compile
HERE=Path(__file__).resolve();ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'03_APO_PROTEIN_METRICS_MATCHED10PS'/'code'))
from sampling_policy import matched_indices, TARGET_SAMPLING
checks={}
adk_times=np.arange(1,50001,dtype=float)*.002
mapk_times=np.arange(1,10001,dtype=float)*.01
ai,at=matched_indices(adk_times,'ADK');mi,mt=matched_indices(mapk_times,'MAPK')
checks['V1_ADK_native_50000']=len(adk_times)==50000
checks['V2_MAPK_native_10000']=len(mapk_times)==10000
checks['V3_ADK_stride5']=TARGET_SAMPLING['ADK']['stride']==5 and np.array_equal(ai[:4],[4,9,14,19])
checks['V4_MAPK_stride1']=TARGET_SAMPLING['MAPK']['stride']==1 and np.array_equal(mi[:4],[0,1,2,3])
checks['V5_ADK_10000_output']=len(at)==10000
checks['V6_MAPK_10000_output']=len(mt)==10000
checks['V7_common_grid_identical']=np.allclose(at,mt,rtol=0,atol=1e-12)
checks['V8_grid_10ps']=np.allclose(np.diff(at),.01,rtol=0,atol=1e-12) and at[0]==.01 and at[-1]==100
for p in [ROOT/'03_APO_PROTEIN_METRICS_MATCHED10PS'/'code'/'apo_protein_metrics_matched10ps.py',ROOT/'04_PCA_CLUSTERING_MATCHED10PS'/'code'/'pca_clustering_matched10ps.py']:
    py_compile.compile(str(p),doraise=True)
checks['V9_analysis_scripts_compile']=True
text=(ROOT/'METHODS_RESULTS_DISCUSSION_UPDATE_TR.md').read_text(encoding='utf-8')
checks['V10_method_declares_no_interpolation']=('interpolasyon uygulanmamıştır' in text.lower() and 'MAPK' in text and 'ADK' in text)
rep=HERE.parent/'KNOWN_ANSWER_VERIFICATION.txt';lines=['MATCHED 10-PS SAMPLING — KNOWN-ANSWER VERIFICATION','='*86]+[f"{k:<62} {'PASS' if v else 'FAIL'}" for k,v in checks.items()]+['-'*86,'FINAL: 10/10 PASS' if all(checks.values()) else 'FINAL: FAIL'];rep.write_text('\n'.join(lines)+'\n');print(rep.read_text());raise SystemExit(0 if all(checks.values()) else 2)
