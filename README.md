# Doktora Tezi — Kaynak Kod ve Yeniden Üretilebilirlik Materyalleri

Bu açık depo, doktora tezi kapsamında kamuya açık tutulmasına karar verilen kaynak kod ve yeniden üretilebilirlik materyallerini içermektedir.

## Bölüm 2 — Moleküler dinamik ve ilişkili hesaplamalı çalışmalar

Raporlanan MAPK/ADK yapısal MD metrikleri için canonical analiz zinciri:
`chapter2_related_computational_work/namd_cordycepin_matched10ps_v2.2/`

Bu zincirde production trajeleri ortak 10 ps zaman çözünürlüğünde analiz edilmektedir: MAPK stride 1, ADK stride 5. Raporlanan temel çıktılar Cα RMSD, backbone RMSD, Cα RMSF, tüm-protein Rg, PCA/RMSIP ve cluster-quality diagnostikleridir. NVT/NPT raporlanan yapısal observabllara dahil edilmez; ham DCD dosyaları değiştirilmez ve interpolasyon, smoothing, filtering, imputation veya outlier removal uygulanmaz.

NAMD MAPK/ADK çalışma yapılandırmaları ve runtime yardımcıları `chapter2_related_computational_work/molecular_dynamics/` altındadır. MAPK NAMD yapılandırmaları doğrulanmış yerel GPU kaynak paketiyle eşitlenmiştir. Bu dizindeki eski standalone `analysis/` yardımcıları provenance amacıyla korunmaktadır; raporlanan final matched-10ps metrikler için canonical kaynak yukarıdaki `namd_cordycepin_matched10ps_v2.2/` dizinidir.

## Bölüm 3 — Cyx-KA CohortMaster ve Cyx-KA Vesseller

Kaynak kodlar yayın öncesi kapalıdır; kamuya açık benchmark/validation materyalleri `chapter3_tools/benchmarks/` altında tutulmaktadır. Bu senkronizasyonda doğrulanmış daha yeni bir Chapter 3 kaynak kod paketi bulunmadığından çekirdek araç kodları eklenmemiş veya değiştirilmemiştir.

Bu depoda açık kaynak lisansı tanımlanmamıştır.
