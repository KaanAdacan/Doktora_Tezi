MAPK moleküler dinamik kodları

namd_config/ dizini NVT, NPT ve 0–100 ns üretim simülasyonu için altı aşamalı NAMD yapılandırma setini içerir. Bu yapılandırmalar doğrulanmış `MAPK_NVT_NPT_GPU(2).zip` kaynak paketiyle birebir senkronize edilmiştir; paket SHA-256 değeri kökteki `SOURCE_ARCHIVE_SHA256.txt` dosyasında kayıtlıdır.

runtime_scripts/ dizini yerel WSL/GPU zincir çalıştırıcısı, doğrulama, durum izleme, canlı izleme ve durdurma yardımcılarını içerir. Kaynak paketin orijinal kullanım notu `README_WSL_GPU.txt`, dosya SHA-256 envanteri `SHA256SUMS_LOCAL.txt` olarak aynı dizinde tutulur.

analysis/ dizinindeki standalone RMSD/RMSF/Rg/kalıntı dinamikleri/snapshot kodları provenance amacıyla korunur. Tezde raporlanan final matched-10ps MAPK/ADK metrikleri için canonical analiz zinciri `../../namd_cordycepin_matched10ps_v2.2/` dizinidir.

Beklenen harici girdiler:
inputs/mapk_ion.psf
inputs/mapk_ion.pdb
inputs/par_water_ions_ale.prm
inputs/par_all36_prot.prm
inputs/par_all36_cgenff.prm

Büyük girdi dosyaları bu kamuya açık kod deposuna bilinçli olarak dahil edilmemiştir.
