# Chapter 3 Benchmark ve Doğrulama Materyalleri

Bu paket Cyx-KA CohortMaster ve Cyx-KA Vesseller için tez kapsamında doğrulanan benchmark materyallerini içerir. Araştırma yazılımlarının kaynak kodları bu pakete dahil değildir.

## CohortMaster
- GEO10: 10 GSE, 11 series-matrix bileşeni; tüm raw dosyalar rapordaki SHA-256 kayıtlarıyla birebir eşleşmektedir.
- TCGA-GBM/PTEN: tez-run kabulü ve bağımsız R karşılaştırması PASS.
- GSE143383/IGF1R: 7/7 feature için tez-run kabulü ve bağımsız R karşılaştırması PASS.

## Vesseller
- 38 özgün hold-out CAM görüntüsü ve 38 frozen reference mask.
- 38/38 per-image değerlendirme ve 38/38 four-way SHA-256 kimlik eşleşmesi.
- Final export kontrolü: 161 PASS, 0 FAIL, 0 MISSING.
- Yazılım sürümü: 1.1.3.2.

## Büyük dosyalar
Raw GEO matrix arşivleri ve doğrulama ZIP'leri Git geçmişine eklenmemeli; `PUBLIC_RELEASE_ASSET_MANIFEST.tsv` içindeki SHA-256 değerleriyle GitHub Release asset olarak saklanmalıdır.
