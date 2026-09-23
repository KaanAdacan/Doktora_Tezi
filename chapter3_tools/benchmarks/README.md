# Chapter 3 Benchmark ve Doğrulama Materyalleri

Bu dizin, Cyx-KA CohortMaster ve Cyx-KA Vesseller için doktora tezi kapsamında doğrulanan benchmark, audit ve provenance özetlerini içerir. Araştırma yazılımlarının kaynak kodları bu public benchmark paketine dahil değildir.

## CohortMaster
- GEO10: 10 GSE, 11 series-matrix bileşeni; tüm raw dosyalar rapordaki SHA-256 kayıtlarıyla birebir eşleşmektedir.
- TCGA-GBM/PTEN: tez-run kabulü ve bağımsız R karşılaştırması PASS.
- GSE143383/IGF1R: 7/7 feature için tez-run kabulü ve bağımsız R karşılaştırması PASS.

## Vesseller — 38 görüntülük locked benchmark
- 38 özgün hold-out CAM görüntüsü ve 38 frozen reference mask.
- 38/38 başarılı analiz; rejected input=0, analysis failure=0.
- 38/38 per-image değerlendirme ve 38/38 four-way SHA-256 kimlik eşleşmesi.
- Final export kontrolü: 161 PASS, 0 FAIL, 0 MISSING.
- Yazılım sürümü: 1.1.3.2.
- Sayısal benchmark özeti: `Vesseller_38_LOCKED_BENCHMARK_SUMMARY.tsv`.
- Public frozen benchmark arşivi: GitHub Release asset `Chapter3_Vesseller_BENCHMARK_PUBLIC.zip`.
- Release asset SHA-256: `b63fbbc9d1eec8bf5163c3b85edffc78de74c6df6b7e1e1dc26d0d4db9e6b055`.

## 560 görüntülük operasyonel benchmark
Bu veri seti 38 görüntülük locked accuracy benchmarkından ayrıdır ve batch tamamlama, kaynak-çıktı izlenebilirliği, QC görünürlüğü ve yeniden çalıştırma gereksinimi açısından değerlendirilir.

- Üç araçlı özet: `CAM_560_OPERATIONAL_BENCHMARK_SUMMARY.tsv`.
- Cyx-KA Vesseller: 560/560 sonuç; QC PASS=423, QC REVIEW=137.
- AngioTool: `angioanaliz.xlsx` 560/560 teknik sonuç satırı üretmiştir; ancak locked durum `QC-run/FAIL` olup görsel/skeleton QC bulguları `AngioTool_LOCKED_RUN_QC.tsv` içinde kaydedilmiştir.
- NEREA: 463 sonuç satırı; 80 eksik kaynak, 34 basename-kaynak kimliği belirsiz dosya ve temiz source-level finalizasyon için 114 yeniden çalıştırma gereksinimi. Audit özeti `NEREA_560_AUDIT_SUMMARY.tsv` içindedir.
- Comparator araçlarda geçerli reference-mask karşılaştırması oluşmadığı durumda accuracy hücreleri yapay olarak 0 verilmez; `N/A` olarak tutulur.

## Public release ve bütünlük
`chapter3-benchmarks-v1.0.0` release'i 38-image Vesseller frozen benchmark materyalini ve Chapter 3 public doğrulama varlıklarını taşır. Büyük raw GEO matrix arşivleri ve doğrulama ZIP'leri Git geçmişine eklenmemeli; `PUBLIC_RELEASE_ASSET_MANIFEST.tsv` içindeki SHA-256 değerleriyle GitHub Release asset olarak saklanmalıdır.

## Tez / Supplement kullanımı
Tezin EK/Tekrarlanabilirlik bölümünde:
1. bu dizinin immutable commit bağlantısı,
2. `chapter3-benchmarks-v1.0.0` GitHub Release bağlantısı,
3. Vesseller public benchmark ZIP dosya adı ve SHA-256 değeri
birlikte verilmelidir.

Bu yaklaşım tez metnindeki özet tabloları sürüm kontrollü kanıta bağlar; kaynak kodun public paylaşımını gerektirmez.
