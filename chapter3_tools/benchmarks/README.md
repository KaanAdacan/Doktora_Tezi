# Chapter 3 Benchmark ve Doğrulama Materyalleri

Bu dizin, Cyx-KA CohortMaster ve Cyx-KA Vesseller için doktora tezi kapsamında doğrulanan benchmark, validation ve provenance materyallerini içerir. Araştırma yazılımlarının kaynak kodları bu public benchmark paketine dahil değildir.

## CohortMaster
- GEO10: 10 GSE, 11 series-matrix bileşeni; raw benchmark arşivleri GitHub Release asset olarak tutulur.
- TCGA-GBM/PTEN: tez-run kabulü ve bağımsız R karşılaştırması PASS.
- GSE143383/IGF1R: 7/7 feature için tez-run kabulü ve bağımsız R karşılaştırması PASS.
- Doğrulama özeti: `CohortMaster_VALIDATION_SUMMARY.tsv`.

## CAM — 38 görüntülük locked benchmark
Bu public CAM accuracy benchmarkı 38 görüntülük bağımsız locked hold-out setini kapsar.

### Cyx-KA Vesseller
- 38 özgün hold-out CAM görüntüsü ve 38 frozen reference mask.
- 38/38 başarılı analiz; rejected input=0, analysis failure=0.
- 38/38 per-image değerlendirme ve 38/38 four-way SHA-256 kimlik eşleşmesi.
- Yazılım sürümü: 1.1.3.2.
- Sayısal benchmark özeti: `Vesseller_38_LOCKED_BENCHMARK_SUMMARY.tsv`.
- Public frozen benchmark arşivi: GitHub Release asset `Chapter3_Vesseller_BENCHMARK_PUBLIC.zip`.
- Release asset SHA-256: `b63fbbc9d1eec8bf5163c3b85edffc78de74c6df6b7e1e1dc26d0d4db9e6b055`.

### AngioTool
- Aynı 38-image locked hold-out seti üzerinde 38/38 sonuç satırı doğrulandı.
- Public quantitative table: `AngioTool_38_LOCKED_RESULTS.tsv`.
- Locked run parameters: `AngioTool_38_RUN_PARAMETERS.tsv`.
- Kullanıcı makinesine ait absolute local path, tarih/saat ve gereksiz workbook metadata alanları public tablodan çıkarılmıştır; bilimsel ölçüm sütunları korunmuştur.

### NEREA
- Aynı 38-image benchmark için sağlanan arşivde 23 directory entry bulunmakta, payload dosyası bulunmamaktadır.
- Comparison-eligible quantitative çıktı oluşmadığı için accuracy hücreleri yapay olarak 0 verilmez; `N/A` olarak tutulur.
- Durum kaydı: `NEREA_38_BENCHMARK_STATUS.tsv`.

### Source provenance
Kullanıcı tarafından sağlanan 38-image comparator arşivlerinin byte-level SHA-256 kayıtları `COMPARATOR_38_SOURCE_ARCHIVE_SHA256.tsv` içinde tutulur. 38 original hold-out görüntüsü Vesseller public benchmark assetinde zaten bulunduğundan ayrı bir ikinci public kopya oluşturulmaz.

## CAM — 560 görüntülük operasyonel benchmark
560 görüntülük geniş CAM koleksiyonu accuracy benchmarkından ayrı olarak operasyonel tamamlama, kaynak-çıktı izlenebilirliği ve QC görünürlüğü açısından değerlendirilmiştir.

- Üç araçlı operasyonel özet: `CAM_560_OPERATIONAL_BENCHMARK_SUMMARY.tsv`.
- Cyx-KA Vesseller: 560/560 sonuç; QC PASS=423, QC REVIEW=137.
- AngioTool: 560/560 teknik sonuç satırı üretti, ancak locked run vasküler segmentasyon QC'sini karşılamadı; özet `AngioTool_LOCKED_RUN_QC.tsv` içindedir.
- NEREA: 463 sonuç satırı oluştu; 80 kaynak dosya eksik, 34 kaynak dosyada basename kimliği belirsiz ve temiz source-level finalizasyon için 114 dosyanın yeniden çalıştırılması gerekmektedir. Audit özeti `NEREA_560_AUDIT_SUMMARY.tsv` içindedir.
- Bu nedenle 560-image koleksiyon için üç araç arasında tamamlanmış comparison-eligible bir benchmark oluşmamıştır. Comparatorlarda geçerli reference-mask karşılaştırması bulunmadığında accuracy değerleri yapay olarak 0 verilmez; `N/A` olarak tutulur.
- 560 ham görüntü koleksiyonu public release asseti olarak paylaşılmamaktadır. Public depo, bu operasyonel çalışmanın yalnızca doğrulanmış özet/audit kayıtlarını taşır.

## Public release ve bütünlük
`chapter3-benchmarks-v1.0.0` release'i Vesseller frozen benchmark materyalini ve CohortMaster doğrulama varlıklarını taşır. Büyük GEO raw arşivleri `PUBLIC_RELEASE_ASSET_MANIFEST.tsv` içindeki SHA-256 değerleriyle GitHub Release asset olarak saklanır.

## Tez / Supplement kullanımı
Tezin EK/Tekrarlanabilirlik bölümünde immutable commit bağlantısı, Chapter 3 release bağlantısı ve ilgili SHA-256 kayıtları birlikte verilmelidir. Bu depo benchmark/doğrulama kanıtını public tutar; araştırma yazılımlarının kaynak kodları tez eklerinde ayrıca sunulmuştur.
