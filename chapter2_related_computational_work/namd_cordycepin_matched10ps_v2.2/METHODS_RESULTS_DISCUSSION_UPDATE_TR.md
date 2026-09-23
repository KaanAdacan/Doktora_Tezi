# MATCHED 10-PS ANALYSIS — TEZ METNİ GÜNCELLEME NOTU

## Metot — değiştirilecek
MAPK ve ADK üretim trajelerinin farklı kayıt aralıkları nedeniyle analizler ortak 10 ps zaman çözünürlüğünde harmonize edilmiştir. MAPK trajesi 10 ps/frame kaydedildiğinden stride 1 ile bütün üretim kareleri kullanılmıştır. ADK trajesi 2 ps/frame kaydedildiğinden her beşinci kayıtlı üretim karesi seçilmiş (stride 5) ve böylece her iki hedef için 0.01–100.00 ns aralığında 10,000 analiz karesi elde edilmiştir. Bu seçim deterministik temporal subsampling olup ham DCD dosyaları değiştirilmemiş; interpolasyon uygulanmamıştır, smoothing, filtering, imputation veya outlier çıkarma yapılmamıştır. RMSD, RMSF, Rg, PCA ve cluster-quality analizleri aynı matched 10-ps frame setleri üzerinden yeniden hesaplanacaktır.

## Şimdiden kesin olarak yeniden türetilebilen per-frame sonuçlar
Mevcut per-frame Stage03 tablolarından yapılan exact 10-ps seçimde MAPK değişmemekte (stride 1, 10,000 frame), ADK ise raw-frame 5,10,15,...,50000 seçimiyle 10,000 frame'e indirilmektedir.

- MAPK Cα RMSD: 2.466786 ± 0.325251 Å; 0–50 ns 2.607836 Å; 50–100 ns 2.325736 Å.
- MAPK Rg: 21.989452 ± 0.118745 Å; 0–50 ns 21.973897 Å; 50–100 ns 22.005008 Å.
- ADK Cα RMSD (matched 10 ps): 3.324799 ± 0.602123 Å; 0–50 ns 2.900193 Å; 50–100 ns 3.749406 Å.
- ADK Rg (matched 10 ps): 20.878410 ± 0.223791 Å; 0–50 ns 20.922989 Å; 50–100 ns 20.833830 Å.

## Yeniden hesaplanması gereken sonuçlar
RMSF ve PCA frame-ensemble bağımlı analizlerdir. Bunların matched 10-ps değerleri raw DCD üzerinde yeniden hesaplanmalıdır; all-native Stage03/04 değerleri final matched-10ps tez metnine doğrudan taşınmamalıdır.

## Tartışma
RMSD/Rg üzerinden ana nitel yorum değişmemektedir: MAPK geç yarıda daha düşük reference-relative RMSD rejimi örneklerken ADK ikinci yarıda T1 referansından daha uzak bir rejime geçmektedir; ADK Rg'nin eşzamanlı monoton genişleme göstermemesi global unfolding yerine internal conformational rearrangement yorumuyla daha uyumludur.
