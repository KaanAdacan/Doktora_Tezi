MAPK — WSL LOCAL NVIDIA GPU / NAMD 3.x
=====================================

HEDEF KLASOR
  ~/NAMD/systems/MAPK/NVT_NPT_GPU

AMAÇ
  MAPK_TRUBA_READY paketindeki bilimsel protokolu yerel WSL GPU'da,
  ADK yerel zincirindeki calistirma mantigiyla kosmak.

KORUNAN BILIMSEL AYARLAR
  - timestep: 2.0 fs
  - rigidBonds: all
  - nonbondedFreq: 1
  - fullElectFrequency: 2
  - stepspercycle: 20
  - cutoff/switch/pairlist: 12.0 / 10.0 / 13.5 A
  - PME grid: 256 x 256 x 256
  - NVT: 298 K, 50,000 minimization + 2,500,000 MD step = 5 ns
  - NPT: 298 K, 5,000,000 step = 10 ns
  - Production: 310 K, 4 x 12,500,000 step = 100 ns
  - Langevin / Langevin piston ayarlari korunmustur.

YEREL WSL ICIN DEGISENLER
  - Slurm/TRUBA .batch dosyalari kullanilmaz.
  - NAMD komutu: namd3 +p8 +setcpuaffinity +devices 0 <config>
  - GPUresident on eklendi.
  - dcdfreq/restartfreq/xstFreq: 1000 -> 5000 step.
    Bu dinamikleri degistirmez; yalniz dosya yazma sikligini azaltir.
    2 fs timestep ile trajectory frame araligi 10 ps olur.
  - outputEnergies/outputPressure = 500 korunmustur.
  - outputTiming = 5000 yapilmistir.
  - Tamamlanan stageler yeniden kosulmaz; zincir stage bazinda devam eder.
  - Yarim kalmis bir stage varsa dosyalari failed_runs/ altina arsivlenir ve
    stage predecessor final state'ten temizden yeniden baslatilir.

KURULUM
  mkdir -p ~/NAMD/systems/MAPK
  cd ~/NAMD/systems/MAPK
  unzip MAPK_NVT_NPT_GPU.zip
  cd NVT_NPT_GPU
  chmod +x *.sh
  bash verify_local.sh

BASLAT
  bash start_mapk.sh

DURUM
  bash status_local.sh

CANLI IZLEME
  bash watch_mapk.sh
  Cikmak icin Ctrl+C. Tek bir watch acman yeterli.

DURDUR
  bash stop_mapk.sh

DEVAM / YENIDEN BASLAT
  bash start_mapk.sh
  Tamamlanmis stageler SKIP/PASS olur. Yarim stage temizden tekrar kosar.

PERFORMANS
  Varsayilan PE=8, ADK yerel GPU akisi ile ayni tercih.
  Degistirmek gerekirse:
    PE=4 bash start_mapk.sh
  GPU secmek gerekirse:
    GPU_DEVICE=0 bash start_mapk.sh

NOTLAR
  1) NVT/NPT 298 K, production 310 K: kaynak MAPK protokolunden aynen korunmustur.
  2) PME 256^3: kaynak paketten aynen korunmustur; performans pahali olabilir.
     Bilimsel/numerik protokolu degistirmemek icin otomatik dusurulmemistir.
  3) 181,249 atom ve dcdfreq=5000 ile 100 ns production toplam DCD boyutu
     kaba olarak ~20-25 GB mertebesinde beklenir (dosya formatina gore degisir).
  4) Local zincir 6 stage: NVT -> NPT -> 0-25 -> 25-50 -> 50-75 -> 75-100 ns.
