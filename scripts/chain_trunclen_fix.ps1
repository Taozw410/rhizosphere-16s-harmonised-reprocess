# truncLen fix chain: run AFTER batch8 chain finishes.
# 1) wait for batch8 "ALL DONE" and for all R pipeline processes to exit
# 2) swap in the patched pipeline (parameter-fingerprint guard), keep a backup
# 3) archive the superseded out_stage2/<proj> dirs (move, never delete)
# 4) set measured truncLen rows in meta/primer_table_stage2.csv
# 5) rerun those projects in two streams (2 cores each)
# 6) refresh inventory / marker gate / summaries / merged report, then write the A/B table
$root  = 'D:\CodexWork\mic-re'
$log   = "$root\out50k\trunclen_chain.log"
$py    = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh  = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref   = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
$projs = @('PRJNA596848', 'PRJNA665860', 'PRJNA996548')
$stamp = '20261008'

"trunclen fix chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# ---- 1) wait for the whole batch8 chain, then for idle R ----
while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\batch8_chain.log") -and ((Get-Content "$root\out50k\batch8_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'ALL DONE')
  if (-not $done) { continue }
  $r = @(Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" -ErrorAction SilentlyContinue |
         Where-Object { $_.CommandLine -like '*run_dada2_parallel.R*' })
  if ($r.Count -eq 0) { break }
  Add-Content -LiteralPath $log -Value "waiting: $($r.Count) R pipeline process(es) still running $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
}
Add-Content -LiteralPath $log -Value "batch8 finished and R idle $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- 2) swap in the patched pipeline ----
$src = "$root\scripts\run_dada2_parallel.R"
$new = "$root\scripts\_bak\run_dada2_parallel.patched.R"
if ((Test-Path $new) -and ((Get-FileHash -LiteralPath $src).Hash -ne (Get-FileHash -LiteralPath $new).Hash)) {
  try {
    Copy-Item -LiteralPath $src -Destination "$root\scripts\_bak\run_dada2_parallel.R.prepatch_$stamp.bak" -Force -ErrorAction Stop
    Copy-Item -LiteralPath $new -Destination $src -Force -ErrorAction Stop
    conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript "$root\scripts\_bak\parse_check.R" *>> $log
    Add-Content -LiteralPath $log -Value "pipeline patched with parameter-fingerprint guard" -Encoding utf8
  } catch {
    Add-Content -LiteralPath $log -Value "pipeline patch SKIPPED: $($_.Exception.Message)" -Encoding utf8
  }
}

# ---- 3) archive superseded results ----
foreach ($p in $projs) {
  $d = "$root\out_stage2\$p"
  if (Test-Path $d) {
    $beforeRoot = "$root\out50k\_archive_trunclen_before"
    New-Item -ItemType Directory -Force -Path $beforeRoot | Out-Null
    $arch = Join-Path $beforeRoot $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath $d -Destination $arch -Force
    Add-Content -LiteralPath $log -Value "archived $p -> $arch" -Encoding utf8
  } else {
    Add-Content -LiteralPath $log -Value "$p has no previous result dir (first run)" -Encoding utf8
  }
}

# ---- 4) measured truncLen rows (read-through: L=292 vs R=301 -> 250/250; long insert with
#         3' tail in the overlap region: L=469 vs R=294 -> 280/250, overlap 61 bp) ----
& $py "$root\scripts\set_trunclen_rows.py" "$root\meta\primer_table_stage2.csv" `
  "PRJNA996548,250,250,trunclen_L292_R301_readthrough" `
  "PRJNA596848,280,250,trunclen_L469_R294_tail_quality" `
  "PRJNA665860,280,250,trunclen_L469_R294_tail_quality" *>> $log

# ---- 5) rerun in two streams ----
$A = 'PRJNA596848,PRJNA665860'
$B = 'PRJNA996548'
$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile', '-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$A' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\trunclen_streamA.log'")
Start-Sleep -Seconds 15
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile', '-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$B' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\trunclen_streamB.log'")
Add-Content -LiteralPath $log -Value "streams launched $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content -LiteralPath $log -Value "streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# ---- 6) gates + merged report + A/B table ----
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$root\out50k\_archive_trunclen_before" "$root\out_stage2" ($projs -join ',') --out "$root\out50k\trunclen_ab_report.md" *>> $log
Add-Content -LiteralPath $log -Value "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

