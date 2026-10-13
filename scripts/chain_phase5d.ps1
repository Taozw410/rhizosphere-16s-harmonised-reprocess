# Phase 5d：PRJNA596848 / PRJNA665860 的"短读长样本组"修正变体（truncLen 244,244）
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\phase5d_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$ref  = 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz'
"phase5d start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while ($true) {
  Start-Sleep -Seconds 300
  $done = (Test-Path "$root\out50k\phase5c_chain.log") -and ((Get-Content "$root\out50k\phase5c_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'PHASE5C ALL DONE')
  if ($done) { break }
}
Add-Content -LiteralPath $log -Value "phase5c finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$before = "$root\out50k\_archive_phase5d_before"
New-Item -ItemType Directory -Force -Path $before | Out-Null
foreach ($p in @('PRJNA596848','PRJNA665860')) {
  if (Test-Path "$root\out_stage2\$p") {
    $arch = Join-Path $before $p
    if (Test-Path $arch) { $arch = "${arch}_$(Get-Date -Format 'HHmmss')" }
    Move-Item -LiteralPath "$root\out_stage2\$p" -Destination $arch -Force
  }
}
$p1 = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel_keepall.R' 'PRJNA596848' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_phase5d.csv' *>> '$root\out50k\f5d_streamA.log'")
Start-Sleep -Seconds 15
$p2 = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel_keepall.R' 'PRJNA665860' '$root\data\raw50k' '$root\out_stage2' '$ref' 2 1e7 6 '$root\meta\primer_table_phase5d.csv' *>> '$root\out50k\f5d_streamB.log'")
Add-Content -LiteralPath $log -Value "streams launched $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$p1.WaitForExit(); $p2.WaitForExit()
Add-Content -LiteralPath $log -Value "streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
& $py "$root\scripts\master_scan.py" *>> $log
& $py "$root\scripts\marker_gene_qc.py" *>> $log
& $py "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& $py "$root\scripts\merged_report.py" *>> $log
& $py "$root\scripts\pick_phase5b_winner.py" *>> $log
& $py "$root\scripts\trunclen_ab_report.py" "$root\out50k\_archive_phase5b_before" "$root\out_stage2" 'PRJNA596848,PRJNA665860' --out "$root\out50k\phase5d_ab_report.md" *>> $log
& $py "$root\scripts\select_final45.py" --out "$root\meta\final45_selection.csv" --report "$root\out50k\选集报告.md" *>> $log
Add-Content -LiteralPath $log -Value "PHASE5D ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
