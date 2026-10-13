# Batch 7 (new 20 projects): wait for the pending batch, then process with two streams, then gate + reports.
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\batch7_chain.log"
"batch7 chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# 1) wait for the pending batch chain (16 projects + long reads + gates + reports)
while ($true) {
  Start-Sleep -Seconds 180
  $done = (Test-Path "$root\out50k\pend_chain.log") -and ((Get-Content "$root\out50k\pend_chain.log" -Raw -ErrorAction SilentlyContinue) -match 'ALL DONE')
  if ($done) { break }
}
Add-Content $log "pending batch finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 2) primer table
if (-not (Test-Path "$root\meta\primer_table_batch7.csv")) {
  & 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\make_primer_table_batch7.py" *>> $log
}
Add-Content $log "primer table ready $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 3) balance two streams by sample count
$projs = Get-Content "$root\meta\batch7_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() }
$info = foreach ($p in $projs) {
  $f = @(Get-ChildItem "$root\data\raw50k\$p" -File -ErrorAction SilentlyContinue)
  [pscustomobject]@{ project = $p; n = $f.Count }
}
$A = @(); $B = @(); $sA = 0; $sB = 0
foreach ($x in ($info | Sort-Object -Property n -Descending)) {
  if ($sA -le $sB) { $A += $x.project; $sA += $x.n } else { $B += $x.project; $sB += $x.n }
}
Add-Content $log "b7 stream A=$($A -join ',') ($sA samples)" -Encoding utf8
Add-Content $log "b7 stream B=$($B -join ',') ($sB samples)" -Encoding utf8

$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile', '-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($A -join ',')' '$root\data\raw50k' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_batch7.csv' *>> '$root\out50k\batch7_streamA.log'")
Start-Sleep -Seconds 15
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile', '-Command',
  "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($B -join ',')' '$root\data\raw50k' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_batch7.csv' *>> '$root\out50k\batch7_streamB.log'")
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content $log "b7 streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 4) inventory + marker-gene gate + summaries
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\master_scan.py" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\marker_gene_qc.py" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\merged_report.py" *>> $log
Add-Content $log "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
