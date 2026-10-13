$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\longread_chain.log"
"longread chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# 1) 等两条短读流结束
while ($true) {
  Start-Sleep -Seconds 120
  $a = (Test-Path "$root\out50k\streamA.log") -and ((Get-Content "$root\out50k\streamA.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  $b = (Test-Path "$root\out50k\streamB.log") -and ((Get-Content "$root\out50k\streamB.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  if ($a -and $b) { break }
}
Add-Content $log "short-read streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 2) 长读项目按样本量均衡分两组
$projs = Get-ChildItem "$root\data\raw50k_longread" -Directory | Select-Object -ExpandProperty Name
$info = foreach ($p in $projs) {
  $f = @(Get-ChildItem "$root\data\raw50k_longread\$p" -File)
  $paired = @($f | Where-Object { $_.Name -match '_1\.fastq\.gz$' }).Count
  $single = @($f | Where-Object { $_.Name -notmatch '_[12]\.fastq\.gz$' }).Count
  [pscustomobject]@{ project = $p; n = $(if ($paired -gt 0) { $paired } else { $single }) }
}
$gA=@(); $gB=@(); $sA=0; $sB=0
foreach ($x in ($info | Sort-Object -Property n -Descending)) {
  if ($sA -le $sB) { $gA += $x.project; $sA += $x.n } else { $gB += $x.project; $sB += $x.n }
}
$A = ($gA -join ','); $B = ($gB -join ',')
Add-Content $log "split: A=$A ($sA)  B=$B ($sB)" -Encoding utf8

# 3) 双流跑长读
$pA = Start-Process -FilePath 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe' -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',"conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$A' '$root\data\raw50k_longread' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\longreadA.log'")
Start-Sleep -Seconds 10
$pB = Start-Process -FilePath 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe' -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',"conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$B' '$root\data\raw50k_longread' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\longreadB.log'")
Add-Content $log "launched longread A=$($pA.Id) B=$($pB.Id) $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pA.WaitForExit(); Add-Content $log "longread A done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
$pB.WaitForExit(); Add-Content $log "longread B done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

# 4) 报告
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
Add-Content $log "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
