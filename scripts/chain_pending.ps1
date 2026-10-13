$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\pend_chain.log"
"pending chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while ($true) {
  Start-Sleep -Seconds 120
  $a = (Test-Path "$root\out50k\pend_streamA.log") -and ((Get-Content "$root\out50k\pend_streamA.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  $b = (Test-Path "$root\out50k\pend_streamB.log") -and ((Get-Content "$root\out50k\pend_streamB.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  if ($a -and $b) { break }
}
Add-Content $log "short done $(Get-Date -Format 'HH:mm:ss'); starting long-read streams" -Encoding utf8
$projs = Get-ChildItem "$root\data\raw50k_longread" -Directory | Select-Object -ExpandProperty Name
$info = foreach ($p in $projs) {
  $f = @(Get-ChildItem "$root\data\raw50k_longread\$p" -File)
  $paired = @($f | Where-Object { $_.Name -match '_1\.fastq\.gz$' }).Count
  $single = @($f | Where-Object { $_.Name -notmatch '_[12]\.fastq\.gz$' }).Count
  [pscustomobject]@{ project=$p; n=$(if ($paired -gt 0) { $paired } else { $single }) }
}
$A=@();$B=@();$sA=0;$sB=0
foreach ($x in ($info | Sort-Object -Property n -Descending)) {
  if ($sA -le $sB) { $A += $x.project; $sA += $x.n } else { $B += $x.project; $sB += $x.n }
}
Add-Content $log "long split A=$($A -join ',') ($sA)  B=$($B -join ',') ($sB)" -Encoding utf8
$pwsh='C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$pA = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',"conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($A -join ',')' '$root\data\raw50k_longread' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\pend_longA.log'")
Start-Sleep -Seconds 10
$pB = Start-Process -FilePath $pwsh -WindowStyle Hidden -PassThru -ArgumentList @('-NoProfile','-Command',"conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\run_dada2_parallel.R' '$($B -join ',')' '$root\data\raw50k_longread' '$root\out_stage2' 'D:/bioenv/ref/silva_nr99_v138.1_train_set.fa.gz' 2 1e7 6 '$root\meta\primer_table_stage2.csv' *>> '$root\out50k\pend_longB.log'")
$pA.WaitForExit(); $pB.WaitForExit()
Add-Content $log "long done $(Get-Date -Format 'HH:mm:ss'); gate + reports" -Encoding utf8
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\master_scan.py" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\merged_report.py" *>> $log
Add-Content $log "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
