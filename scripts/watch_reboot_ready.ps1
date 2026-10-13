$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\reboot_ready_watcher.log"
$flag = "$root\REBOOT_READY.txt"
Remove-Item $flag -Force -ErrorAction SilentlyContinue
"watcher start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while ($true) {
  Start-Sleep -Seconds 120
  $a = (Test-Path "$root\out50k\streamA.log") -and ((Get-Content "$root\out50k\streamA.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  $b = (Test-Path "$root\out50k\streamB.log") -and ((Get-Content "$root\out50k\streamB.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  if ($a -and $b) { break }
}
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\merged_report.py" *>> $log
$done = @(Get-ChildItem "$root\out_stage2" -Directory | Where-Object { Test-Path (Join-Path $_.FullName 'asv_taxonomy.csv') }).Count
"SHORT-READ BATCH DONE at $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss'); projects=$done; READY FOR REBOOT" | Set-Content -LiteralPath $flag -Encoding utf8
Add-Content $log "READY FOR REBOOT $(Get-Date -Format 'HH:mm:ss') done=$done" -Encoding utf8
