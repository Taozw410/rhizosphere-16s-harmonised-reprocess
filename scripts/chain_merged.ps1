$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\merged_report_chain.log"
"merged-report watcher start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
while ($true) {
  Start-Sleep -Seconds 120
  $lr = "$root\out50k\longread_chain.log"
  if ((Test-Path $lr) -and ((Get-Content $lr -Raw -ErrorAction SilentlyContinue) -match 'ALL DONE')) { break }
}
Add-Content $log "longread chain done $(Get-Date -Format 'HH:mm:ss'); generating merged report" -Encoding utf8
# 再刷一次两阶段的三分类，然后出合并报告
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out50k" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out50k" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\merged_report.py" *>> $log
Add-Content $log "merged report done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
