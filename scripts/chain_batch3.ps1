$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\batch3_chain.log"
"batch3 chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
# 等两条流结束
while ($true) {
  Start-Sleep -Seconds 120
  $a = (Test-Path "$root\out50k\batch3_streamA.log") -and ((Get-Content "$root\out50k\batch3_streamA.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  $b = (Test-Path "$root\out50k\batch3_streamB.log") -and ((Get-Content "$root\out50k\batch3_streamB.log" -Raw -ErrorAction SilentlyContinue) -match 'all done')
  if ($a -and $b) { break }
}
Add-Content $log "streams done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
# 标记基因闸门 + 主清单
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\master_scan.py" *>> $log
Add-Content $log "master scan done $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
# 汇总与报告
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\summarize_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\classify_w1.py" "$root\out_stage2" *>> $log
& 'D:\bioenv\envs\riz16s\python.exe' "$root\scripts\merged_report.py" *>> $log
Add-Content $log "ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
