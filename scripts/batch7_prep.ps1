$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\batch7_prep.log"
"batch7 prep start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
$py = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'

# 1) wait until no curl is left and both probe pools are populated
$stable = 0
while ($true) {
  Start-Sleep -Seconds 60
  $nc = @(Get-Process curl -ErrorAction SilentlyContinue).Count
  $n5 = @(Get-ChildItem "$root\data\probe5" -Directory -ErrorAction SilentlyContinue).Count
  $n7 = @(Get-ChildItem "$root\data\probe7" -Directory -ErrorAction SilentlyContinue).Count
  if ($nc -eq 0 -and $n5 -ge 106 -and $n7 -ge 10) { $stable++ } else { $stable = 0 }
  Add-Content $log "wait: curl=$nc probe5=$n5 probe7=$n7 stable=$stable $(Get-Date -Format 'HH:mm:ss')"
  if ($stable -ge 2) { break }
}
Add-Content $log "downloads settled $(Get-Date -Format 'HH:mm:ss')"

# 2) retry stragglers
& $pwsh -NoProfile -ExecutionPolicy Bypass -File "$root\scripts\retry_probes.ps1" *>> $log
Add-Content $log "retry pass done $(Get-Date -Format 'HH:mm:ss')"

# 3) strict screening
& $py "$root\scripts\screen_batch7.py" *>> $log
Add-Content $log "screening done $(Get-Date -Format 'HH:mm:ss')"

# 4) final selection + primer table
& $py "$root\scripts\select_final_b7.py" *>> $log
& $py "$root\scripts\make_primer_table_batch7.py" *>> $log
Add-Content $log "PREP DONE $(Get-Date -Format 'HH:mm:ss')"
