$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\batch7_fetch.log"
"batch7 fetch start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
$py = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'

while ($true) {
  Start-Sleep -Seconds 120
  if ((Test-Path "$root\out50k\batch7_prep.log") -and ((Get-Content "$root\out50k\batch7_prep.log" -Raw) -match 'PREP DONE')) { break }
}
$projs = Get-Content "$root\meta\batch7_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() }
Add-Content $log "selected projects = $($projs.Count): $($projs -join ',')"

& "$root\scripts\download_ena_parallel.ps1" `
  -ManifestPath "$root\meta\batch7_download_manifest.csv" `
  -RawRoot "$root\data\raw" -Projects $projs -PerProject 30 -Throttle 12 -MaxTime 2400 `
  -LogPath "$root\out50k\batch7_dl.log" *>> $log
Add-Content $log "download done $(Get-Date -Format 'HH:mm:ss')"

& $py "$root\scripts\verify_gzip_raw.py" "$root\data\raw" ($projs -join ',') --workers 2 --out "$root\out50k\gzip_check_batch7.csv" *>> $log
$gz = $LASTEXITCODE
& $py "$root\scripts\subsample_fastq.py" "$root\data\raw" "$root\data\raw50k" 50000 20261008 ($projs -join ',') *>> $log
$sub = $LASTEXITCODE
if ($gz -eq 0 -and $sub -eq 0) {
  Add-Content $log "FETCH DONE $(Get-Date -Format 'HH:mm:ss')"
} else {
  Add-Content $log "FETCH FAILED (gzip_check=$gz subsample=$sub) $(Get-Date -Format 'HH:mm:ss')"
}
