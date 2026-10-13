# Repair chain for batch7/batch8 data: wait for the re-download, verify gzip integrity (with up to
# 3 bounded re-download rounds), re-subsample only the incomplete projects, then check file counts.
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out50k\repair_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
"repair chain start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

while ($true) {
  Start-Sleep -Seconds 120
  if ((Test-Path "$root\out50k\repair_dl.log") -and ((Get-Content "$root\out50k\repair_dl.log" -Raw -ErrorAction SilentlyContinue) -match 'REPAIR DOWNLOAD DONE')) { break }
}
Add-Content -LiteralPath $log -Value "first repair download finished $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8

$b7 = @(Get-Content "$root\meta\batch7_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() })
$b8 = @(Get-Content "$root\meta\batch8_selected.txt" -Encoding UTF8 | Where-Object { $_.Trim() } | ForEach-Object { $_.Trim() })
$all = @($b7 + $b8)
$csv = "$root\out50k\gzip_check_repair.csv"

$ok = $false
for ($round = 1; $round -le 3; $round++) {
  & $py "$root\scripts\verify_gzip_raw.py" "$root\data\raw" ($all -join ',') --workers 2 --out $csv *>> $log
  $gz = $LASTEXITCODE
  Add-Content -LiteralPath $log -Value "round $round gzip check exit=$gz $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
  if ($gz -eq 0) { $ok = $true; break }

  $bad = @(Import-Csv $csv | Where-Object { $_.status -ne 'OK' })
  $badProjects = @($bad | Select-Object -ExpandProperty project -Unique)
  Add-Content -LiteralPath $log -Value "round $round bad files=$($bad.Count) projects=$($badProjects -join ',')" -Encoding utf8
  $dest = "$root\_recycle\corrupt_fastq_round$round"
  foreach ($b in $bad) {
    $dd = Join-Path $dest $b.project
    New-Item -ItemType Directory -Force -Path $dd | Out-Null
    $srcf = "$root\data\raw\$($b.project)\$($b.file)"
    if (Test-Path $srcf) { Move-Item -LiteralPath $srcf -Destination (Join-Path $dd $b.file) -Force }
  }
  $p7 = @($badProjects | Where-Object { $b7 -contains $_ })
  $p8 = @($badProjects | Where-Object { $b8 -contains $_ })
  if ($p7.Count -gt 0) {
    & "$root\scripts\download_ena_parallel.ps1" -ManifestPath "$root\meta\batch7_download_manifest.csv" -RawRoot "$root\data\raw" -Projects $p7 -PerProject 30 -Throttle 4 -MaxTime 14400 -LogPath "$root\out50k\repair_batch7_dl.log" *>> $log
  }
  if ($p8.Count -gt 0) {
    & "$root\scripts\download_ena_parallel.ps1" -ManifestPath "$root\meta\batch8_maize_manifest.csv" -RawRoot "$root\data\raw" -Projects $p8 -PerProject 30 -Throttle 4 -MaxTime 14400 -LogPath "$root\out50k\repair_batch8_dl.log" *>> $log
  }
}
if (-not $ok) {
  Add-Content -LiteralPath $log -Value "REPAIR FAILED: gzip still broken after retries $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
  exit 1
}

$b7todo = @($b7 | Where-Object { $_ -ne 'PRJDB19884' })
& $py "$root\scripts\subsample_fastq.py" "$root\data\raw" "$root\data\raw50k" 50000 20261008 ($b7todo -join ',') *>> $log
$s1 = $LASTEXITCODE
& $py "$root\scripts\subsample_fastq.py" "$root\data\raw" "$root\data\raw50k" 50000 20261008 'PRJNA427813,PRJNA593708' *>> $log
$s2 = $LASTEXITCODE
Add-Content -LiteralPath $log -Value "subsample exit codes: batch7=$s1 batch8=$s2" -Encoding utf8

& $py "$root\scripts\check_subsample_counts.py" "$root\data\raw" "$root\data\raw50k" ($all -join ',') --out "$root\out50k\subsample_counts_repair.csv" *>> $log
$cc = $LASTEXITCODE
if ($s1 -eq 0 -and $s2 -eq 0 -and $cc -eq 0) {
  Add-Content -LiteralPath $log -Value "REPAIR ALL DONE $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
} else {
  Add-Content -LiteralPath $log -Value "REPAIR INCOMPLETE (subsample=$s1,$s2 counts=$cc) $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
}
