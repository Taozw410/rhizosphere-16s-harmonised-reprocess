$root = 'D:\CodexWork\mic-re'
$log = "$root\out50k\probe_retry.log"
"retry start $(Get-Date -Format 'HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'

$man = @()
$man += Import-Csv "$root\meta\batch6_manifest.csv" | Select-Object project, crop, run, reads, fastq
$man += Import-Csv "$root\meta\batch7_maize_manifest.csv" | Select-Object project, crop, run, reads, fastq

foreach ($pool in @(@{dir='probe5'; base='probe5'}, @{dir='probe7'; base='probe7'})) {
  $base = "$root\data\$($pool.base)"
  $bad = Get-ChildItem $base -Directory | Where-Object {
    $fs = Get-ChildItem $_.FullName -File -ErrorAction SilentlyContinue
    ($fs.Count -eq 0) -or (($fs | Where-Object { $_.Length -lt 1000 }).Count -gt 0)
  }
  Add-Content $log "pool $($pool.base): $($bad.Count) project(s) need retry"
  foreach ($d in $bad) {
    $proj = $d.Name
    Add-Content $log "  $proj"
    # clear placeholders
    Get-ChildItem $d.FullName -File | Where-Object { $_.Length -lt 1000 } | Remove-Item -Force -ErrorAction SilentlyContinue
    $runs = @($man | Where-Object { $_.project -eq $proj } | Sort-Object { -[int]$_.reads })
    $ok = $false
    foreach ($r in $runs) {
      if ($ok) { break }
      $urls = @($r.fastq -split ';') | ForEach-Object { $_.Trim() } | Where-Object { $_ }
      $got = 0
      foreach ($u in $urls) {
        $fn = Split-Path $u -Leaf
        $out = Join-Path $d.FullName $fn
        & curl.exe -L --ssl-no-revoke --fail --silent --show-error --retry 5 --retry-delay 5 --max-time 900 -C - -o $out "https://$u" 2>> $log
        if ((Test-Path $out) -and ((Get-Item $out).Length -gt 1000)) { $got++ }
      }
      if ($got -eq $urls.Count -and $got -gt 0) { $ok = $true; Add-Content $log "    ok via $($r.run) ($got files)" }
      else {
        Get-ChildItem $d.FullName -File | Where-Object { $_.Length -lt 1000 } | Remove-Item -Force -ErrorAction SilentlyContinue
        Add-Content $log "    failed via $($r.run)"
      }
    }
    if (-not $ok) { Add-Content $log "    GIVE UP $proj" }
  }
}
Add-Content $log "retry done $(Get-Date -Format 'HH:mm:ss')"
