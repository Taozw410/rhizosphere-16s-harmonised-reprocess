param(
  [Parameter(Mandatory=$true)][string]$ManifestPath,
  [Parameter(Mandatory=$true)][string]$RawRoot,
  [string[]]$Projects = @(),
  [int]$PerProject = 30,
  [int]$Throttle = 12,
  [int]$MaxTime = 1800,
  [string]$LogPath = ""
)

$man = Import-Csv $ManifestPath
if ($Projects.Count -gt 0) { $man = $man | Where-Object { $Projects -contains $_.project } }

function Log([string]$m) {
  $line = "{0}  {1}" -f (Get-Date -Format 'HH:mm:ss'), $m
  Write-Host $line
  if ($LogPath) { Add-Content -LiteralPath $LogPath -Value $line -Encoding utf8 }
}

$jobs = New-Object System.Collections.ArrayList
foreach ($p in ($man | Group-Object project)) {
  $runs = $p.Group | Sort-Object run
  $n = $runs.Count
  if ($n -le $PerProject) {
    $take = $runs
  } else {
    # 等距抽样：覆盖整个 accession 区间，避免只落到某一处理组
    $idx = @()
    for ($i = 0; $i -lt $PerProject; $i++) { $idx += [int][math]::Floor($i * $n / $PerProject) }
    $take = $runs[$idx]
  }
  Log ("{0}: {1} runs available -> taking {2}" -f $p.Name, $n, $take.Count)
  foreach ($t in $take) {
    foreach ($u in ($t.fastq -split ';')) {
      $u = $u.Trim()
      if ($u) { [void]$jobs.Add([pscustomobject]@{ Project = $t.project; Url = 'https://' + $u }) }
    }
  }
}
if ($jobs.Count -eq 0) { Log 'no jobs'; exit 1 }

Log ("TOTAL jobs = {0}  throttle = {1}" -f $jobs.Count, $Throttle)
$sw = [System.Diagnostics.Stopwatch]::StartNew()

$results = $jobs | ForEach-Object -Parallel {
  $raw = $using:RawRoot
  $mt  = $using:MaxTime
  $dir = Join-Path $raw $_.Project
  if (-not (Test-Path $dir)) { New-Item -ItemType Directory -Force -Path $dir | Out-Null }
  $fn  = Split-Path $_.Url -Leaf
  $out = Join-Path $dir $fn
  if ((Test-Path $out) -and ((Get-Item $out).Length -gt 1000)) {
    return [pscustomobject]@{ file = $fn; status = 'skip'; bytes = (Get-Item $out).Length }
  }
  # 注意：只看文件大小会把"被 --max-time 截断"的下载当成成功（本项目实测：11 个 fastq 因此损坏）。
  # 这里改为看 curl 退出码 + 最小体积；失败则用更大的超时续传一次。
  & curl.exe -L --ssl-no-revoke --fail --silent --show-error --retry 3 --retry-delay 5 --max-time $mt -C - -o $out $_.Url 2>$null
  $code = $LASTEXITCODE
  $sz = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }
  if ($code -ne 0 -or $sz -lt 10000) {
    & curl.exe -L --ssl-no-revoke --fail --silent --show-error --retry 3 --retry-delay 5 --max-time 14400 -C - -o $out $_.Url 2>$null
    $code = $LASTEXITCODE
    $sz = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }
  }
  if ($code -eq 0 -and $sz -gt 10000) {
    return [pscustomobject]@{ file = $fn; status = 'ok'; bytes = $sz }
  } else {
    return [pscustomobject]@{ file = $fn; status = 'FAIL'; bytes = $sz }
  }
} -ThrottleLimit $Throttle

$sw.Stop()
$ok = @($results | Where-Object { $_.status -eq 'ok' })
$sk = @($results | Where-Object { $_.status -eq 'skip' })
$bad = @($results | Where-Object { $_.status -eq 'FAIL' })
Log ("elapsed    : {0:N1} min" -f $sw.Elapsed.TotalMinutes)
Log ("downloaded : {0}   skipped: {1}   failed: {2}" -f $ok.Count, $sk.Count, $bad.Count)
Log ("bytes      : {0:N1} MB   ({1:N2} MB/s)" -f (($results|Measure-Object bytes -Sum).Sum/1MB), ((($results|Measure-Object bytes -Sum).Sum/1MB)/($sw.Elapsed.TotalSeconds)))
if ($bad.Count -gt 0) { $bad | ForEach-Object { Log ('  FAILED ' + $_.file) } }
