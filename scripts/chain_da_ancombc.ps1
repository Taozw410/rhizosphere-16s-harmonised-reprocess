# 等代理恢复 -> 补装 microbiome（ANCOMBC2 的输入转换依赖）-> 跑完 ANCOMBC2 -> 写 da_report.md
$root = 'D:\CodexWork\mic-re'
$log  = "$root\out_analysis\da_ancombc_chain.log"
$py   = 'D:\bioenv\envs\riz16s\python.exe'
$pwsh = 'C:\Users\hthome\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\powershell\pwsh.exe'
$px   = 'http://127.0.0.1:7888'
function Say { param([string]$m) Add-Content -LiteralPath $log -Value "$((Get-Date -Format 'HH:mm:ss'))  $m" -Encoding utf8 }
"da ancombc chain start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Set-Content -LiteralPath $log -Encoding utf8

# 1) 等代理
$env:https_proxy = $px; $env:http_proxy = $px
$up = $false
for ($i = 1; $i -le 240; $i++) {
  & curl.exe -L --ssl-no-revoke --silent --show-error --max-time 12 -o NUL "https://bioconductor.org/" 2>$null
  if ($LASTEXITCODE -eq 0) { $up = $true; Say "proxy is up (attempt $i)"; break }
  Start-Sleep -Seconds 60
}
if (-not $up) { Say 'proxy never came back within 4 h -> ABORT'; exit 1 }

# 2) 计算 microbiome 的依赖闭包
& $pwsh -NoProfile -Command "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\install_plan3.R'" *>> $log
if ($LASTEXITCODE -ne 0) { Say 'install_plan3 failed'; exit 2 }

# 3) 用 curl 抓 zip（重试）
$rows = Import-Csv -LiteralPath "$root\out_analysis\install_plan3.csv"
$dir  = "$root\out_analysis\pkgzips"
New-Item -ItemType Directory -Force -Path $dir | Out-Null
Say "planned packages: $($rows.Count)"
foreach ($r in $rows) {
  $f = Join-Path $dir ($r.pkg + '_' + $r.ver + '.zip')
  if ((Test-Path $f) -and ((Get-Item $f).Length -gt 20000)) { continue }
  $u = $r.url.TrimEnd('/') + '/' + $r.pkg + '_' + $r.ver + '.zip'
  for ($k = 1; $k -le 12; $k++) {
    & curl.exe -L --ssl-no-revoke --silent --show-error --connect-timeout 30 --max-time 900 -C - -o $f $u 2>$null
    if ((Test-Path $f) -and ((Get-Item $f).Length -gt 20000)) { break }
    Start-Sleep -Seconds 6
  }
  Say ("  {0}: {1} KB" -f $r.pkg, $(if (Test-Path $f) { [math]::Round((Get-Item $f).Length/1KB,0) } else { 0 }))
}

# 4) 本地安装
& $pwsh -NoProfile -Command "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\install_local_zips.R'" *>> $log
Say 'local install attempted'

# 5) 跑 ANCOMBC2（da_aldex_ancombc.R 会复用已缓存的 da_aldex2.rds，只补 ANCOMBC2 与报告）
& $pwsh -NoProfile -Command "conda run --no-capture-output -p 'D:\bioenv\envs\riz16s' Rscript '$root\scripts\da_aldex_ancombc.R'" *>> "$root\out_analysis\da.log"
Say "da script exit=$LASTEXITCODE"
Add-Content -LiteralPath $log -Value "$(Get-Date -Format 'HH:mm:ss')  DA ANCOMBC DONE" -Encoding utf8
