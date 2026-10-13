# 等汇总生成后，自动重跑三分类清单
$log = 'D:\CodexWork\mic-re\out50k\batch2.log'
while ($true) {
  Start-Sleep -Seconds 120
  if (Test-Path $log) {
    $c = Get-Content $log -Raw -ErrorAction SilentlyContinue
    if ($c -match 'summary generated') { break }
  }
}
& 'D:\bioenv\envs\riz16s\python.exe' 'D:\CodexWork\mic-re\scripts\summarize_w1.py' *>> 'D:\CodexWork\mic-re\out50k\summary.log'
& 'D:\bioenv\envs\riz16s\python.exe' 'D:\CodexWork\mic-re\scripts\classify_w1.py' *>> 'D:\CodexWork\mic-re\out50k\summary.log'
Add-Content -LiteralPath $log -Value "tiers regenerated $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
