# 等第二批跑完，自动生成汇总
$log = 'D:\CodexWork\mic-re\out50k\batch2.log'
while ($true) {
  Start-Sleep -Seconds 120
  if (Test-Path $log) {
    $c = Get-Content $log -Raw -ErrorAction SilentlyContinue
    if ($c -match 'batch2 done') { break }
  }
}
& 'D:\bioenv\envs\riz16s\python.exe' 'D:\CodexWork\mic-re\scripts\summarize_w1.py' *>> 'D:\CodexWork\mic-re\out50k\summary.log'
Add-Content -LiteralPath $log -Value "summary generated $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
