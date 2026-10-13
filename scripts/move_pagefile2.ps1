$log = 'D:\CodexWork\mic-re\out50k\pagefile_move.log'
Add-Content -LiteralPath $log -Value "retry with uint32 $(Get-Date -Format 'HH:mm:ss')" -Encoding utf8
try {
    $pf = New-CimInstance -ClassName Win32_PageFileSetting -Property @{
        Name        = 'D:\pagefile.sys'
        InitialSize = [uint32]8192
        MaximumSize = [uint32]24576
    } -ErrorAction Stop
    Add-Content $log "created: $($pf.Name) init=$($pf.InitialSize) max=$($pf.MaximumSize)" -Encoding utf8
} catch {
    Add-Content $log "New-CimInstance failed: $($_.Exception.Message)" -Encoding utf8
    # 退回：直接写注册表
    try {
        Set-ItemProperty -Path 'HKLM:\SYSTEM\CurrentControlSet\Control\Session Manager\Memory Management' `
            -Name PagingFiles -Value @('D:\pagefile.sys 8192 24576') -Type MultiString -ErrorAction Stop
        Add-Content $log "registry fallback written: D:\\pagefile.sys 8192 24576" -Encoding utf8
    } catch {
        Add-Content $log "registry fallback failed: $($_.Exception.Message)" -Encoding utf8
    }
}
Add-Content $log "--- verify ---" -Encoding utf8
Get-CimInstance Win32_PageFileSetting -ErrorAction SilentlyContinue | ForEach-Object { Add-Content $log "  WMI: $($_.Name) init=$($_.InitialSize) max=$($_.MaximumSize)" -Encoding utf8 }
(Get-ItemProperty 'HKLM:\SYSTEM\CurrentControlSet\Control\Session Manager\Memory Management' -Name PagingFiles -ErrorAction SilentlyContinue).PagingFiles | ForEach-Object { Add-Content $log "  REG: $_" -Encoding utf8 }
$c2 = Get-CimInstance Win32_ComputerSystem
Add-Content $log "AutomaticManagedPagefile = $($c2.AutomaticManagedPagefile)" -Encoding utf8
Add-Content $log "DONE" -Encoding utf8
