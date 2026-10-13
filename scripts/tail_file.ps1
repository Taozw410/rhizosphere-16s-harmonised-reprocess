param([Parameter(Mandatory=$true)][string]$Path, [int]$Lines = 25)
if (-not (Test-Path -LiteralPath $Path)) { "MISSING $Path"; exit 0 }
$fs = [System.IO.File]::Open($Path, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::ReadWrite)
try {
  $len = $fs.Length
  $take = [Math]::Min($len, 65536)
  $fs.Seek(-$take, [System.IO.SeekOrigin]::End) | Out-Null
  $buf = New-Object byte[] $take
  $read = $fs.Read($buf, 0, $take)
  $txt = [System.Text.Encoding]::UTF8.GetString($buf, 0, $read)
} finally { $fs.Close() }
$arr = $txt -split "`r?`n"
$n = [Math]::Min($Lines, $arr.Count)
$arr[($arr.Count - $n)..($arr.Count - 1)] -join "`n"
