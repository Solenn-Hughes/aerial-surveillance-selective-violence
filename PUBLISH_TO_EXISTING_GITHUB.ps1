param([Parameter(Mandatory=$true)][string]$RepoUrl)
$ErrorActionPreference = 'Stop'
if (-not (Get-Command git -ErrorAction SilentlyContinue)) { throw 'Git is not installed. Use the GitHub website Upload files option instead.' }
$src = $PSScriptRoot
if (-not (Test-Path (Join-Path $src 'RUN_VERIFIED_14.R'))) { throw 'Run this script from the extracted ZIP folder.' }
$dest = Join-Path $env:TEMP 'jpr-replication-github-upload'
if (Test-Path $dest) { throw "Temporary target folder already exists: $dest. Remove/rename it first to avoid overwriting work." }
& git clone -- $RepoUrl $dest
if ($LASTEXITCODE -ne 0) { throw 'Git clone failed. Check URL and sign in.' }
$items = @('README.md','RUN_VERIFIED_14.R','PREPARE_LOCAL_DATA.R','RUN_RAW_UPSTREAM.R','RUN_APPENDIX.R','RUN_ALL.R','.gitignore','code','documentation','reference','GITHUB_UPLOAD_INSTRUCTIONS.txt')
foreach ($i in $items) { Copy-Item -LiteralPath (Join-Path $src $i) -Destination $dest -Recurse -Force }
Push-Location $dest
try {
  & git add README.md RUN_VERIFIED_14.R PREPARE_LOCAL_DATA.R RUN_RAW_UPSTREAM.R RUN_APPENDIX.R RUN_ALL.R .gitignore code documentation reference GITHUB_UPLOAD_INSTRUCTIONS.txt
  if ($LASTEXITCODE -ne 0) { throw 'git add failed' }
  & git status --short
  Write-Host 'Review the files above. Confirm your repository is PRIVATE and no restricted data are staged.'
  $ok = Read-Host 'Type PUSH to commit and push'
  if ($ok -cne 'PUSH') { Write-Host 'Cancelled. Nothing pushed.'; exit 0 }
  & git commit -m 'Add verified model replication code and provenance notes'
  if ($LASTEXITCODE -ne 0) { throw 'git commit failed; configure Git username/email if needed.' }
  & git push
  if ($LASTEXITCODE -ne 0) { throw 'git push failed. Check authentication/branch permissions.' }
  Write-Host 'Pushed to GitHub.'
} finally { Pop-Location }
