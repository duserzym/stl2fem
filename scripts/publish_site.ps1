# Publish docs/ (landing page + grain viewer with its data bundle) as a single-commit gh-pages branch, force-pushed.
# The data bundle is rebuilt by scripts/build_grain_viewer.py; publishing this way keeps every rebuild out of the
# repository history. GitHub Pages serves the gh-pages branch root.
param([string]$Remote = "origin", [string]$Branch = "gh-pages")
$repo = Split-Path $PSScriptRoot -Parent
$git = (Get-Command git -ErrorAction SilentlyContinue).Source
if (-not $git) { $git = (Get-ChildItem "$env:LOCALAPPDATA\GitHubDesktop\app-*\resources\app\git\cmd\git.exe" | Sort-Object FullName | Select-Object -Last 1).FullName }
function Invoke-Git { & $git @args 2>&1 | ForEach-Object { "$_" }; if ($LASTEXITCODE -ne 0) { throw "git $args failed" } }
if (-not (Test-Path "$repo\docs\grain-viewer\data\index.json")) { throw "docs/grain-viewer/data/index.json missing: run scripts/build_grain_viewer.py first" }

$stage = Join-Path ([IO.Path]::GetTempPath()) ("stl2fem-publish-" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory $stage | Out-Null
try {
    Copy-Item -Recurse -Force "$repo\docs\*" $stage
    Copy-Item -Force "$repo\docs\.nojekyll" $stage -ErrorAction SilentlyContinue
    if (-not (Test-Path "$stage\.nojekyll")) { New-Item -ItemType File "$stage\.nojekyll" | Out-Null }
    Set-Content -Encoding ascii (Join-Path $stage ".gitattributes") "*.gz binary"
    $url = & $git -C $repo remote get-url $Remote
    $source = & $git -C $repo rev-parse HEAD
    Invoke-Git -C $stage init -q -b $Branch
    Invoke-Git -C $stage config core.autocrlf false
    Invoke-Git -C $stage config user.name (& $git -C $repo config user.name)
    Invoke-Git -C $stage config user.email (& $git -C $repo config user.email)
    Invoke-Git -C $stage add -A
    Invoke-Git -C $stage commit -q -m "Publish docs/ from $($source.Substring(0,10)) with the current grain-viewer data"
    Invoke-Git -C $stage push --force $url "${Branch}:${Branch}"
    Write-Host "Published $Branch from $source"
} finally {
    Remove-Item -Recurse -Force $stage
}
