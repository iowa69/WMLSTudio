param([string]$SchemeSource = "", [switch]$SkipTests)
$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
Push-Location $projectRoot
try {
    if (-not $IsWindows -and $env:OS -ne "Windows_NT") {
        throw "Build the Windows package on Windows 11 or a Windows GitHub runner."
    }
    & uv sync --locked --group dev
    if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed" }
    if ($SchemeSource) {
        & uv run python studio_scripts/stage_schemes.py --source $SchemeSource
    } else {
        & uv run python studio_scripts/stage_schemes.py --download
    }
    if ($LASTEXITCODE -ne 0) { throw "Scheme staging failed" }
    & uv run python studio_packaging/stage_notices.py
    if ($LASTEXITCODE -ne 0) { throw "License text staging failed" }
    if (-not $SkipTests) {
        $env:QT_QPA_PLATFORM = "offscreen"
        & uv run pytest
        if ($LASTEXITCODE -ne 0) { throw "Tests failed" }
    }
    & uv run pyinstaller --noconfirm studio_packaging/wmlstudio.spec
    if ($LASTEXITCODE -ne 0) { throw "Freezing failed" }
    $env:QT_QPA_PLATFORM = "offscreen"
    & uv run python studio_packaging/check_frozen.py dist/WMLSTudio --output dist/frozen-check.json
    if ($LASTEXITCODE -ne 0) { throw "Frozen reference typing and desktop checks failed" }
    Compress-Archive -Path dist/WMLSTudio -DestinationPath dist/WMLSTudio-Windows-x64.zip -Force
    Get-FileHash -Algorithm SHA256 dist/WMLSTudio-Windows-x64.zip
} finally {
    Pop-Location
}
