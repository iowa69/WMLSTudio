param([string]$SchemeSource = "", [string]$HydraSource = "", [string]$SkesaSource = "", [switch]$SkipTests)
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
    if ($HydraSource) {
        & uv run python studio_packaging/stage_bio_tools.py --platform windows-x64 --hydra-source $HydraSource
    } else {
        & uv run python studio_packaging/stage_bio_tools.py --platform windows-x64 --download-hydra-starter
    }
    if ($LASTEXITCODE -ne 0) { throw "Native BLAST tool staging failed" }
    $skesaArguments = @("studio_packaging/stage_bio_tools.py", "--platform", "windows-x64", "--require-skesa")
    if ($SkesaSource) { $skesaArguments += @("--skesa-source", $SkesaSource) }
    & uv run python @skesaArguments
    if ($LASTEXITCODE -ne 0) { throw "A verified native SKESA bundle is required. Supply -SkesaSource from the Windows SKESA workflow artifact." }
    & uv run python studio_packaging/stage_notices.py
    if ($LASTEXITCODE -ne 0) { throw "License text staging failed" }
    if (-not $SkipTests) {
        $env:QT_QPA_PLATFORM = "offscreen"
        & uv run pytest
        if ($LASTEXITCODE -ne 0) { throw "Tests failed" }
    }
    & uv run pyinstaller --noconfirm studio_packaging/wmlstudio.spec
    if ($LASTEXITCODE -ne 0) { throw "Freezing failed" }
    & uv run python studio_packaging/audit_windows_runtime.py dist/WMLSTudio --output dist/windows-runtime-closure.json
    if ($LASTEXITCODE -ne 0) { throw "Portable Windows dependency closure failed; system-installed redistributables are not accepted" }
    $env:QT_QPA_PLATFORM = "offscreen"
    & uv run python studio_packaging/check_frozen.py dist/WMLSTudio --output dist/frozen-check.json
    if ($LASTEXITCODE -ne 0) { throw "Frozen reference typing and desktop checks failed" }
    Compress-Archive -Path dist/WMLSTudio -DestinationPath dist/WMLSTudio-Windows-x64.zip -Force
    Get-FileHash -Algorithm SHA256 dist/WMLSTudio-Windows-x64.zip
} finally {
    Pop-Location
}
