param([string]$SchemeSource = "", [string]$HydraSource = "", [string]$SkesaSource = "", [string]$SkaSource = "", [switch]$SkipTests)
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
    & uv run python studio_scripts/stage_characterization.py
    if ($LASTEXITCODE -ne 0) { throw "Characterization reference staging failed" }
    & uv run python studio_packaging/stage_fastqc.py --platform windows-x64
    if ($LASTEXITCODE -ne 0) { throw "FastQC and private Java staging failed" }
    $env:WMLSTUDIO_TEST_FASTQC_ROOT = Join-Path $projectRoot "src/wmlstudio/resources/tools/fastqc"
    $skaArguments = @("studio_packaging/stage_ska.py", "--platform", "windows-x64")
    if ($SkaSource) { $skaArguments += @("--source", $SkaSource) }
    & uv run python @skaArguments
    if ($LASTEXITCODE -ne 0) { throw "A verified native SKA2 bundle is required. Supply -SkaSource from the native Windows SKA2 workflow artifact." }
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
    # The download carries its version, so a saved ZIP can still be identified
    # months later and two builds cannot be confused for one another.
    $version = & uv run python -c "import wmlstudio; print(wmlstudio.__version__)"
    if ($LASTEXITCODE -ne 0) { throw "Could not read the application version" }
    $archive = "dist/WMLSTudio-$version-Windows-x64.zip"
    Compress-Archive -Path dist/WMLSTudio -DestinationPath $archive -Force
    Get-FileHash -Algorithm SHA256 $archive
} finally {
    Pop-Location
}
