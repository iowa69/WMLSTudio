#!/usr/bin/env bash
# Native PE build; never invokes WSL, Docker, Linux binaries or emulation.
set -euo pipefail
if [[ "${MSYSTEM:-}" != UCRT64 ]]; then
    echo "Run this script in an MSYS2 UCRT64 shell." >&2
    exit 2
fi
recipe_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$recipe_dir/.." && pwd)"
commit=c1413581e4f37211892d3c4310d01f3d9a9b3490
destination="${1:-$repo_dir/artifacts/Tools/skesa}"
if [[ -e "$destination" ]]; then
    echo "Destination already exists: $destination. Choose a new directory." >&2
    exit 2
fi
mkdir -p -- "$(dirname -- "$destination")"
build_dir="$(mktemp -d -t wmlstudio-skesa-build-XXXXXXXX)"
# Keep failed build/source/logs for diagnosis. Never recursively delete user paths.
echo "Build directory: $build_dir"
git init "$build_dir/source"
git -C "$build_dir/source" remote add origin https://github.com/ncbi/SKESA.git
git -C "$build_dir/source" fetch --depth=1 origin "$commit"
git -C "$build_dir/source" checkout --detach FETCH_HEAD
[[ "$(git -C "$build_dir/source" rev-parse HEAD)" == "$commit" ]]
git -C "$build_dir/source" apply --check "$recipe_dir/skesa_windows.patch"
git -C "$build_dir/source" apply "$recipe_dir/skesa_windows.patch"
mkdir "$build_dir/stage"
# NO_NGS omits online SRA. No -march=native: retain upstream's SSE4.2 baseline.
g++ -std=c++14 -DNO_NGS -DBOOST_ALL_DYN_LINK -O3 -msse4.2 -pthread -Wall -Wno-format-y2k \
    -I"$build_dir/source" "$build_dir/source/skesa.cpp" "$build_dir/source/glb_align.cpp" \
    -o "$build_dir/stage/skesa.exe" \
    -lboost_program_options-mt -lboost_iostreams-mt -lboost_timer-mt \
    -lboost_chrono-mt -lz 2>&1 | tee "$build_dir/compile.log"
python "$recipe_dir/stage_skesa_runtime.py" \
    --executable "$build_dir/stage/skesa.exe" --runtime-bin /ucrt64/bin \
    --source "$build_dir/source" --recipe-dir "$recipe_dir" --commit "$commit"
cp -- "$build_dir/compile.log" "$build_dir/stage/compile.log"
python "$recipe_dir/check_skesa.py" "$build_dir/stage/skesa.exe" \
    --output "$build_dir/stage/smoke-test.json" --adapter-source "$repo_dir"
# Publish only after a real assembly and validation of the runtime DLL closure.
mv -T -- "$build_dir/stage" "$destination"
echo "Native SKESA 2.4.0 ready: $destination"
