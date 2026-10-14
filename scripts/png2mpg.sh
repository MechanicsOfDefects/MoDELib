#!/bin/bash
# Encode consecutive img_0.png, img_1.png, ... beside this script.
# Usage: bash png2mpg.sh [output.mpg] [fps]
# Default: video.mpg at 25 fps. Existing output files are never overwritten.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ $# -gt 2 || ${1:-} == '--help' ]]; then
    echo "Usage: $0 [output.mpg] [fps: 24, 25, 30, 50, or 60]"
    exit 0
fi
output=${1:-video.mpg}
fps=${2:-25}
case "$fps" in
    24|25|30|50|60) ;;
    *) echo "Unsupported MPEG frame rate: $fps" >&2; exit 1 ;;
esac
if ! command -v ffmpeg >/dev/null 2>&1; then
    echo "ffmpeg is required but was not found in PATH." >&2
    exit 1
fi
if [[ -e "$output" || -L "$output" ]]; then
    echo "Output already exists: $output. Choose a different output name." >&2
    exit 1
fi
count=0
for filename in img_*.png; do
    if [[ -f "$filename" && "$filename" =~ ^img_([0-9]+)\.png$ ]]; then
        count=$((count + 1))
    fi
done
if [[ $count -eq 0 ]]; then
    echo "No img_<number>.png images found." >&2
    exit 1
fi
for ((k=0; k<count; ++k)); do
    if [[ ! -f "img_${k}.png" ]]; then
        echo "Missing img_${k}.png. Run renamePNG.sh on the completed image set first." >&2
        exit 1
    fi
done
staging=$(mktemp -d .png2mpg.XXXXXX)
trap 'rm -f -- "$staging/video.mpg"; rmdir -- "$staging"' EXIT
# Pad odd image dimensions by at most one pixel for YUV 4:2:0 encoding.
# qscale 1 gives high quality; MPEG-2 is explicitly selected for the MPG container.
ffmpeg -hide_banner -nostdin -n -framerate "$fps" -start_number 0 -i 'img_%d.png' \
    -frames:v "$count" -an -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' \
    -c:v mpeg2video -qscale:v 1 -pix_fmt yuv420p -f mpeg "$staging/video.mpg"
mv -n -- "$staging/video.mpg" "$output"
if [[ -e "$staging/video.mpg" ]]; then
    echo "Output appeared during encoding; refusing to overwrite: $output" >&2
    exit 1
fi
printf 'Created %s from %s frames at %s fps.\n' "$output" "$count" "$fps"
