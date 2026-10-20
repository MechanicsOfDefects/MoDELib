#!/bin/bash
# Example (from the repository root):
# bash scripts/png2mp4.sh tutorials/inclusions/config
#
# Encode img_<runID>.png in the supplied folder in numeric runID order.
# Relative output paths are resolved inside folderContainingImages.
# IDs may be sparse or start above zero; original images are never renamed.
# Each image occupies one frame at the selected fps, regardless of runID gaps.
# Usage: bash png2mp4.sh folderContainingImages [output.mp4] [fps] [crf]
# Default: video.mp4 at 25 fps, CRF 20 (lower = higher quality/larger file). Existing output files are never overwritten.
set -euo pipefail
if [[ ${1:-} == '--help' ]]; then
    echo "Usage: bash png2mp4.sh folderContainingImages [output.mp4] [fps] [crf]"
    exit 0
fi
if [[ $# -lt 1 || $# -gt 4 ]]; then
    echo "Usage: bash png2mp4.sh folderContainingImages [output.mp4] [fps] [crf]" >&2
    exit 1
fi
if [[ ! -d "$1" ]]; then
    echo "Image folder does not exist: $1" >&2
    exit 1
fi
cd -- "$1"
output=${2:-video.mp4}
fps=${3:-25}
crf=${4:-20}
if [[ ! "$crf" =~ ^([0-9]|[1-4][0-9]|5[01])$ ]]; then
    echo "CRF must be an integer from 0 (lossless) to 51 (lowest quality)." >&2
    exit 1
fi
case "$fps" in
    24|25|30|50|60) ;;
    *) echo "Unsupported frame rate: $fps" >&2; exit 1 ;;
esac
if ! command -v ffmpeg >/dev/null 2>&1; then
    echo "ffmpeg is required but was not found in PATH." >&2
    exit 1
fi
if [[ -e "$output" || -L "$output" ]]; then
    echo "Output already exists: $output. Choose a different output name." >&2
    exit 1
fi
# Extract the whole numeric ID; lexical filename sorting puts 100 before 20.
files=()
while IFS=$'\t' read -r runID filename; do
    files+=("$filename")
done < <(
    for filename in img_*.png; do
        if [[ -f "$filename" && "$filename" =~ ^img_([0-9]+)\.png$ ]]; then
            printf '%s\t%s\n' "${BASH_REMATCH[1]}" "$filename"
        fi
    done | LC_ALL=C sort -k1,1n -k2,2
)
count=${#files[@]}
if [[ $count -eq 0 ]]; then
    echo "No img_<runID>.png images found." >&2
    exit 1
fi
staging=$(mktemp -d .png2mp4.XXXXXX)
trap 'rm -f -- "$staging"/img_*.png "$staging/video.mp4"; rmdir -- "$staging"' EXIT
# Sequential temporary symlinks let FFmpeg read sparse IDs at a constant fps.
for ((k=0; k<count; ++k)); do
    ln -s -- "../${files[k]}" "$staging/img_${k}.png"
done
# Pad odd image dimensions by at most one pixel for YUV 4:2:0 encoding.
# H.264 CRF 20 preserves detail; slow improves compression without reducing resolution.
# faststart places the MP4 index at the beginning for quicker playback.
ffmpeg -hide_banner -nostdin -n -framerate "$fps" -start_number 0 -i "$staging/img_%d.png" \
    -frames:v "$count" -an -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' \
    -c:v libx264 -preset slow -crf "$crf" -pix_fmt yuv420p -movflags +faststart -f mp4 "$staging/video.mp4"
mv -n -- "$staging/video.mp4" "$output"
if [[ -e "$staging/video.mp4" ]]; then
    echo "Output appeared during encoding; refusing to overwrite: $output" >&2
    exit 1
fi
printf 'Created %s from %s frames at %s fps.\n' "$output" "$count" "$fps"
