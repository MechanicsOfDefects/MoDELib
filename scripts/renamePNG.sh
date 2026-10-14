#!/bin/bash
# Renumber img_<runID>.png in numeric runID order, in this script's directory.
# Use --dry-run to preview. A successful run saves the old/new filename mapping.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--dry-run" ) ]]; then
    echo "Usage: $0 [--dry-run]" >&2
    exit 1
fi

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

if [[ ${#files[@]} -eq 0 ]]; then
    echo "No img_<runID>.png files found."
    exit 0
fi

for ((k=0; k<${#files[@]}; ++k)); do
    if [[ -d "img_${k}.png" ]]; then
        echo "Cannot replace directory img_${k}.png" >&2
        exit 1
    fi
    printf '%s -> img_%s.png\n' "${files[k]}" "$k"
done
if [[ ${1:-} == "--dry-run" ]]; then
    exit 0
fi

# Stage ALL sources before assigning final names, preventing collisions.
staging=$(mktemp -d .renamePNG.XXXXXX)
trap 'echo "Rename interrupted; remaining images and mapping are in $staging" >&2' ERR
printf 'original\trenamed\n' > "$staging/mapping.tsv"
for ((k=0; k<${#files[@]}; ++k)); do
    printf '%s\timg_%s.png\n' "${files[k]}" "$k" >> "$staging/mapping.tsv"
    mv -- "${files[k]}" "$staging/img_${k}.png"
done
for ((k=0; k<${#files[@]}; ++k)); do
    mv -- "$staging/img_${k}.png" "img_${k}.png"
done
mapping="renamePNG-map-${staging##*.}.tsv"
mv -- "$staging/mapping.tsv" "$mapping"
rmdir -- "$staging"
trap - ERR
printf 'Renamed %s images. Original names recorded in %s\n' "${#files[@]}" "$mapping"
