#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="$ROOT/manifests/samples_basic.tsv"
COHORT="article1_basic_full_raw"
OUTDIR="$ROOT/raw/$COHORT"
FORCE=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [--manifest PATH] [--cohort NAME] [--outdir PATH] [--force]

Download FASTQ files listed in a sample manifest and concatenate lane-split runs
into one merged FASTQ per sample label.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --cohort) COHORT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "$OUTDIR/runs" "$OUTDIR/merged"

echo "Manifest: $MANIFEST"
echo "Cohort:   $COHORT"
echo "Outdir:   $OUTDIR"

awk -F '\t' -v cohort="$COHORT" 'NR>1 && $1==cohort && $2=="yes" {print}' "$MANIFEST" |
while IFS=$'\t' read -r cohort selected sample_label condition replicate library cell_line run_accessions fastq_urls note; do
  IFS=';' read -r -a runs <<< "$run_accessions"
  IFS=';' read -r -a urls <<< "$fastq_urls"

  if [[ ${#runs[@]} -ne ${#urls[@]} ]]; then
    echo "Mismatched run/url counts for $sample_label" >&2
    exit 1
  fi

  files=()
  for i in "${!runs[@]}"; do
    run="${runs[$i]}"
    url="${urls[$i]}"
    [[ "$url" =~ ^(ftp|https?):// ]] || url="https://$url"
    dest="$OUTDIR/runs/${run}.fastq.gz"
    if [[ $FORCE -eq 1 || ! -s "$dest" ]]; then
      echo "[download] $run -> $dest"
      curl -L --retry 5 --retry-delay 5 --continue-at - -o "$dest" "$url"
    else
      echo "[skip] $dest"
    fi
    files+=("$dest")
  done

  merged="$OUTDIR/merged/${sample_label}.fastq.gz"
  if [[ $FORCE -eq 1 || ! -s "$merged" ]]; then
    echo "[merge] ${sample_label}"
    cat "${files[@]}" > "$merged"
  else
    echo "[skip] $merged"
  fi

  echo "[check] gzip -t $merged"
  gzip -t "$merged"
done
