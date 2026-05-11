#!/usr/bin/env bash
set -euo pipefail

INPUT=""
OUTPUT=""
COLUMN="noncodingRNA_fam"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT="$2"; shift 2;;
    --output)
      OUTPUT="$2"; shift 2;;
    --column)
      COLUMN="$2"; shift 2;;
    -h|--help)
      echo "Usage: $0 --input INPUT --output OUTPUT [--column NAME]"; exit 0;;
    *)
      echo "Unknown arg: $1"; exit 1;;
  esac
done

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
  echo "--input and --output are required" >&2
  exit 2
fi

if [[ ! -f "$INPUT" ]]; then
  echo "Input file not found: $INPUT" >&2
  exit 3
fi

mkdir -p "$(dirname "$OUTPUT")"

# find 1-based column index
COLUMN_NUMBER=$(head -n 1 "$INPUT" | tr '\t' '\n' | nl -ba | awk -v col="$COLUMN" '$2 == col {print $1}')

if [[ -z "$COLUMN_NUMBER" ]]; then
  echo "Error: column '$COLUMN' not found in header of $INPUT" >&2
  exit 4
fi

# write header then sorted body using a stable locale (reproducible across systems)
(head -n 1 "$INPUT" && tail -n +2 "$INPUT" | LC_ALL=C sort -t $'\t' -k "${COLUMN_NUMBER},${COLUMN_NUMBER}") > "$OUTPUT"

exit 0
