#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/target/release/rendogbs_finder" \
  --enzymes-run "$SCRIPT_DIR/enzymes_run.csv" \
  --ref "/Users/eliskakorbova/Desktop/biopython/rust3/genome.fasta" \
  --out-dir "/Users/eliskakorbova/Desktop/biopython/rust3/results" \
  --parallel 4
