#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "RendoGBS pipeline started"
echo "Working dir: $PWD"
echo "Binary dir: $DIR"

# defaults
WORKDIR=""
REF=""
PARALLEL="2"
SIZE=""
CHROMS=""
DPI=""
COMBINATIONS="fast"
COMBINATIONS_FILE=""
TE=""
ANNOTATION=""
SKIP_PLOTS=0

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir)
      WORKDIR="${2:-}"
      shift 2
      ;;
    --ref)
      REF="${2:-}"
      shift 2
      ;;
    --parallel)
      PARALLEL="${2:-}"
      shift 2
      ;;
    --size)
      SIZE="${2:-}"
      shift 2
      ;;
    --chroms)
      CHROMS="${2:-}"
      shift 2
      ;;
    --dpi)
      DPI="${2:-}"
      shift 2
      ;;
    --te)
      TE="${2:-}"
      shift 2
      ;;
    --annotation)
      ANNOTATION="${2:-}"
      shift 2
      ;;
    --combinations)
      COMBINATIONS="${2:-}"
      shift 2
      ;;
    --combinations-file)
      COMBINATIONS_FILE="${2:-}"
      shift 2
      ;;
    --skip-plots)
      SKIP_PLOTS=1
      shift
      ;;
    -h|--help)
      cat <<'EOF'
Placeholder
EOF
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

# check

if [[ -z "$WORKDIR" ]]; then
  echo "[ERROR] --workdir is required" >&2
  exit 1
fi

if [[ -z "$REF" ]]; then
  echo "[ERROR] --ref is required" >&2
  exit 1
fi

if [[ -z "$PARALLEL" ]]; then
  echo "[ERROR] --parallel is required" >&2
  exit 1
fi

if [[ -z "$COMBINATIONS" ]]; then
  echo "[ERROR] --combinations is required" >&2
  exit 1
fi

if [[ "$COMBINATIONS" != "fast" && "$COMBINATIONS" != "custom" ]]; then
  echo "[ERROR] --combinations must be 'fast' or 'custom'" >&2
  exit 1
fi

if [[ "$COMBINATIONS" == "custom" && -z "$COMBINATIONS_FILE" ]]; then
  echo "[ERROR] --combinations-file is required when --combinations custom" >&2
  exit 1
fi

if [[ -z "$SIZE" ]]; then
  echo "[ERROR] --size is required" >&2
  exit 1
fi

if [[ "$SKIP_PLOTS" -eq 0 ]]; then
  if [[ -z "$CHROMS" ]]; then
    echo "[ERROR] --chroms is required unless --skip-plots is used" >&2
    exit 1
  fi
fi

mkdir -p "$WORKDIR"

# 1) combinations_processing.py
echo "Step 1: combination_processing.py"

if [[ "$COMBINATIONS" == "custom" ]]; then
  python3 "$DIR/combination_processing.py" \
    --workdir "$WORKDIR" \
    --combinations custom \
    --combinations-file "$COMBINATIONS_FILE"
else
  python3 "$DIR/combination_processing.py" \
    --workdir "$WORKDIR" \
    --combinations fast
fi

# 2) rendogbs_finder (Rust)
echo "Step 2: rendogbs_finder"

"$DIR/rendogbs_finder" \
  --out-dir "$WORKDIR/results/cuts/" \
  --ref "$REF" \
  --parallel "$PARALLEL" \
  --enzymes-run "$WORKDIR/results/enzymes_run.csv"

# 3) cut_merge (Rust)
echo "Step 3: cut_merge"

"$DIR/cut_merge" \
  --out-dir "$WORKDIR/results/" \
  --parallel "$PARALLEL" \
  --combinations "$WORKDIR/results/combinations.csv" \
  --cuts-dir "$WORKDIR/results/cuts"

# 4) fragment_generation.py
echo "Step 4: fragment_generation.py"

python3 "$DIR/fragment_generation.py" \
  --workdir "$WORKDIR" \
  --parallel "$PARALLEL" \
  --size "$SIZE"

# 5) fragment_generation (Rust)
echo "Step 5: fragment_generation"

"$DIR/fragment_generation" \
  --out-dir "$WORKDIR/results/" \
  --parallel "$PARALLEL" \
  --combinations "$WORKDIR/results/combinations.csv" \
  --size "$SIZE"

# 6) postprocess_metrics.py
echo "Step 6: postprocess_metrics.py"

python3 "$DIR/postprocess_metrics.py" \
  --workdir "$WORKDIR" \
  --ref "$REF" \
  --parallel "$PARALLEL" \
  --size "$SIZE"

# 7) annotation.py
if [[ -n "$ANNOTATION" || -n "$TE" ]]; then
  echo "Step 7: annotation.py"

  ann_args=(--workdir "$WORKDIR" --parallel "$PARALLEL")

  if [[ -n "$TE" ]]; then
    ann_args+=(--te "$TE")
  fi

  if [[ -n "$ANNOTATION" ]]; then
    ann_args+=(--annotation "$ANNOTATION")
  fi

  python3 "$DIR/annotation.py" "${ann_args[@]}"
else
  echo "Step 7: annotation.py skipped (no --te / --annotation provided)"
fi

# 8) summary.py
echo "Step 8: summary.py"

python3 "$DIR/summary.py" \
  --workdir "$WORKDIR"

# 9) plots.py
if [[ "$SKIP_PLOTS" -eq 1 ]]; then
  echo "Step 9: plots.py skipped (--skip-plots)"
else
  echo "Step 9: plots.py"

  plot_args=(--workdir "$WORKDIR" --size "$SIZE" --chroms "$CHROMS")
  if [[ -n "$DPI" ]]; then
    plot_args+=(--dpi "$DPI")
  fi

  python3 "$DIR/plots.py" "${plot_args[@]}"
fi

echo "Pipeline finished successfully"