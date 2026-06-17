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
    --fast-combinations)
      cat <<'EOF'
EcoRI      + MseI
EcoRI      + MboI
EcoRI      + NlaIII
PstI       + MspI
PstI       + MseI
PstI       + TaqI
PstI       + MboI
PstI       + HpaII
SbfI       + MspI
SbfI       + MseI
ApeKI      + MseI
BamHI      + MseI
HindIII    + MseI
HindIII    + NlaIII
NheI       + MboI
NheI       + MseI
XbaI       + MboI
SpeI       + MboI
KpnI       + MboI
NcoI       + MseI
NcoI       + MboI
BglII      + MseI
ClaI       + MboI
EOF
      exit 0
      ;;
    --enzymes)
      cat <<'EOF'
AanI        AvaII       BseYI       BstYI       Eco72I      Hpy99I      NcoI        PvuI        SphI
AatII       AvrII       BsiEI       BstZ17I     EcoRI       HpyCH4III   NdeI        PvuII       SrfI
Acc65I      BaeGI       BsiHKAI     Bsu36I      EcoRV       HpyCH4IV    NgoMIV      RsrII       SspI
AccI        BamHI       BsiWI       BtgI        EheI        HpyCH4V     NheI        RsaI        StyD4I
AciI        BanI        BsoBI       ClaI        FatI        KasI        NlaIII      SacI        StyI
AclI        BanII       Bsp1286I    CpoI        Fnu4HI      KpnI        NotI        SacII       SwaI
AfeI        BbvCI       BspDI       CsiI        FseI        MboI        NruI        SalI        TaiI
AflII       BclI        BspEI       CviKI-1     FspAI       MluCI       NsiI        Sau3AI      TaqI
AflIII      BfaI        BspHI       CviQI       FspI        MluI        NspI        Sau96I      TasI
AgeI        BfoI        BspOI       DdeI        HaeII       MreI        PacI        SbfI        TatI
AluI        BglII       BsrBI       DpnI        HaeIII      MscI        PaeR7I      ScaI        TauI
ApaI        BlpI        BsrFI       DpnII       HhaI        MseI        PciI        ScrFI       TfiI
ApaLI       BmgBI       BsrGI       DraI        HinP1I      MspA1I      PluTI       SexAI       TseI
ApeKI       Bpu10I      BssHII      EaeI        HincII      MspI        PpuMI       SfcI        Tsp45I
ApoI        BsaAI       BssSI       EagI        HindIII     MssI        PsiI        SgrAI       TspMI
AscI        BsaHI       BstBI       Ecl136II    HinfI       MunI        PspGI       SgrDI       XbaI
AseI        BsaJI       BstEII      Eco105I     HpaI        NaeI        PspOMI      SmaI        XhoI
AsiSI       BsaWI       BstNI       Eco147I     HpaII       NarI        PspXI       SmlI        XmaI
AvaI        BseSI       BstUI       Eco53kI     Hpy188I     NciI        PstI        SpeI        ZraI
EOF
      exit 0
      ;;

    -h|--help)
      cat <<'EOF'
rendogbs  -  In-silico ddRAD/GBS library pipeline

Runs diges and analysis, then generates plots
(TSV summary + figures) sequentially.  Plots start only after the pipeline
exits successfully.  All output is written to <workdir>/results/.

CONTAINER TYPE (default: --docker)
  --docker                  Run the docker container (default)
  --singularity             Run the singularity container rendogbs-v1.sif in this
                            directory
  --image          <image>  Docker image name or singularity image file name.

PBS SUPPORT
  --qsub                    Do not run directly but submit as PBS job using qsub
  --limits|-l     <limits>  Specify arbitrary PBS job limits (typically mem=XXgb)
  --interactive|-I          Run as interactive PBS job
  --name|-N         <name>  Specify PBS job name (defaults to rendogbs)

REQUIRED
  --ref            <file>   Reference FASTA (.fa / .fasta / .fa.gz)
  --workdir        <dir>    Working directory (created if absent)
  --parallel       <int>    Combinations processed in parallel per batch, don't more threads than combinations
  --size           <range>  Fragment size window, e.g. 200-400 (both inclusive).
                            Standard 100 bp bins (0-99 .. 900-999 + >=1000) are
                            always reported; filtered.csv retains only fragments
                            within --size-range.

COMBINATIONS  (optional, default: fast)
  --combinations   fast     23 built-in common ddRAD/GBS pairs (default)
                   custom   read from --combinations-file
  --combinations-file <f>   CSV for custom mode (header line: enzyme_a,enzyme_b),
                            please reffer to --enzymes to choose enzymes from the list 
                            of 171 available enzymes.

ANNOTATION  (optional)
  --annotation     <file>   GFF3/GFF/GTF gene annotation
  --te             <file>   RepeatMasker .out TE annotation

PLOT OPTIONS  (optional)
  --chroms         <int>    Longest N contigs treated as chromosomes in
                            per-chromosome plots (default: 10)
  --dpi            <int>    Figure resolution in DPI (default: 300)
  --skip-plots              Run pipeline only, skip plot generation

OTHER
  -h, --help                Show this help and exit
  --enzymes                 Show available enzymes and exit
  --fast-combinations       Show fast combinations and exit

OUTPUT STRUCTURE
  workdir/
    results/
      contig_lengths.txt          contigs sorted descending by length
      summary.tsv                 Excel-ready: one row per combination,
                                  all bin counts + GC metrics + annotation %
      EcoRI_MseI/
        cuts.csv                  all cut sites (accession, position, enzyme)
        fragments.csv             all adjacent pairs of different enzymes
        filtered.csv              fragments within --size window and with resolved multi-cutter regions
        distribution.csv          standard 100 bp bins + custom window row
        gc_metrics.csv            GC statistics for filtered fragments
        annotation_summary.csv    TE / gene coverage (only with --te / --annotation)
        uncertain_cuts.csv        fragments with ambiguous cut sites
        statistics_uncertain.csv  statistics for ambiguous cut sites
        statistics_cutting.csv    statistics for all cut sites
      AciI_HindIII/
        ...
      plots/
        heatmap_fragment_lengths.png    fragment count heatmap (10 bp bins, log scale)
        heatmap_chrom_distribution.png  filtered fragments per chromosome (heatmap)
        bar_chrom_distribution.png      filtered fragments per chromosome (bar chart)
        gc_distribution.png             GC% per combination: mean +/- SD
        annotation_coverage.png         annotation category coverage stacked bar
        size_distributions.png          100 bp bin line plot, user window highlighted

EXAMPLES
  # Don't use more threads than combinations, pipeline won't finish faster
  # Basic run, fast preset, 200-400 bp window
  rendogbs.sh --ref genome.fasta --workdir ./run1 --parallel 11 --size 200-400

  # All combinations, with annotation, 11 chromosomes in plots
  rendogbs.sh --ref genome.fasta --workdir ./run1 --parallel 11 --size 150-350 \
              --combinations fast \
              --annotation genes.gff3.gz --te repeats.out \
              --chroms 11

  # Custom enzyme pairs, pipeline only
  rendogbs.sh --ref genome.fa --workdir ./run1 --parallel 4 --size 200-500 \
              --combinations custom --combinations-file my_pairs.csv \
              --skip-plots

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