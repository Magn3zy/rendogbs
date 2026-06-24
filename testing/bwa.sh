#!/bin/bash

REF=""
OUTDIR=""
DATADIR=""

mkdir -p "$OUTDIR"

SAMPLES=(
"ERR11741110"
"ERR11741153"
"ERR11741120"
"ERR11741148"
"ERR11741113"
)

for s in "${SAMPLES[@]}"; do
(
    R1="$DATADIR/${s}_1.fastq"
    R2="$DATADIR/${s}_2.fastq"

    echo "[$(date)] Starting $s"

    bwa mem \
        -t 4 \
        -R "@RG\tID:${s}\tSM:${s}\tPL:ILLUMINA" \
        "$REF" "$R1" "$R2" \
    | samtools sort \
        -@ 2 \
        -o "$OUTDIR/${s}.bam"

    samtools index "$OUTDIR/${s}.bam"

    echo "[$(date)] Finished $s"
) &
done

wait

echo "All alignments completed."
