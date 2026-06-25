#!/bin/bash

REF=""
OUTDIR=""

mkdir -p "$OUTDIR"

FILES=(
"/home/raidfolder/rendogbs/data/SRR5359996_1.fastq"
"/home/raidfolder/rendogbs/data/SRR5360014_1.fastq"
"/home/raidfolder/rendogbs/data/SRR5360050_1.fastq"
"/home/raidfolder/rendogbs/data/SRR5360054_1.fastq"
"/home/raidfolder/rendogbs/data/SRR5360061_1.fastq"
"/home/raidfolder/rendogbs/data/SRR5360070_1.fastq"
"/home/raidfolder/rendogbs/data/SRR5360078_1.fastq"
"/home/raidfolder/rendogbs/data/SRR5360088_1.fastq"
)

for f in "${FILES[@]}"; do
(
    sample=$(basename "$f" .fastq)

    echo "[$(date)] Starting $sample"

    bwa mem \
        -t 4 \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
        "$REF" "$f" \
    | samtools sort \
        -@ 2 \
        -o "$OUTDIR/${sample}.bam"

    samtools index "$OUTDIR/${sample}.bam"

    echo "[$(date)] Finished $sample"
) &
done

wait

echo "All alignments completed."