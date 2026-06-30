# ddRAD rendogbs – in-silico fragment prediction validation

This repository contains the validation pipeline for [rendogbs], a tool that predicts in-silico 
the fragments produced by ddRAD/GBS library preparation (restriction enzyme digestion + size selection). 
Main job: check how well the predicted fragments agree with real sequencing data: 
WGS reads aligned to the reference genome are compared against the predicted fragment set to see how many 
fragments are actually supported by real reads, and how recovery success relates to fragment GC content.

The validation was run across **20 species** (fish, birds, insects, mollusks, plants). 
For each species, 10 SRA runs (accessions) were randomly selected, downloaded, aligned to the species
reference genome, and compared against the rendogbs prediction.

## Repository structure

```
.
├── <Species_1>/
│   ├── <accession_1>.csv        # validation summary output for one sample
│   ├── <accession_2>.csv
│   ├── ...
│   └── metadata.txt              # data source, reference genome, rendogbs parameters for this species
├── <Species_2>/
│   └── ...
├── ...
├── bwa_pe.sh                     # paired-end alignment (bwa mem + samtools sort/index)
├── bwa_se.sh                     # single-end alignment
├── validation_pe.sh              # validation script for paired-end BAMs
├── validation_se.sh              # validation script for single-end BAMs
├── enzyme_combs_count.py         # counts predicted fragments per enzyme combination
├── enzyme_combs_counts.csv       # output of enzyme_combs_count.py
├── enzyme_combs_counts_barplot.png
├── generate_fasta.py             # builds FASTA of predicted fragments for downstream use
└── rendogbs.R                    # R script for generating plot of GC content
```

Species directories include: *Anas platyrhynchos, Brassica napus, 
Camellia sinensis var. assamica, Coffea arabica, Crassostrea virginica, 
Fragaria × ananassa, Halyomorpha halys, Labeo rohita, Oncorhynchus mykiss, 
Oryza sativa, Populus tremula, Prunus persica, Quercus cerris, Quercus ilex, 
Quercus rubra, Salmo trutta, Scatophagus argus, Sesamum indicum, 
Solanum lycopersicum, Sparus aurata*.

Each species folder has its own `metadata.txt` documenting the source publication, 
the `fasterq-dump` command used to download the data, the reference genome source, 
and the exact `rendogbs.sh` invocation (restriction enzymes, size-selection range, 
number of chromosomes processed, etc.) used for that species.

## Workflow overview

```
SRA accessions → fasterq-dump → FASTQ
                                   │
              reference genome  ───┼── bwa-mem and samtools (bwa_pe.sh / bwa_se.sh) → sorted + indexed BAM
                                   │
        rendogbs.sh prediction  ───┘
        (filtered fragments .csv)
                                   │
                    validation_pe.sh / validation_se.sh
                                   │
        ┌──────────────────────────┴───────────────────────────┐
        matched_fragments.csv   unmatched_fragments.csv   gc_bins.tsv   fragment_coverage.tsv
                                                                   │
                                                    per-sample summary .csv
                                                    (the files inside each species folder)
```

## Files desription

### `metadata.txt` (per species)

Record for each species directory:
- citation of the source publication the samples/protocol were taken from,
- the randomly selected SRA accessions and the exact `fasterq-dump` command used to retrieve them,
- the reference genome source and how it was obtained (e.g. `wget` from NCBI),
- the exact `rendogbs.sh` command used to generate the fragment predictions for 
that species (`--ref`, `--size`, `--combinations`/`--combinations-file`, `--parallel`, `--chroms`, `--workdir`).

### `bwa_pe.sh`

Aligns **paired-end** FASTQ files (`<sample>_1.fastq` / `<sample>_2.fastq`) to a reference genome.

- Loops over a hardcoded `SAMPLES` array; 
each sample is processed as a background job (`&` ... `wait`)
- `bwa mem -t 4 -R "@RG\tID:...\tSM:...\tPL:ILLUMINA"` 
against `$REF`, piped directly into `samtools sort -@ 2` to produce a coordinate-sorted BAM.
- `samtools index` is run on the resulting BAM.
- `$REF`, `$OUTDIR`, `$DATADIR` and the `SAMPLES` array are left blank/placeholder 
in the script and must be filled in (or parameterized) before running it for a given species.

### `bwa_se.sh`

Same logic as `bwa_pe.sh` but for **single-end** reads: 
takes a hardcoded list of absolute FASTQ paths (`FILES` array) instead of a paired sample list, 
derives the sample name from the filename, and runs `bwa mem` with a single FASTQ input. 
Output handling (sort, index, parallel background jobs) is identical to the PE version.

### `validation_pe.sh` / `validation_se.sh`

Core validation script. Compares a set of predicted ddRAD fragments 
against an aligned BAM file and reports how many fragments are recovered by real reads, 
broken down by GC content.

```bash
./validation_pe.sh filtered.csv sample.bam genome.fa
# optional tolerance overrides for specific trimming methods:
START_TOL=10 END_TOL=10 ./validation_pe.sh filtered.csv sample.bam genome.fa
```

Inputs:
- `filtered.csv` – predicted fragments from rendogbs, columns: 
`accession, start_pos, start_enzyme, end_enzyme, fragment_length`
- `<bam>` – reads aligned to the reference (from `bwa_pe.sh`/`bwa_se.sh`)
- `<genome.fa>` – the reference genome (FASTA)

Steps performed:

1. **Indexing** – creates `genome.fa.fai` and `<bam>.bai` if missing.

2. **CSV → BED** (`awk`) – converts the predicted-fragment CSV into a sorted BED file (0-based coordinates, 
strand is set from `start_enzyme`).

3. **BAM → BED** (`samtools view -F 2308` + `awk`) – filters out unmapped/secondary/supplementary alignments, 
computes each read's reference-consumed length from its CIGAR string, and writes sorted read intervals.

4. **Coverage** (`bedtools coverage`) – fragment intervals 
(expanded by `START_TOL`/`END_TOL`, default 10 bp on each side) are intersected with the read BED 
to count reads per fragment and to report the overall on-target vs. off-target read fraction.

5. **Containment matching** (`bedtools intersect -wa -wb`) – a fragment is called **matched** 
if at least one read is fully contained within the fragment ± tolerance 
(`read_start ≥ frag_start − START_TOL` and `read_end ≤ frag_end + END_TOL`); otherwise it is **unmatched**.
6. **GC analysis** (embedded Python, via `bedtools getfasta`) – extracts each fragment's sequence, 
computes GC%, bins fragments into 5%-wide GC bins, and computes the recovery rate (matched / total) within each bin.

Outputs used mainly for debugging:

| file | content |
|---|---|
| `matched_fragments.csv` | fragments supported by reads (+ GC%, read count, status) |
| `unmatched_fragments.csv` | fragments with no supporting reads |
| `fragment_coverage.tsv` | all fragments sorted by read count |
| `gc_bins.tsv` | recovery rate per 5% GC bin |

The per-sample `.csv` files inside each species folder (e.g. `Labeo_rohita/SRR19358299.csv`) 
are a condensed summary derived from this output: total predicted fragments, matched/unmatched counts 
and percentages, total/on-target/off-target read counts, followed by the full `gc_bins.tsv` table.

`validation_pe.sh` and `validation_se.sh` share essentially the same matching/coverage/GC logic, 
the only difference is how the input BAM was produced (paired-end vs. single-end alignment) 
the validation itself works on coordinate intervals from the BAM regardless of library type.

### `enzyme_combs_count.py` / `enzyme_combs_counts.csv` / `enzyme_combs_counts_barplot.png`

Counts how many predicted fragments fall into each restriction-enzyme combination tested by rendogbs, 
exports the counts to `enzyme_combs_counts.csv`, and renders a bar plot (`enzyme_combs_counts_barplot.png`) 
summarizing fragment counts per enzyme combination.

### `generate_fasta.py`

Generates a FASTA file of the predicted fragment sequences 
(can be used as input for downstream steps such as GC computation or extraction outside of the 
`validation_*.sh` pipeline).

### `rendogbs.R`

TODO - Vašek

## Requirements

- `bwa`, `samtools`
- `bedtools`
- `python3`
- `sra-tools` (`fasterq-dump`) for data download
- `R` (for `rendogbs.R`)
- bash

## Reproducing the pipeline for one species

```bash
# 1. download data per that species' metadata.txt
fasterq-dump <accessions...> -O data --split-files

# 2. index reference fasta
bwa index xxx.fasta

# 2. align reads (fill in REF/OUTDIR/DATADIR and the sample list inside the script first)
./bwa_pe.sh    # or bwa_se.sh for single-end data

# 3. predict fragments (rendogbs, see metadata.txt for the exact parameters used)
sh ./rendogbs.sh --ref genome.fna --size <min-max> --combinations ...

# 4. validate predictions against the alignment
./validation_pe.sh filtered.csv sample.bam genome.fa   # or validation_se.sh
```