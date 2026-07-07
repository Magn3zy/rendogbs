# rendogbs — Core Architecture Documentation

## Table of Contents
1. [combinations_processing.py](#1-combinations_processingpy)
2. [rendogbs_finder (Rust)](#2-rendogbs_finder--rust-binary)
3. [cut_merge (Rust)](#3-cut_merge--rust-binary)
4. [fragment_generation (Rust)](#4-fragment_generation--rust-binary)
5. [fragment_generation.py](#5-fragment_generationpy)
6. [postprocess_metrics.py](#6-postprocess_metricspy)
7. [annotation.py](#7-annotationpy)
8. [summary.py](#8-summarypy)
9. [plots.py](#9-plotspy)

---

# combinations_processing.py 

## Summary
First step of the rendogbs pipeline. Validates and resolves enzyme combinations
before any computationally expensive steps run. Exits early on any invalid input
with a clear message.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory. `results/` is created. |
| `--combinations` | no (default: `fast`) | `fast` uses the 23 predefined combinations. `custom` reads from `--combinations-file`. |
| `--combinations-file` | only if `--combinations custom` | Path to user-supplied CSV with two columns: `enzyme_a`, `enzyme_b`. |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Input files
Both files are bundled with the pipeline (in `scripts/`) and copied automatically
by the wrapper:

- `enzymes.csv` — master enzyme database with expanded IUPAC sequences and cut offsets.
  Contains all supported enzymes (~171 entries) for nonpalindromic ones contains reverse complement.
- `allowed_pairs.csv` — precomputed list of valid enzyme pairs. A pair is allowed only
  if their recognition sites do not fully overlap, which would make simulation impossible.

## Validation logic
Each combination goes through two independent checks:

**1) Enzyme presence check** — both enzymes must exist in `enzymes.csv`.
If a user specifies an enzyme name with a typo or an unsupported enzyme, it is
skipped with a warning rather than crashing the pipeline. We use the most 
common name for each enzyme refer to --enzymes as part of help.

**2) Allowed pairs check** — the pair must exist in `allowed_pairs.csv`.
This check rejects combinations where the recognition sites of both enzymes fully
overlap. Such combinations cannot be simulated reliably because it is ambiguous
which enzyme cuts at a given position. This is a biological constraint, not a
software limitation.

Both checks are intentionally non-fatal (warnings, not errors) so that a batch
of mostly valid combinations is not completely rejected because of one bad entry.
The pipeline exits only if zero valid combinations remain after filtering.

## Fast combinations
The 23 hardcoded `FAST_COMBOS` represent commonly used ddRAD/GBS enzyme pairs
from the literature. They are provided as a convenience. Users who want to 
explore other combinations use `--combinations custom` we support 14 040 combinations.

## Output files (written to `<workdir>/results/`)
- `enzymes_run.csv` — subset of `enzymes.csv` containing only the unique enzymes
  needed for this run. Passed to `rendogbs_finder`.
- `combinations.csv` — validated pairs in `enzyme_a,enzyme_b` format. Used by
  all downstream scripts and binaries as the authoritative list of what to process.

---

# rendogbs_finder (Rust binary)

## Summary
Finds all recognition site cut positions for each enzyme across the entire reference
genome. Output is one CSV per enzyme containing every cut site position. This is the
most computationally intensive step in the pipeline.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--enzymes-run` | yes | Path to `enzymes_run.csv` produced by `combinations_processing.py` |
| `--ref` / `-r` | yes | Reference genome FASTA (plain or gzipped) |
| `--out-dir` / `-o` | yes | Output directory for per-enzyme CSV files (`results/cuts/`) |
| `--parallel` / `-p` | no (default: 2) | Total threads. One for the writer thread, rest search. |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Forward strand only
Recognition sites of restriction enzymes used in ddRAD/GBS are mostly palindromic —
the sequence on the reverse strand is the reverse complement of the forward strand
and identical in terms of cut position. Scanning only the forward strand finds
all cut sites without duplication and halves the memory and compute requirements.
For 6 nonpalindromic enzymes we use reverse complement from expanded `enzymes.csv`.

## Aho-Corasick multi-pattern search
All enzyme recognition sequences are loaded as patterns into a single Aho-Corasick
automaton. This means the genome is scanned **once** regardless of how many enzymes
are in the run — scanning 23 enzymes costs the same as scanning 1.

`find_overlapping_iter` is used rather than `find_iter` because two different enzyme
recognition sites can start at the same position — overlapping matches must all be
reported.

## Threading model: search threads + one dedicated writer thread
The binary uses `N-1` Rayon threads for parallel genome search (one contig per thread)
and one dedicated OS thread for writing results to disk.

**Dedicated writer thread:**

- Multiple search threads writing to the same file simultaneously would require locking
  
- The writer owns all file handles exclusively — no synchronization needed on the write side
  
- Hits are sent via `mpsc::channel` (multi-producer single-consumer) which is lock-free
  on the send side

## Buffered writing with flush threshold
The writer accumulates rows in a per-enzyme in-memory buffer and flushes to disk
when the buffer exceeds 64 MB. This avoids a syscall per row (which would be the
bottleneck on large genomes) while keeping memory usage bounded.

## Output files (`results/cuts/`)
One CSV per enzyme: `<enzyme_name>.csv`

| Column | Description |
|---|---|
| `accession` | Contig/chromosome name from FASTA header (first whitespace-delimited token) |
| `motif_start` | 0-based start position of the recognition site on the forward strand |
| `cut_position` | Absolute position of the cut, computed as `motif_start + cut_offset` |

## Post-write integrity check
After all hits are written, the writer thread counts newlines in each output file
and compares against the number of hits sent. A mismatch indicates a write error
or filesystem issue and is reported explicitly. This check exists to prevent silent
data loss.

---

# cut_merge (Rust binary)

## Summary
Merges per-enzyme cut site CSVs (output of `rendogbs_finder`) into per-combination
`cuts.csv` files. Each output file contains the cut sites of both enzymes for one
combination, sorted by accession and position, ready for fragment generation.


## Arguments
| Argument | Required | Description |
|---|---|---|
| `--cuts-dir` | yes | Directory containing per-enzyme CSVs from `rendogbs_finder` (`results/cuts/`) |
| `--combinations` | yes | `combinations.csv` produced by `combinations_processing.py` |
| `--out-dir` | yes | Root output directory — one subdir per combination created |
| `--parallel` / `-p` | no (default: 2) | Number of Rayon threads for parallel combination processing |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Enzyme cache: load once, use for all combinations
Binary builds a `cuts_cache` upfront — a `HashMap<enzyme_name, Vec<Cut>>` that loads
each unique enzyme file exactly once. This matters when the same enzyme appears
in many combinations (e.g. `MseI` appears in 8 of the 23 fast combinations) —
without the cache, `MseI` would be read from disk 8 times.

The cache is built before the parallel section and shared read-only across all
worker threads — no locking needed.

## Per-combination processing
For each combination:
1. Both enzyme cut vectors are retrieved from the cache
2. Merged into a single vector of references (`&Cut`) — no data is copied
3. Sorted by `(accession, position)` — this sort order is required by
   `fragment_generation` which assumes sorted input for its vectorized block processing
4. Written to `<out_dir>/<EnzA>_<EnzB>/cuts.csv`

Combinations are processed in parallel via `par_iter()` — each combination is
independent so no synchronization is needed between workers.

## Output files (`results/<EnzA_EnzB>/cuts.csv`)
| Column | Description |
|---|---|
| `accession` | Contig/chromosome name |
| `cut_position` | Cut position on the forward strand |
| `enzyme` | Enzyme name |

## BufWriter
Each `cuts.csv` is written row by row. `BufWriter` accumulates rows 
in an 8 KB kernel buffer (default value) and flushes in larger chunks.

---

# fragment_generation (Rust binary)

## Summary
Generates all adjacent different-enzyme fragment pairs from merged cut sites for
each combination. Output is `fragments.csv` — the complete unfiltered set used
exclusively for **cut statistics and visualisation**. No library selection logic
is applied here, both ends must be different enzymes.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--out-dir` | yes | Root results directory — reads `<combo>/cuts.csv`, writes `<combo>/fragments.csv` |
| `--combinations` | yes | `combinations.csv` from `combinations_processing.py` |
| `--parallel` / `-p` | no (default: 2) | Number of Rayon threads |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Fragment generation logic (find_fragments)
The input `cuts.csv` is already sorted by `(accession, position)` — guaranteed by
`cut_merge`. Fragment generation is a single pass with `windows(2)`:

For each adjacent pair of cuts:
- **Different accession** → skip (cuts on different contigs)
- **Same enzyme** → skip
- **Different enzyme, same accession** → emit fragment

Brute-force enumeration of all cut pairs.

## Arc for the cuts cache
The cuts cache (`HashMap<combo_key, Vec<Cut>>`) is built once in the main thread
and shared read-only across all Rayon worker threads via `Arc`. No locking needed
because no thread mutates the cache.

## Output format (results/<EnzA_EnzB>/fragments.csv)
| Column | Description |
|---|---|
| `accession` | Contig/chromosome name |
| `start_pos` | Cut position of the left enzyme |
| `start_enzyme` | Name of the enzyme at the left end |
| `end_enzyme` | Name of the enzyme at the right end |
| `fragment_length` | `end_pos - start_pos` in base pairs |

## Buffered streaming output
Fragments are written directly to disk using `BufWriter`.
Rows are streamed as they are generated, avoiding allocation of an intermediate
fragment vector or a large in-memory string. `BufWriter` accumulates writes in
memory and flushes them in large blocks, substantially reducing the number of
system calls compared with writing each row individually.

## Output verification
After writing is complete, the output buffer is explicitly flushed and the
generated file is read back once. The number of newline characters is compared 
with the number of fragment rows written during generation. If the counts differ, 
the program terminates with an error, preventing incomplete or corrupted output 
files from silently propagating through the pipeline. Similar to logic 
applied in `rendogbs_finder`.

---

# fragment_generation.py 

## Summary
Core library prediction script. Takes merged cut sites for each combination and
produces the predicted ddRAD/GBS library — fragments that would be
selectively enriched given the user-specified size window. Mainly overlap detection, 
cluster/hotspot handling, uncertain cut tracking, and size filtering.

## Our definition of clusters

**1) Hotspot** - defined as overlap of 2 cut sites 

**Simplified scheme of cut sites**

A ------- B/A ------- B 

A ------- A/B ------- B

**1) Clusters** - defined as overlap of multiple cut sites

**Simplified scheme of cut sites**

A ------- ABABABABABA ------- B 

A ------- AAAAAABBBBB ------- B

> **Note:** Given logic used for emiting possible fragments derived from these
overlaping recognition sites, we use the same logic for hotspots and clusters.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory - paths derived |
| `--size` | yes | Size window in `LOW-HIGH` format, e.g. `150-350` (inclusive both ends) |
| `--parallel` | no (default: 2) | Number of worker processes |

> **Note:** All required arguments are passed automatically by the shell wrapper.


## Worker initialisation: enzyme dict loaded once per process
`ProcessPoolExecutor` spawns N worker processes. Each process calls `init_worker`
exactly once at startup, loading `enzymes.csv` into the global `_GLOBAL_EDICT`.

## Size window parsing (parse_size_range)
Validated early in `main()`. If the format is wrong the pipeline exits
with a clear message.

## Enzyme motif intervals (compute_pair_intervals)
For overlap detection it's needed the full extent of each recognition
site on the genome, not just the cut position. For each cut site:

```
motif_start = cut_position - cut_offset
motif_end   = cut_position + (motif_length - cut_offset)
```

This reconstructs the recognition site interval from the cut position and the
enzyme dictionary. The computation is fully vectorized over all cuts at once using
NumPy boolean masks per enzyme — not loop over individual rows.

## Per-accession block processing (process_accession_block)
The sorted cut list is split into contiguous blocks of the same accession. Each
block is processed independently — fragments cannot span contig boundaries.

### Overlap detection
Two adjacent recognition site intervals overlap if they share at least 1 bp
(touching edges are **not** overlap):

```
overlap = max(s1, s2) < min(e1, e2)
```

A cluster/hotspot is a run of consecutive overlapping pairs. Detected via
run-length encoding on the boolean overlap array — vectorized, not loop.

### Why clusters are not discarded
At a cluster of overlapping recognition sites it is ambiguous which enzyme cuts
at which position. However, discarding the entire cluster would silently lose
fragments that may be real library members. Instead:

**1)** All cuts inside a cluster are flagged as `uncertain_cuts`

**2)** The boundary of the cluster is examined for fragments that pass the size filter

**3)** Up to 2 outermost cuts per side of the cluster boundary are tried

**4)** The first candidate that passes the size filter is written, then the search stops

We present results of clusters via `statistics_cutting.csv`. And results are also
included in final `summary.tsv`.

### Normal fragments (outside clusters)
For cuts outside clusters: adjacent pairs of different enzymes on the same 
accession that pass the size filter are emitted. Fully vectorized via
boolean mask.

## Output files
| File | Description |
|---|---|
| `filtered.csv` | Predicted library fragments passing the size window |
| `uncertain_cuts.csv` | All cut sites inside overlap clusters |
| `statistics_uncertain.csv` | Per-combination counts: total cuts, uncertain cuts, % uncertain, filtered fragments, cluster count |
| `statistics_cutting.csv` | Same content — duplicate written for compatibility with downstream scripts |
| `batch_summary.csv` | One row per combination with all counts — written to `results/` root |

## Path resolution in main()
All paths are derived from `--workdir`:
- `enzymes.csv` — located relative to the script file itself
- `combinations.csv` — `<workdir>/results/combinations.csv`
- `cuts.csv` per combination — `<workdir>/results/<EnzA_EnzB>/cuts.csv`
- All outputs — `<workdir>/results/<EnzA_EnzB>/`

---

# postprocess_metrics.py

## Summary
Computes per-combination metrics that require the reference genome sequence:
fragment length distribution and GC content statistics. Takes both the unfiltered
`fragments.csv` (from Rust binary, for distribution counts) and the library-selected
`filtered.csv` (from `fragment_generation.py`, for GC calculation) as input.

Also produces `contig_lengths.txt` — a genome-wide contig length table sorted
longest-first, used by `plots.py` for chromosome-level visualisation.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory |
| `--ref` | yes | Reference FASTA (plain or gzipped) |
| `--size` | yes | Size window e.g. `150-350` — used to annotate distribution bins |
| `--parallel` | no (default: 2) | Number of worker processes |


## Reference genome indexing
The reference genome is accessed using `pyfaidx`. If a FASTA index (`.fai`) does not already 
exist inside `results/`, it is created automatically before processing begins.

Storing the index inside the pipeline output directory avoids modifying the
original reference genome location while allowing fast random access to
individual fragments


## Reference loading (init_worker)
Each worker process opens its own `pyfaidx.Fasta` object during startup via
`ProcessPoolExecutor(initializer=init_worker)`.

The reference object is stored in the module-level `_WORKER_REF` variable,
allowing every combination processed by that worker to reuse the same indexed
reference without reopening the FASTA file.

Because `pyfaidx` performs indexed random access, only the requested sequence
regions are read from disk rather than loading the entire genome into memory.

## contig_lengths.txt
After the FASTA index has been created, the main process writes
`results/contig_lengths.txt`.

Contig names and lengths are obtained directly from the indexed reference via
`pyfaidx`. Contigs are sorted by decreasing length and written.

`plots.py` uses this table to determine which contigs should be treated as
chromosomes (argument `--chroms N` selects the N longest contigs), providing
approximation of chromosome order for well-assembled genomes.

## Parallel processing
Each enzyme combination is processed independently using `ProcessPoolExecutor`.

**Every worker:**

- loads `fragments.csv`

- loads `filtered.csv`

- computes fragment size distribution

- computes GC statistics

- writes `distribution.csv`

- writes `gc_metrics.csv`

Because combinations are completely independent, processing scales efficiently.

## Distribution (build_distribution_rows)
Fragment lengths are grouped into standard 100 bp bins
(0–99, 100–199, ..., ≥1000 bp) using the **unfiltered**
`fragments.csv` generated by the Rust binary.

These bins represent the complete fragment population before library size
selection. An additional row named `custom_LOW-HIGH` is appended to
`distribution.csv`. Its count corresponds to the number of fragments retained
after size selection (`filtered.csv`).

Percentages are calculated relative to the total number of unfiltered
fragments, allowing direct comparison between the selected library and the
complete fragment population. Standard bins overlapping the user-selected 
size range are marked with the `selected` label.


## GC content (compute_gc_stats)
GC statistics are calculated only from fragments contained in `filtered.csv`.
For every filtered fragment, the corresponding sequence is retrieved directly
from the indexed reference genome using:

`reference[accession][start : start + fragment_length]`

GC content is calculated as:

`(G + C) / fragment_length`

The following statistics are computed using NumPy:

- mean

- median

- minimum

- maximum

- standard deviation


Results are reported as percentages. If no filtered fragments are available, 
all GC statistics are reported as `n/a`.

## Output files
| File | Description |
|---|---|
| `results/contig_lengths.txt` | Tab-separated accession + length, sorted longest first |
| `results/<EnzA_EnzB>/distribution.csv` | Fragment counts per 100 bp bin + custom window row |
| `results/<EnzA_EnzB>/gc_metrics.csv` | GC statistics for filtered fragments: mean, median, min, max, std |

---

# annotation.py

## Summary
Optional step that intersects the predicted library (`filtered.csv`) with genomic
annotations (GFF/GFF3/GTF and/or RepeatMasker TE output) and reports how many
bases of the predicted library overlap each annotation category. Output is
`annotation_summary.csv` per combination.

This script is only invoked by wrapper script if the user provides `--te` or 
`--annotation` (or both). If neither is given the pipeline skips this step entirely. 
Bedtools is required dependency only for annotation.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory |
| `--annotation` | optional | GFF3/GFF/GTF gene annotation (plain or gzipped) |
| `--te` | optional | RepeatMasker `.out` TE annotation (plain or gzipped) |
| `--parallel` | no (default: 4) | Number of combinations processed in parallel |

> **Note:** At least one of `--annotation` or `--te` must be provided, if is wrapper 
script will start `annotation.py`.

## Bedtools
Interval arithmetic on genomic coordinates — intersect, subtract, sort, merge. 
It's in image, not a dependency that user has to download.

## Annotation cache (_annotation_cache/)
GFF and TE files are parsed and converted to sorted, merged per-category BED files
**once** and stored in `results/_annotation_cache/`. All combinations then work
against these cached BED files rather than re-parsing the raw annotation files
per combination.

## GFF processing (build_gff_sources)
Each GFF feature type becomes a separate annotation category (`gff:CDS`,
`gff:exon`, etc.). `region` and `chromosome` entries are skipped as they are
coordinate system entries.

Categories are ordered by `GFF_PRIORITY`:
```
CDS → five_prime_UTR → three_prime_UTR → exon → mRNA → gene
```
This order matters for the exclusive base counting. Any feature types
not in the priority list are appended alphabetically after the priority features.

GFF coordinates are 1-based inclusive — converted to 0-based half-open BED
coordinates by `start - 1`.

## TE processing (build_te_sources)
RepeatMasker `.out` format — first 3 lines are header and skipped. Each TE entry
becomes a category labelled `te:Class_Family` (e.g. `te:LINE_L1`). Coordinates
are already 1-based — converted to 0-based by `start - 1`.

## Exclusive base counting (count_unique_category_bases)
Each genomic base is counted in **at most one** category, 
assigned to the highest-priority category that covers it.

1. Pre-intersect every category with the library BED — work only on
   library-sized intervals from this point forward
2. Walk categories in priority order
3. For each category subtract already-claimed bases (`bedtools subtract`)
4. Count remaining (exclusive) bases for this category
5. Add exclusive bases to the claimed set
6. Move to the next category

This means if a base is covered by both `CDS` and `exon`, it is counted only
under `CDS` (higher priority). `exon` gets credit only for bases not already
claimed by `CDS`. The total across all categories therefore never exceeds
`total_filtered_bases`.

## Subprocess pipeline management
All bedtools calls are run as subprocesses. Chained operations (sort | merge,
intersect | sort | merge) are wired as Unix pipes using `subprocess.Popen` with
`stdout=subprocess.PIPE` — no intermediate temp files for the pipe stages.

`_check_procs` waits for all processes in a pipeline and raises if any exited
non-zero.

## Temporary directories per combination
Each combination's bedtools work happens inside a `tempfile.TemporaryDirectory`
that is deleted automatically when the combination finishes. The annotation cache
(shared across combinations) is kept in `results/_annotation_cache/` and is
presented after run finishes.

##  Worker processes
Combinations are processed in parallel via `ProcessPoolExecutor`. `sources` 
(parsed annotation BED files) and `bedtools` path are loaded once per worker 
process at startup via `_init_annotation_worker` initializer.

## Output (results/<EnzA_EnzB>/annotation_summary.csv)
| Column | Description |
|---|---|
| `category` | Annotation label (`gff:CDS`, `te:LINE_L1`, `total_filtered_bases`) |
| `bases_overlapping` | Number of library bases exclusively in this category |
| `pct_of_library` | Percentage of total library bases |

First data row is always `total_filtered_bases` — the denominator for all
percentages, representing the total non-redundant base count of the predicted
library after `bedtools merge`.

---

# summary.py

## Summary
Combines per-combination CSV outputs from all previous pipeline steps into a
single `summary.tsv`. One row per combination, one column per metric. Intended
for direct import into Excel for comparison across combinations.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory. All paths derived from this. |

## Dynamic column accumulation
Each combination may have a different set of annotation columns depending on
what GFF features or TE families were found in its library. `all_cols` grows as
each combination is processed — any key not yet seen is appended to the column
list. Combinations processed before a new column appears get an empty string for
that column via `row.get(col, "")`.

## Input files per combination (all optional)
| File | Loader | What it contributes |
|---|---|---|
| `statistics_cutting.csv` | `load_stats` | First row only: total_cuts, uncertain_cuts, filtered_fragments, etc. |
| `gc_metrics.csv` | `load_gc` | `metric → value` pairs: gc_mean_pct, gc_median_pct, .... |
| `distribution.csv` | `load_distribution` | `length_range → count` pairs: one column per size bin |
| `annotation_summary.csv` | `load_annotation` | `category → pct_of_library` pairs — absent if annotation.py did not run |

All loaders return `{}` if the file does not exist — missing files produce empty
cells in the output.

---

# plots.py

## Summary
Final step of the pipeline. Reads per-combination CSVs and produces six
publication-ready figures comparing all combinations at once.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory |
| `--size` | yes | Size window used for the run, e.g. `150-350` |
| `--chroms` | yes | How many contigs to treat as chromosomes (longest first from `contig_lengths.txt`) |
| `--dpi` | no (default: 300) | Output resolution |
| `--parallel` | no (default: 2) | Number of worker processes used for parallel plot generation (maximum 5) |

## Combination discovery (discover_combos)
Combinations are discovered from subdirectories of `results/` that contain
`fragments.csv`. This means the plot script is self-contained — it does not
depend on `combinations.csv` and will pick up any combination directory that
exists, including merging multiple runs into one figure - the usage of this 
script is done so advanced users can generate graphs independently from the 
rest of the pipeline. Refer to argument guide, the script will help you to
get the right arguments, but you have to have installed all dependencies 
(numpy and matplotlib).

## Parallel
The five independent plotting functions are executed concurrently using
`ProcessPoolExecutor`. Each worker generates one complete figure independently, 
allowing CPU-intensive tasks (CSV loading, calculations and rendering) 
to run simultaneously. The number of worker processes is controlled by 
`--parallel`, but is internally limited to five because only five independent 
plot-generation tasks exist. If any plotting process fails, the exception 
is propagated immediately and the pipeline terminates with an error.

## Data loading
Large CSV files are read using `numpy.loadtxt()` whenever only a single column
is required, it reduces memory usage and significantly improves loading speed 
for large fragment datasets.

Specifically:

- `fragments.csv` → reads only `fragment_length`
- `filtered.csv` → reads only `fragment_length` or `accession`
- loaded arrays are processed directly with NumPy functions such as
  `histogram()`, `minimum()` and vectorised operations.

Metadata tables (`distribution.csv`, `gc_metrics.csv`,
`annotation_summary.csv`) use `csv.DictReader` because they contain
only a small number of rows and use another layer of logic for separation of TE 
and gene annotation, detailed description is in each graph function bellow.


## Output figures (results/plots/)

### 1. heatmap_fragment_lengths.png
Source: `fragments.csv` (all fragments, unfiltered — from Rust binary)
Shows the full fragment length distribution for every combination as a heatmap
with 10 bp bins up to 1000 bp and log-scale colour. Uses all fragments 
so the user can see the complete picture and evaluate whether the chosen
size window captures the main peak or is offset. Overlaping cut sites end up in 0-9 bp
bin and should be ignored if for specific combination exist large number of clusters
and overlapping motifs, refer to `summary.tsv`.

### 2. heatmap_chrom_distribution.png
Source: `filtered.csv` (library-selected fragments)
Per-chromosome fragment counts for the N longest contigs (user-specified via
`--chroms`). Contig order follows `contig_lengths.txt` — longest contig first.
Highlights whether any combination produces uneven chromosomal coverage,
which is even expected in lot of cases.

### 3. bar_chrom_distribution.png
Source: `filtered.csv`
Same data as the heatmap but as a grouped bar chart — to compare exact
counts between combinations on specific chromosomes.

### 4. gc_distribution.png
Source: `gc_metrics.csv`
Mean GC ± 1 SD per combination with whiskers showing min/max. Allows rapid
visual comparison of GC bias across combinations.

### 5. size_distributions.png (line plot)
Source: `distribution.csv`
Fragment count per 100 bp bin for each combination as overlapping line plots.
The user-specified size window is highlighted.

## Annotation plot (annotation_coverage.png)
Source: `annotation_summary.csv` (absent if annotation.py did not run)
Stacked horizontal bar chart showing what fraction of each combination's library
overlaps GFF features and TE categories. Two panels side by side — GFF left,
TE right — rendered only if the respective data exists.

`total_filtered_bases` and `gff:region` are excluded from the plot
(`_SKIP_CATS`) — these are denominator/coordinate entries, not biological
feature categories.

Categories are split by prefix: `gff:*` uses `tab20` colormap, `te:*` uses
`Set2` — visually distinct palettes.

## Figure sizing: dynamic not fixed
All figures scale their dimensions based on the number of combinations and
chromosomes.

## Agg backend (matplotlib.use("Agg"))
The script runs on servers and in containers that have no display. `Agg` is a
non-interactive backend that renders directly to file without requiring a display
server.
