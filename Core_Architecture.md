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

## Purpose
First step of the rendogbs pipeline. Validates and resolves enzyme combinations
before any computationally expensive steps run. All heavy work (cut site finding,
fragment generation) happens downstream, exits early on any invalid input.

## Separate validation step
Rust binaries downstream (`rendogbs_finder`, `cut_merge`, `fragment_generation`)
expect clean, validated input. Rather than adding validation logic into each binary,
all input resolution happens here once. If something is wrong with the
user's combinations, the pipeline fails immediately with a clear message.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory. `results/` is created here. |
| `--combinations` | no (default: `fast`) | `fast` uses the 23 predefined combinations. `custom` reads from `--combinations-file`. |
| `--combinations-file` | only if `--combinations custom` | Path to user-supplied CSV with two columns: `enzyme_a`, `enzyme_b`. |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Input files
Both files are bundled with the pipeline (in `scripts/`) and copied automatically
by the wrapper — the user never needs to handle them manually:

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
from the literature. They are provided as a convenience for users who do not have
a specific combination in mind. Users who want to explore other combinations use
`--combinations custom` we support 14 040 combinations.

## Output files (written to `<workdir>/results/`)
- `enzymes_run.csv` — subset of `enzymes.csv` containing only the unique enzymes
  needed for this run. Passed to `rendogbs_finder` (Rust) so it does not load the
  full database unnecessarily.
- `combinations.csv` — validated pairs in `enzyme_a,enzyme_b` format. Used by
  all downstream scripts and binaries as the authoritative list of what to process.

---

# rendogbs_finder (Rust binary)

## Purpose
Finds all recognition site cut positions for each enzyme across the entire reference
genome. Output is one CSV per enzyme containing every cut site position. This is the
most computationally intensive step in the pipeline — all design decisions here
prioritize throughput on large genomes.


## Arguments
| Argument | Required | Description |
|---|---|---|
| `--enzymes-run` | yes | Path to `enzymes_run.csv` produced by `combinations_processing.py` |
| `--ref` / `-r` | yes | Reference genome FASTA (plain or gzipped) |
| `--out-dir` / `-o` | yes | Output directory for per-enzyme CSV files (`results/cuts/`) |
| `--parallel` / `-p` | no (default: 2) | Total threads. One is always reserved for the writer thread, rest search. |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Forward strand only
Recognition sites of restriction enzymes used in ddRAD/GBS are palindromic —
the sequence on the reverse strand is the reverse complement of the forward strand
and identical in terms of cut position. Scanning only the forward strand finds
all cut sites without duplication and halves the memory and compute requirements.
This is a deliberate biological simplification, not an oversight. For 6 
nonpalindromic enzymes we use reverse complement from expanded `enzymes.csv`.

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
- Disk I/O is slow and would stall search threads if done inline
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
or filesystem issue and is reported explicitly. This check exists because silent
data loss during buffered I/O would be problematic.

---

# cut_merge (Rust binary)

## Purpose
Merges per-enzyme cut site CSVs (output of `rendogbs_finder`) into per-combination
`cuts.csv` files. Each output file contains the cut sites of both enzymes for one
combination, sorted by accession and position, ready for fragment generation.


## Arguments
| Argument | Required | Description |
|---|---|---|
| `--cuts-dir` | yes | Directory containing per-enzyme CSVs from `rendogbs_finder` (`results/cuts/`) |
| `--combinations` | yes | `combinations.csv` produced by `combinations_processing.py` |
| `--out-dir` | yes | Root output directory — one subdir per combination is created here |
| `--parallel` / `-p` | no (default: 2) | Number of Rayon threads for parallel combination processing |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Enzyme cache: load once, use for all combinations
Rather than reading each enzyme CSV once per combination that uses it, the binary
builds a `cuts_cache` upfront — a `HashMap<enzyme_name, Vec<Cut>>` that loads
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
| `cut_position` | Absolute cut position on the forward strand |
| `enzyme` | Enzyme name — used downstream to distinguish which end of a fragment belongs to which enzyme |

## BufWriter
Each `cuts.csv` is written row by row. Without buffering this would be one syscall
per row — on a large genome with millions of cut sites this would be the bottleneck.
`BufWriter` accumulates rows in an 8 KB kernel buffer (default value) and flushes 
in larger chunks.

---

# fragment_generation (Rust binary)

## Purpose
Generates all adjacent different-enzyme fragment pairs from merged cut sites for
each combination. Output is `fragments.csv` — the complete unfiltered set used
exclusively for **cut statistics and visualisation**. No library selection logic
is applied here — this binary has no knowledge of size windows, clusters, or
adapter compatibility beyond the basic rule that both ends must be different enzymes.

## What this binary is NOT
This binary is not part of the library prediction logic. It does not interact with
`fragment_generation.py`. The two generate different outputs from the same `cuts.csv`:

- `fragment_generation` (Rust) → `fragments.csv` — all fragments, for graphs cut statistics,
  total cut counts and raw fragment length distributions
- `fragment_generation.py` (Python) → `filtered.csv` — library-selected fragments
  after size filtering and cluster handling, for library size prediction

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
- **Different accession** → skip (cuts on different contigs, no fragment between them)
- **Same enzyme** → skip (same adapter on both ends, not relevant for cut statistics)
- **Different enzyme, same accession** → emit fragment

Intentionally no further logic — this is a brute-force enumeration of all cut pairs.

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

## Write strategy: String buffer then single fs::write
All rows are accumulated into a single `String` and written in one `fs::write` call
per combination — avoids repeated syscalls for what is a bulk write.

---

# fragment_generation.py 

## Purpose
Core library prediction script. Takes merged cut sites for each combination and
produces the predicted ddRAD/GBS library — fragments that would actually be
selectively enriched given the user-specified size window. This is where the
biological complexity lives: overlap detection, cluster/hotspot handling, uncertain
cut tracking, and size filtering.

This script is explicitly separate from the Rust `fragment_generation` binary —
the Rust binary produces raw fragment counts for statistics, this script produces
the predicted library with full biological logic applied.

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
| `--workdir` | yes | Root working directory. All paths are derived from this. |
| `--size` | yes | Size window in `LOW-HIGH` format, e.g. `150-350` (inclusive both ends) |
| `--parallel` | no (default: 2) | Number of worker processes |

> **Note:** All required arguments are passed automatically by the shell wrapper.

## Why Python
The cluster/hotspot logic requires conditional branching, boundary candidate
evaluation, and uncertain cut tracking that is significantly easier to implement
and validate correctly in Python with NumPy than in Rust. The vectorized NumPy
operations on per-accession blocks give sufficient performance — the bottleneck
on large genomes is disk I/O, not computation.

## Worker initialisation: enzyme dict loaded once per process
`ProcessPoolExecutor` spawns N worker processes. Each process calls `init_worker`
exactly once at startup, loading `enzymes.csv` into the global `_GLOBAL_EDICT`.
This means the enzyme dictionary is read from disk N times total (once per worker).

## Size window parsing (parse_size_range)
Validated early in `main()` before any worker is spawned. If the format is wrong
the pipeline exits immediately with a clear message.

## Enzyme motif intervals (compute_pair_intervals)
For overlap detection the script needs to know the full extent of each recognition
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
This split also enables correct cluster detection which is local to each contig.

### Overlap detection
Two adjacent recognition site intervals overlap if they share at least 1 bp
(touching edges are **not** overlap):

```
overlap = max(s1, s2) < min(e1, e2)
```

A cluster/hotspot is a run of consecutive overlapping pairs. Detected via
run-length encoding on the boolean overlap array — vectorized, not loop.

### Why clusters are not simply discarded
At a cluster of overlapping recognition sites it is ambiguous which enzyme cuts
at which position. However, discarding the entire cluster would silently lose
fragments that may be real library members. Instead:

**1)** All cuts inside a cluster are flagged as `uncertain_cuts`

**2)** The boundary of the cluster is examined for fragments that pass the size filter

**3)** Up to 2 outermost cuts per side of the cluster boundary are tried

**4)** The first candidate that passes the size filter is written, then the search stops

**5)** This is conservative: one fragment per cluster boundary side maximum

This approach is transparent — the user sees exactly how many uncertain cuts and
clusters exist per combination via `statistics_uncertain.csv`. And results are also
included in final `summary.tsv`.

### Normal fragments (outside clusters)
For cuts outside clusters the logic is simple: adjacent pairs of different enzymes
on the same accession that pass the size filter are emitted. Fully vectorized via
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
- `enzymes.csv` — located relative to the script file itself (bundled with pipeline)
- `combinations.csv` — `<workdir>/results/combinations.csv`
- `cuts.csv` per combination — `<workdir>/results/<EnzA_EnzB>/cuts.csv`
- All outputs — `<workdir>/results/<EnzA_EnzB>/`

## load_combinations: format tolerance
The combinations CSV parser accepts multiple column name conventions
(`enzyme_a/enzyme_b`, `enzyme1/enzyme2`, `combo`, `combination`, or first two
columns). This is defensive — the authoritative format in this pipeline is
`enzyme_a,enzyme_b` from `combinations_processing.py`, but the tolerance means
the script can also be run as standalone with a manually created CSV in any reasonable
format, mainly this feature was used for testing.

---

# postprocess_metrics.py

## Purpose
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

## A separate step from fragment_generation.py
GC content calculation requires extracting sequences from the reference genome for
every filtered fragment. Loading a large reference genome (e.g. wheat ~14 Gb) into
memory in every worker that processes combinations would be prohibitive. This script
loads the reference once per worker process via `init_worker` and then processes all
combinations assigned to that worker against the already-loaded reference.

## Reference genome loading (init_worker)
The reference FASTA is loaded once per worker process at startup via
`ProcessPoolExecutor(initializer=init_worker)`. The result is stored in the
module-level `_WORKER_REF_SEQS` dict — a plain `{accession: sequence}` mapping.
All combinations processed by that worker share the same in-memory reference
without copying.

FASTA parsing handles both plain and gzipped files, strips whitespace from
sequence lines, takes only the first whitespace-delimited token from each header
line as the accession name (consistent with `rendogbs_finder`), and uppercases
all sequences so GC counting does not need case handling.

## contig_lengths.txt
Contigs are sorted longest-first in the main process before workers start.
This file is written once to `results/contig_lengths.txt` and is used by
`plots.py` to determine which contigs to treat as chromosomes (`--chroms N`
takes the N longest). Sorting by length is a proxy for chromosome order in
well-assembled genomes where chromosomes are typically the longest sequences.

## Distribution (build_distribution_rows)
Fragment counts are binned into standard 100 bp bins (0-99, 100-199, ... >=1000)
using `fragments.csv` — the **unfiltered** set from the Rust binary. This gives
the full picture of where cuts fall across all fragment sizes, not just the
selected window.

The user-specified size window is recorded as an additional `custom_LOW-HIGH` row
at the end of `distribution.csv` with the count from `filtered.csv`. Bins that
overlap the size window are marked with `note="selected"` so downstream scripts
and users can identify them without re-parsing the size argument.

The percentage column uses total unfiltered fragments as the denominator — so the
`custom` row percentage shows what fraction of all possible fragments fall inside
the selected window. This is a key number for evaluating whether a combination
will produce enough loci for the intended study.

## GC content (compute_gc_stats)
Computed only on `filtered.csv` — the library-selected fragments. Computing GC
on unfiltered fragments would mix in fragments that will never be sequenced and
would not reflect the actual GC bias of the library.

For each filtered fragment the sequence is extracted directly from the reference:
`ref_seqs[accession][start : start + fragment_length]`. GC content is then
`(G + C) / length`. Statistics reported: mean, median, min, max, standard deviation
— all as percentages.

Fragments with no sequence in the reference (accession not found) are silently
skipped. If no filtered fragments exist the output is `n/a` for all metrics.

## Output files
| File | Description |
|---|---|
| `results/contig_lengths.txt` | Tab-separated accession + length, sorted longest first. Written once per run. |
| `results/<EnzA_EnzB>/distribution.csv` | Fragment counts per 100 bp bin + custom window row |
| `results/<EnzA_EnzB>/gc_metrics.csv` | GC statistics for filtered fragments: mean, median, min, max, std |

---

# annotation.py

## Purpose
Optional step that intersects the predicted library (`filtered.csv`) with genomic
annotations (GFF/GFF3/GTF and/or RepeatMasker TE output) and reports how many
bases of the predicted library overlap each annotation category. Output is
`annotation_summary.csv` per combination.

This script is only invoked by wrapper script if the user provides `--te` or 
`--annotation` (or both). If neither is given the pipeline skips this step entirely. 
Bedtools is not a required dependency for the rest of the pipeline, only for annotation.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory |
| `--annotation` | optional | GFF3/GFF/GTF gene annotation (plain or gzipped) |
| `--te` | optional | RepeatMasker `.out` TE annotation (plain or gzipped) |
| `--parallel` | no (default: 4) | Number of combinations processed in parallel |

> **Note:** At least one of `--annotation` or `--te` must be provided, if is wrapper 
script will start `annotation.py`.

## Bedtools (external dependency)
Interval arithmetic on genomic coordinates — intersect, subtract, sort, merge —
is exactly what bedtools is designed for. Reimplementing this correctly in Python
would be substantial work with no benefit. Bedtools is standard in any genomics
environment and is bundled in the Docker/Singularity containers, so it is never
a user-facing dependency.

## Annotation cache (_annotation_cache/)
GFF and TE files are parsed and converted to sorted, merged per-category BED files
**once** and stored in `results/_annotation_cache/`. All combinations then work
against these cached BED files rather than re-parsing the raw annotation files
per combination. This is useful when there are many combinations — parsing a large
GFF3 once as opposed to once per combination.

## GFF processing (build_gff_sources)
Each GFF feature type becomes a separate annotation category (`gff:CDS`,
`gff:exon`, etc.). `region` and `chromosome` entries are skipped as they are
coordinate system entries, not biological features.

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
The core methodological decision: each genomic base is counted in **at most one**
category, assigned to the highest-priority category that covers it.

The algorithm:
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

This is a deliberate design choice — additive counting (where one base can be
counted in multiple categories) would make the percentages sum to more than 100%
and be misleading for library composition interpretation.

## Subprocess pipeline management
All bedtools calls are run as subprocesses. Chained operations (sort | merge,
intersect | sort | merge) are wired as Unix pipes using `subprocess.Popen` with
`stdout=subprocess.PIPE` — no intermediate temp files for the pipe stages.

`_check_procs` waits for all processes in a pipeline and raises if any exited
non-zero. This is necessary because a non-zero exit in the middle of a pipe
does not automatically propagate — it must be checked explicitly.

## Temporary directories per combination
Each combination's bedtools work happens inside a `tempfile.TemporaryDirectory`
that is deleted automatically when the combination finishes. The annotation cache
(shared across combinations) is kept in `results/_annotation_cache/` and is
persistent after run finishes.

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

## Purpose
Aggregates per-combination CSV outputs from all previous pipeline steps into a
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

This means column order reflects the order combinations were processed and the
order keys appeared, not a fixed schema. The first combination that has
`gff:CDS` determines its position in the header — all subsequent combinations
either have a value for it or get an empty cell.

## Input files per combination (all optional)
| File | Loader | What it contributes |
|---|---|---|
| `statistics_cutting.csv` | `load_stats` | First row only: total_cuts, uncertain_cuts, filtered_fragments, etc. |
| `gc_metrics.csv` | `load_gc` | `metric → value` pairs: gc_mean_pct, gc_median_pct, etc. |
| `distribution.csv` | `load_distribution` | `length_range → count` pairs: one column per size bin |
| `annotation_summary.csv` | `load_annotation` | `category → pct_of_library` pairs — absent if annotation.py did not run |

All loaders return `{}` if the file does not exist — missing files produce empty
cells in the output, not errors. This means summary.py runs correctly regardless
of whether annotation.py was run.

## No parallelism
Each combination reads four small CSV files. The bottleneck is disk I/O on small
files, not computation — parallelism would add overhead without benefit here.

---

# plots.py

## Purpose
Final step of the pipeline. Reads per-combination CSVs and produces six
publication-ready figures comparing all combinations at once. No computation
happens here — only reading already-produced CSVs and rendering plots.
The script is intentionally decoupled from the rest of the pipeline and can
be re-run independently at any time (e.g. with different `--chroms` or `--dpi`)
without re-running any upstream steps.

## Arguments
| Argument | Required | Description |
|---|---|---|
| `--workdir` | yes | Root working directory |
| `--size` | yes | Size window used for the run, e.g. `150-350` — used to highlight the selected window on the size distribution plot |
| `--chroms` | yes | How many contigs to treat as chromosomes in the chromosome distribution plots (longest first from `contig_lengths.txt`) |
| `--dpi` | no (default: 300) | Output resolution — 300 dpi is print quality, lower values for faster preview |

## Combination discovery (discover_combos)
Combinations are discovered from subdirectories of `results/` that contain
`fragments.csv`. This means the plot script is self-contained — it does not
depend on `combinations.csv` and will pick up any combination directory that
exists, including partial runs.

## Output figures (results/plots/)

### 1. heatmap_fragment_lengths.png
Source: `fragments.csv` (all fragments, unfiltered — from Rust binary)
Shows the full fragment length distribution for every combination as a heatmap
with 10 bp bins up to 1000 bp and log-scale colour. Uses all fragments (not just
filtered) so the user can see the complete picture and evaluate whether the chosen
size window captures the main peak or is offset. Overlaping cut sites end up in 0-9 bp
bin and should be ignored.

### 2. heatmap_chrom_distribution.png
Source: `filtered.csv` (library-selected fragments)
Per-chromosome fragment counts for the N longest contigs (user-specified via
`--chroms`). Contig order follows `contig_lengths.txt` — longest contig first,
which typically corresponds to chromosome order in well-assembled genomes.
Highlights whether any combination produces uneven chromosomal coverage,
which is even expected in lot of cases.

### 3. bar_chrom_distribution.png
Source: `filtered.csv`
Same data as the heatmap but as a grouped bar chart — easier to compare exact
counts between combinations on specific chromosomes. Both chromosome plots are
produced from the same data pass.

### 4. gc_distribution.png
Source: `gc_metrics.csv`
Mean GC ± 1 SD per combination with whiskers showing min/max. Allows rapid
visual comparison of GC bias across combinations. Combinations with
`gc_mean_pct == "n/a"` (no filtered fragments) are silently skipped.

### 5. size_distributions.png (line plot)
Source: `distribution.csv`
Fragment count per 100 bp bin for each combination as overlapping line plots.
The user-specified size window is highlighted as a shaded region so the user
can immediately see how much of the distribution falls inside vs. outside the
selected window.

## Annotation plot (annotation_coverage.png)
Source: `annotation_summary.csv` (absent if annotation.py did not run)
Stacked horizontal bar chart showing what fraction of each combination's library
overlaps GFF features and TE categories. Two panels side by side — GFF left,
TE right — rendered only if the respective data exists.

**Exclusive base accounting** is already done in `annotation.py` — the
percentages in `annotation_summary.csv` are already non-overlapping. The plot
simply stacks them. The remainder to 100% is rendered as "unannotated" in grey.

`total_filtered_bases` and `gff:region` are excluded from the plot
(`_SKIP_CATS`) — these are denominator/coordinate entries, not biological
feature categories.

Categories are split by prefix: `gff:*` uses `tab20` colormap, `te:*` uses
`Set2` — visually distinct palettes so GFF and TE panels are immediately
distinguishable.

## Figure sizing: dynamic not fixed
All figures scale their dimensions based on the number of combinations and
chromosomes:
- More combinations → taller heatmaps, more bar groups
- More chromosomes → wider bar/heatmap plots
- Font sizes are also scaled down for large numbers of combinations to prevent
  label overlap

## Agg backend (matplotlib.use("Agg"))
The script runs on servers and in containers that have no display. `Agg` is a
non-interactive backend that renders directly to file without requiring a display
server.

## Global style
Set once via `plt.rcParams` at module level — DejaVu Sans for font (universally
available, no LaTeX dependency), consistent axis label sizes across all plots.
DPI is runtime-configurable via `--dpi` and overrides the module-level default.
