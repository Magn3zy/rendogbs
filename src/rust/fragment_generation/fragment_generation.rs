// fragment_generation.rs
// Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
//
// Fragment generation for statistics only, one combo per thread
// Output: enzyme specific fragments.csv

use std::{
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

use clap::Parser;
use csv::ReaderBuilder;
use rayon::prelude::*;

struct Combination {
    enzyme_a: String,
    enzyme_b: String,
}

#[derive(Clone, Debug)]
struct Cut {
    accession: String,
    pos:       i64,
    cut_position:    i64,
}

#[derive(Parser)]
#[command(name = "fragment_generation")]
struct Cli {
    #[arg(long)]
    out_dir: PathBuf,

    #[arg(long)]
    combinations: PathBuf,

    #[arg(long, short = 'p', default_value_t = 2)]
    parallel: usize,
}

fn load_cuts(path: &Path) -> Vec<Cut> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .unwrap_or_else(|_| panic!("Cannot open cuts CSV: {:?}", path));

    rdr.records()
        .filter_map(|r| {
            let r = r.expect("Invalid CSV row");
            if r.len() < 3 {
                return None;
            }
            Some(Cut {
                accession:    r[0].to_string(),
                pos:          r[1].parse().expect("Invalid motif_start"),
                cut_position: r[2].parse().expect("Invalid cut_position"),
            })
        })
        .collect()
}

fn load_combinations(path: &Path) -> Vec<Combination> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .unwrap_or_else(|_| panic!("Cannot open combinations CSV: {:?}", path));

    rdr.records()
        .filter_map(|r| {
            let r = r.expect("Invalid CSV row");
            if r.len() < 2 {
                return None;
            }
            Some(Combination {
                enzyme_a: r[0].to_string(),
                enzyme_b: r[1].to_string(),
            })
        })
        .collect()
}

fn build_cuts_cache(out_dir: &Path, combos: &[Combination]) -> HashMap<String, Vec<Cut>> {
    let mut cache = HashMap::with_capacity(combos.len());

    for combo in combos {
        let key = format!("{}_{}", combo.enzyme_a, combo.enzyme_b);
        let cuts_path = out_dir.join(&key).join("cuts.csv");
        assert!(cuts_path.exists(), "Missing cuts.csv: {:?}", cuts_path);
        cache.insert(key, load_cuts(&cuts_path));
    }

    cache
}

fn write_fragments_and_verify(path: &Path, cuts: &[Cut]) -> u64 {
    // stream fragments directly to BufWriter — no intermediate Vec<Fragment>
    let file = fs::File::create(path)
        .unwrap_or_else(|_| panic!("Cannot create {:?}", path));
    let mut out = BufWriter::new(file); // flushes in chunks not per row

    writeln!(out, "accession,start_pos,start_enzyme,end_enzyme,fragment_length")
        .expect("Failed to write CSV header");

    let mut written: u64 = 0; // count rows as we write, not after

    for window in cuts.windows(2) {
        let a = &window[0];
        let b = &window[1];

        if a.accession != b.accession {
            continue;
        }
        if a.cut_position == b.cut_position {
            continue;
        }

        writeln!(
            out,
            "{},{},{},{},{}",
            a.accession,
            a.pos,
            a.cut_position,
            b.cut_position,
            b.pos - a.pos,
        )
        .expect("Failed to write CSV row");

        written += 1;
    }

    // explicit flush before verify — ensures all buffered bytes hit the file
    out.flush().expect("Failed to flush BufWriter");

    // integrity check: count newlines in written file vs expected row count
    let content = fs::read(path)
        .unwrap_or_else(|_| panic!("Cannot read back {:?}", path));
    let actual = content.iter().filter(|&&b| b == b'\n').count() as u64 - 1;

    if actual != written {
        eprintln!(
            "[MISMATCH] {:?}: expected {} fragment rows, got {} newlines",
            path, written, actual
        );
        std::process::exit(1); // error 1
    }

    written
}

fn process_combo(
    combo:      &Combination,
    cuts_cache: &HashMap<String, Vec<Cut>>,
    out_dir:    &Path,
) {
    let key = format!("{}_{}", combo.enzyme_a, combo.enzyme_b);

    let cuts = cuts_cache
        .get(&key)
        .unwrap_or_else(|| panic!("Missing combo in cache: {}", key));

    let out_path = out_dir.join(&key).join("fragments.csv");
    let n = write_fragments_and_verify(&out_path, cuts);

    println!("[OK] {} -> fragments: {}", key, n);
}

fn main() {
    let cli = Cli::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.parallel)
        .build_global()
        .expect("Cannot build Rayon thread pool");

    let combos = load_combinations(&cli.combinations);

    println!("[INFO] combinations: {}", combos.len());
    println!("[INFO] threads:      {}", cli.parallel);

    let cuts_cache = Arc::new(build_cuts_cache(&cli.out_dir, &combos));

    combos.par_iter().for_each(|combo| {
        process_combo(combo, &cuts_cache, &cli.out_dir);
    });

    println!("[OK] done.");
}
