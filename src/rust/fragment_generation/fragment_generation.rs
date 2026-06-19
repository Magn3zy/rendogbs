// fragment_generation.rs
// Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
//
// Fragment generation for statistics only, one combo per thread
// Output: enzyme specific fragments.csv

use std::{
    collections::HashMap,
    fs,
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
    enzyme:    String,
}

#[derive(Clone, Debug)]
struct Fragment {
    accession:       String,
    start_pos:       i64,
    start_enzyme:    String,
    end_enzyme:      String,
    fragment_length: i64,
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
                accession: r[0].to_string(),
                pos:       r[1].parse().expect("Invalid cut_position"),
                enzyme:    r[2].to_string(),
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

fn find_fragments(cuts: &[Cut]) -> Vec<Fragment> {
    let mut fragments = Vec::new();

    for window in cuts.windows(2) {
        let a = &window[0];
        let b = &window[1];

        if a.accession != b.accession {
            continue;
        }
        if a.enzyme == b.enzyme {
            continue;
        }

        fragments.push(Fragment {
            accession:       a.accession.clone(),
            start_pos:       a.pos,
            start_enzyme:    a.enzyme.clone(),
            end_enzyme:      b.enzyme.clone(),
            fragment_length: b.pos - a.pos,
        });
    }

    fragments
}

fn write_fragments_csv(path: &Path, frags: &[Fragment]) {
    let mut out = String::with_capacity(frags.len() * 64);
    out.push_str("accession,start_pos,start_enzyme,end_enzyme,fragment_length\n");
    for f in frags {
        out.push_str(&format!(
            "{},{},{},{},{}\n",
            f.accession, f.start_pos, f.start_enzyme, f.end_enzyme, f.fragment_length
        ));
    }
    fs::write(path, &out).unwrap_or_else(|_| panic!("Cannot write {:?}", path));
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

    let fragments = find_fragments(cuts);

    write_fragments_csv(&out_dir.join(&key).join("fragments.csv"), &fragments);

    println!("[OK] {} -> fragments: {}", key, fragments.len());
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