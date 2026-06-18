use std::{
    collections::{HashMap, HashSet},
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
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
    pos: i64,
    enzyme: String,
}

#[derive(Parser)]
#[command(name = "cuts_merge")]
struct Cli {
    #[arg(long)]
    cuts_dir: PathBuf,

    #[arg(long)]
    combinations: PathBuf,

    #[arg(long)]
    out_dir: PathBuf,

    #[arg(long, short = 'p', default_value_t = 2)]
    parallel: usize,
}

fn load_cuts(path: &Path, enzyme: &str) -> Vec<Cut> {
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
                pos: r[2].parse().expect("Invalid cut_position"),
                enzyme: enzyme.to_string(),
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

fn write_cuts_csv(path: &Path, cuts: &[&Cut]) {
    let file = fs::File::create(path)
        .unwrap_or_else(|_| panic!("Cannot write {:?}", path));
    let mut out = BufWriter::new(file);

    writeln!(out, "accession,cut_position,enzyme")
        .expect("Failed to write CSV header");

    for c in cuts {
        writeln!(out, "{},{},{}", c.accession, c.pos, c.enzyme)
            .expect("Failed to write CSV row");
    }
}

fn build_cuts_cache(cuts_dir: &Path, combos: &[Combination]) -> HashMap<String, Vec<Cut>> {
    let mut enzymes = HashSet::new();
    for combo in combos {
        enzymes.insert(combo.enzyme_a.clone());
        enzymes.insert(combo.enzyme_b.clone());
    }

    let mut cache = HashMap::with_capacity(enzymes.len());

    for enzyme in enzymes {
        let path = cuts_dir.join(format!("{}.csv", enzyme));
        assert!(path.exists(), "Missing cuts file: {:?}", path);

        let cuts = load_cuts(&path, &enzyme);
        cache.insert(enzyme, cuts);
    }

    cache
}

fn process_combo(
    combo: &Combination,
    cuts_cache: &HashMap<String, Vec<Cut>>,
    out_dir: &Path,
) {
    let ea = &combo.enzyme_a;
    let eb = &combo.enzyme_b;

    let cuts_a = cuts_cache
        .get(ea)
        .unwrap_or_else(|| panic!("Missing enzyme in cache: {}", ea));
    let cuts_b = cuts_cache
        .get(eb)
        .unwrap_or_else(|| panic!("Missing enzyme in cache: {}", eb));

    let mut merged: Vec<&Cut> = Vec::with_capacity(cuts_a.len() + cuts_b.len());
    merged.extend(cuts_a.iter());
    merged.extend(cuts_b.iter());

    merged.sort_by(|x, y| {
        x.accession
            .cmp(&y.accession)
            .then_with(|| x.pos.cmp(&y.pos))
    });

    let combo_dir = out_dir.join(format!("{}_{}", ea, eb));
    fs::create_dir_all(&combo_dir)
        .unwrap_or_else(|_| panic!("Cannot create dir: {:?}", combo_dir));

    write_cuts_csv(&combo_dir.join("cuts.csv"), &merged);

    println!("[OK] {}_{} -> cuts: {}", ea, eb, merged.len());
}

fn main() {
    let cli = Cli::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.parallel)
        .build_global()
        .expect("Cannot build Rayon thread pool");

    fs::create_dir_all(&cli.out_dir)
        .unwrap_or_else(|_| panic!("Cannot create out_dir: {:?}", cli.out_dir));

    let combos = load_combinations(&cli.combinations);

    println!("[INFO] combinations: {}", combos.len());
    println!("[INFO] threads:      {}", cli.parallel);

    let cuts_cache = build_cuts_cache(&cli.cuts_dir, &combos);

    combos.par_iter().for_each(|combo| {
        process_combo(combo, &cuts_cache, &cli.out_dir);
    });

    println!("[OK] done.");
}