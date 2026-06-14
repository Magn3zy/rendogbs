use std::{
    fs,
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

#[derive(Clone, Debug)]
struct Fragment {
    accession: String,
    start_pos: i64,
    start_enzyme: String,
    end_enzyme: String,
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

// load cuts.csv (accession, cut_position, enzyme)
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
                pos: r[1].parse().expect("Invalid cut_position"),
                enzyme: r[2].to_string(),
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

// Cuts order: (accession + pos) — iter: [i, i+1]
fn find_fragments(cuts: &[Cut]) -> Vec<Fragment> {
    let mut fragments: Vec<Fragment> = Vec::new();

    for window in cuts.windows(2) {
        let a = &window[0];
        let b = &window[1];

        // skip if different accession
        if a.accession != b.accession {
            continue;
        }

        // skip if same enzyme
        if a.enzyme == b.enzyme {
            continue;
        }

        fragments.push(Fragment {
            accession: a.accession.clone(),
            start_pos: a.pos,
            start_enzyme: a.enzyme.clone(),
            end_enzyme: b.enzyme.clone(),
            fragment_length: b.pos - a.pos,
        });
    }

    fragments
}

// write fragments.csv
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

// one combo processing: READ - results/<ea>_<eb>/cuts.csv, WRITE - results/<ea>_<eb>/fragments.csv
fn process_combo(combo: &Combination, out_dir: &Path) {
    let ea = &combo.enzyme_a;
    let eb = &combo.enzyme_b;

    let combo_dir = out_dir.join(format!("{}_{}", ea, eb));
    let cuts_path = combo_dir.join("cuts.csv");

    // read cuts.csv
    let cuts = load_cuts(&cuts_path);

    // generate fragments from all cuts
    let fragments = find_fragments(&cuts);

    // write fragments.csv
    write_fragments_csv(&combo_dir.join("fragments.csv"), &fragments);

    println!(
        "[OK] {}_{} -> fragments: {}",
        ea,
        eb,
        fragments.len(),
    );
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

    let out_dir = cli.out_dir;

    combos.par_iter().for_each(|combo| {
        process_combo(combo, &out_dir);
    });

    println!("[OK] done.");
}
