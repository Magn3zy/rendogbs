use std::{
    fs,
    path::{PathBuf, Path},
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
                pos:       r[2].parse().expect("Invalid cut_position"),
                enzyme:    enzyme.to_string(),
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

// cuts.csv - accession, cut_position, enzyme

fn write_cuts_csv(path: &Path, cuts: &[Cut]) {
    let mut out = String::with_capacity(cuts.len() * 48);
    out.push_str("accession,cut_position,enzyme\n");
    for c in cuts {
        out.push_str(&format!("{},{},{}\n", c.accession, c.pos, c.enzyme));
    }
    fs::write(path, &out).unwrap_or_else(|_| panic!("Cannot write {:?}", path));
}

// one combination - READ-ONLY - cuts/<ea>.csv and cuts/<eb>.csv, WRITE - results/<ea>_<eb>/cuts.csv
fn process_combo(
    combo:    &Combination,
    cuts_dir: &Path,
    out_dir:  &Path,
) {
    let ea = &combo.enzyme_a;
    let eb = &combo.enzyme_b;

    // read only
    let path_a = cuts_dir.join(format!("{}.csv", ea));
    let path_b = cuts_dir.join(format!("{}.csv", eb));
    assert!(path_a.exists(), "Missing cuts file: {:?}", path_a);
    assert!(path_b.exists(), "Missing cuts file: {:?}", path_b);

    let cuts = {
        let mut merged = load_cuts(&path_a, ea);
        let mut cuts_b = load_cuts(&path_b, eb);
        merged.append(&mut cuts_b);
        merged.sort_by(|x, y| {
            x.accession.cmp(&y.accession)
                .then_with(|| x.pos.cmp(&y.pos))
        });
        merged
    };

    // write cuts.csv into new dir
    let combo_dir = out_dir.join(format!("{}_{}", ea, eb));
    fs::create_dir_all(&combo_dir)
        .unwrap_or_else(|_| panic!("Cannot create dir: {:?}", combo_dir));

    write_cuts_csv(&combo_dir.join("cuts.csv"), &cuts);

    println!(
        "[OK] {}_{} -> cuts: {}",
        ea, eb, cuts.len()
    );

    drop(cuts);
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

    let cuts_dir = Arc::new(cli.cuts_dir);
    let out_dir  = Arc::new(cli.out_dir);

    combos.par_iter().for_each(|combo| {
        process_combo(combo, &cuts_dir, &out_dir);
    });

    println!("[OK] done.");
}
