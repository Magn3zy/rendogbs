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

#[derive(Clone, Debug)]
struct Fragment {
    accession:       String,
    start_pos:       i64,
    start_enzyme:    String,
    end_enzyme:      String,
    fragment_length: i64,
}

#[derive(Clone, Debug)]
struct SizeRange {
    low:  i64,
    high: i64,
}

fn parse_size_range(s: &str) -> Result<SizeRange, String> {
    let (low_str, high_str) = s
        .split_once('-')
        .ok_or_else(|| format!("--size must be LOW-HIGH, e.g. 150-350, you wrote: '{}'", s))?;

    let low: i64 = low_str
        .parse()
        .map_err(|_| format!("Invalid number: {}", low_str))?;

    let high: i64 = high_str
        .parse()
        .map_err(|_| format!("Invalid number: {}", high_str))?;

    if low > high {
        return Err(format!(
            "LOW ({}) must be lower than HIGH ({})",
            low, high
        ));
    }

    Ok(SizeRange { low, high })
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

    #[arg(long, value_parser = parse_size_range)]
    size: SizeRange,
}

// load cuts_clean.csv (accession, cut_position, enzyme)
fn load_cuts_clean(path: &Path) -> Vec<Cut> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .unwrap_or_else(|_| panic!("Cannot open cuts_clean CSV: {:?}", path));

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
            accession:    a.accession.clone(),
            start_pos:    a.pos,
            start_enzyme: a.enzyme.clone(),
            end_enzyme:   b.enzyme.clone(),
            fragment_length: b.pos - a.pos,
        });
    }

    fragments
}

// write fragments.csv - one thread per combination
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

// one combo processing: READ-ONLY - results/<ea>_<eb>/cuts_clean.csv, WRITE - results/<ea>_<eb>/fragments.csv/filtered.csv
fn process_combo(
    combo:   &Combination,
    out_dir: &Path,
    size:    &SizeRange,
) {
    let ea = &combo.enzyme_a;
    let eb = &combo.enzyme_b;

    let combo_dir = out_dir.join(format!("{}_{}", ea, eb));
    let cuts_clean_path = combo_dir.join("cuts.csv");

    assert!(
        cuts_clean_path.exists(),
        "Missing cuts_clean.csv: {:?}",
        cuts_clean_path
    );

    // read cuts_clean.csv
    let cuts = load_cuts_clean(&cuts_clean_path);

    // lenght of all fragments
    let all = find_fragments(&cuts);
    let all_count = all.len();

    // filter by size range
    let filtered: Vec<Fragment> = all
        .iter()
        .filter(|f| f.fragment_length >= size.low && f.fragment_length <= size.high)
        .cloned()
        .collect();
    let filtered_count = filtered.len();

    // write csv fragments.csv a filtered.csv into results/<ea>_<eb>/ per thread
    write_fragments_csv(&combo_dir.join("fragments.csv"), &all);
    write_fragments_csv(&combo_dir.join("filtered.csv"), &filtered);

    println!(
        "[OK] {}_{} -> all: {}, filtered: {}",
        ea, eb, all_count, filtered_count
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
    println!("[INFO] size window:  {}-{}", cli.size.low, cli.size.high);

    let out_dir = Arc::new(cli.out_dir);
    let size    = Arc::new(cli.size);

    combos.par_iter().for_each(|combo| {
        process_combo(combo, &out_dir, &size);
    });

    println!("[OK] done.");
}
