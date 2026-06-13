use std::{
    collections::HashMap,
    fs,
    io::Write,
    path::{PathBuf, Path},
    thread,
}

use clap::Parser;
use csv::ReaderBuilder;
use rayon::prelude::*;

struct Cut {
    accession: String,
    pos: i64,
    enzyme: String,
}

struct Fragment {
    accession: String,
    start_pos: u64,
    start_enzyme: String,
    end_enzyme: String,
    fragment_lenght: u64,
}

#[derive(Parser)]
#[command(name = "cut_adjuscent")]
struct Cli {
    #[arg(long)]
    combinations: PathBuf,

    #[arg(long)]
    cuts: PathBuf,

    #[arg(long)]
    out_dir: PathBuf,
    
    #[arg(long, short = 'p', default_value_t = 2)]
    parallel: usize,
}

fn load_cuts(path: &Path, enzyme_name: &str) -> Vec<Cut> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .expect("Cannot open cuts CSV");

    rdr.records()
        .filter_map(|r| {
            let r = r.expect("Invalid CSV row in cuts file");
            if r.len() != 3 {
                return None;
            }
            Cut {
                accession: r[0].to_string(),
                pos: r[1].parse().expect("Invalid cut_position"),
                enzyme: r[2].to_string(),
            }
        })
        .collect()
}

fn load_combinations(path: &Path) -> Vec<Fragment> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .expect("Cannot open combinations CSV");

    rdr.records()
        .filter_map(|r| {
            let r = r.expect("Invalid CSV row in combinations file");
            if r.len() !=5 {
                return None;
            }
            Fragment {
                accession: r[0].to_string(),
                start_pos: r[1].parse().expect("Invalid start_pos"),
                start_enzyme: r[2].to_string(),
                end_enzyme: r[3].to_string(),
                fragment_length: r[4].parse().expect("Invalid fragment_length"),
            }
        })
        .collect()
}

fn group_by_contig(fragment: Vec<Cut>) -> HashMap<String, Vec<Cut>> {
    let mut map: HashMap<String, Vec<Cut>> = HashMap::new();

    for c in fragment {
        map.entry(c.accession.clone()).or_default().push(c);
    }

    map
}

fn main() {

}
