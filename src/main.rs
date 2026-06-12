// Input:  --enzymes-run  enzymes_run.csv
//         --ref    genome.fa / .fa.gz
//         --out-dir      output directory
//         --parallel      N contigs processed in parallel (writers use 1, searchers use parallel-1)
// Output: cuts/EcoRI.csv, cuts/MseI.csv ...
//        header: accession,motif_start,cut_position

use std::{
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::mpsc,
    thread,
};

use aho_corasick::AhoCorasick;
use clap::Parser;
use csv::ReaderBuilder;
use rayon::prelude::*;

#[derive(Debug, Clone)]
struct Enzyme {
    name: String,
    sequence: String,
    cut_offset: i64,
}

struct Hit {
    enzyme_idx:  u32,
    acc_idx:     u32,
    motif_start: u64,
    cut_position: i64,
}

#[derive(Parser)]
#[command(name = "rendogbs_finder")]
struct Cli {
    #[arg(long)]
    enzymes_run: PathBuf,

    #[arg(long, short = 'r')]
    r#ref: PathBuf,

    #[arg(long, short = 'o')]
    out_dir: PathBuf,

    #[arg(long, short = 'p', default_value_t = 2)]
    parallel: usize,
}

fn load_enzymes(path: &Path) -> Vec<Enzyme> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .expect("Cannot open enzymes_run.csv");

    rdr.records()
        .map(|rec| {
            let rec = rec.expect("Invalid CSV row");
            Enzyme {
                name: rec[0].to_string(),
                sequence: rec[1].to_string(),
                cut_offset: rec[2].parse().expect("Invalid cut_offset"),
            }
        })
        .collect()
}

fn load_contigs(path: &Path) -> Vec<(String, Vec<u8>)> {
    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = needletail::parse_fastx_file(path).expect("Cannot read FASTA");

    while let Some(rec) = reader.next() {
        let rec = rec.expect("FASTA parse error");
        let acc = std::str::from_utf8(rec.id())
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap_or("?")
            .to_string();
        contigs.push((acc, rec.seq().to_vec()));
    }
    contigs
}

fn main() {
    let cli = Cli::parse();

    // at least 1 writer + N-1 searchers
    let search_threads = cli.parallel.saturating_sub(1).max(1);
    rayon::ThreadPoolBuilder::new()
        .num_threads(search_threads)
        .build_global()
        .unwrap();

    fs::create_dir_all(&cli.out_dir).unwrap();

    let enzymes = load_enzymes(&cli.enzymes_run);
    eprintln!("[INFO] {} patterns loaded", enzymes.len());

    let patterns: Vec<&str> = enzymes.iter().map(|e| e.sequence.as_str()).collect();
    let ac = AhoCorasick::builder()
        .ascii_case_insensitive(true)
        .build(&patterns)
        .expect("AhoCorasick build failed");

    let contigs = load_contigs(&cli.r#ref);
    eprintln!(
        "[INFO] {} contigs | {} search threads + 1 writer thread",
        contigs.len(),
        search_threads
    );

    let (tx, rx) = mpsc::sync_channel::<Hit>(65536);

    // writer thread outside of rayon pool
    let out_dir = cli.out_dir.clone();
    let enzyme_names: Vec<String> = enzymes.iter().map(|e| e.name.clone()).collect();
    let acc_names_owned: Vec<String> = contigs.iter().map(|(acc, _)| acc.clone()).collect();
    let writer_handle = thread::spawn(move || {
        let mut writers: HashMap<u32, BufWriter<fs::File>> = HashMap::new();

        for hit in rx {
            let w = writers.entry(hit.enzyme_idx).or_insert_with(|| {
                let name = &enzyme_names[hit.enzyme_idx as usize];
                let path = out_dir.join(format!("{}.csv", name));
                let mut bw = BufWriter::with_capacity(
                    1 << 20, // 1 MB buffer
                    fs::File::create(&path).unwrap(),
                );
                writeln!(bw, "accession,motif_start,cut_position").unwrap();
                eprintln!("  -> {}.csv (created)", name);
                bw
            });

            writeln!(
                w,
                "{},{},{}",
                acc_names_owned[hit.acc_idx as usize],
                hit.motif_start,
                hit.cut_position,
            )
            .unwrap();
        }

        // channel closed — flush writer
        for (_, mut w) in writers {
            w.flush().unwrap();
        }
    });

    // parallel search
    contigs.par_iter().enumerate().for_each_with(tx, |tx, (acc_idx, (_acc, seq))| {
        for m in ac.find_overlapping_iter(seq) {
            let idx = m.pattern().as_usize();
            let cut = m.start() as i64 + enzymes[idx].cut_offset;

            tx.send(Hit {
                enzyme_idx: idx as u32,
                acc_idx: acc_idx as u32,
                motif_start: m.start() as u64,
                cut_position: cut,
            })
            .expect("Writer thread died unexpectedly");
        }
    });
    // tx dropped here — rx loop inside writer ends

    writer_handle.join().expect("Writer thread panicked");
    eprintln!("[OK] done.");
}