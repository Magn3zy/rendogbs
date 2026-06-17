use std::{
    collections::HashMap,
    fs,
    io::Write,
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
    enzyme_idx:   u32,
    acc_idx:      u32,
    motif_start:  u64,
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
            .trim_end_matches('\r')
            .to_string();
        contigs.push((acc, rec.seq().to_vec()));
    }
    contigs
}

fn main() {
    let cli = Cli::parse();

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

    let out_dir = cli.out_dir.clone();
    let enzyme_names: Vec<String> = enzymes.iter().map(|e| e.name.clone()).collect();
    let acc_names_owned: Vec<String> = contigs.iter().map(|(acc, _)| acc.clone()).collect();

    let writer_handle = thread::spawn(move || {
        // klíč = název enzymu — více IUPAC variant stejného enzymu → jeden soubor
        let mut files:    HashMap<String, fs::File> = HashMap::new();
        let mut buffers:  HashMap<String, Vec<u8>>  = HashMap::new(); // cant be torn appart during flush
        let mut expected: HashMap<String, u64>       = HashMap::new();

        const FLUSH_BYTES: usize = 64 * 1024 * 1024; // 64kB passed test 64 MB for production

        for hit in rx {
            let name = enzyme_names[hit.enzyme_idx as usize].clone();

            if !files.contains_key(&name) {
                let path = out_dir.join(format!("{}.csv", name));
                let mut f = fs::File::create(&path).unwrap();
                f.write_all(b"accession,motif_start,cut_position\n").unwrap();
                eprintln!("  -> {}.csv (created)", name);
                files.insert(name.clone(), f);
                buffers.insert(name.clone(), Vec::with_capacity(FLUSH_BYTES + 256));
                expected.insert(name.clone(), 0);
            }

            let buf = buffers.get_mut(&name).unwrap();
            write!(
                buf,
                "{},{},{}\n",
                acc_names_owned[hit.acc_idx as usize],
                hit.motif_start,
                hit.cut_position,
            ).unwrap();

            *expected.get_mut(&name).unwrap() += 1;

            if buf.len() >= FLUSH_BYTES {
                files.get_mut(&name).unwrap().write_all(buf).unwrap();
                buf.clear();
            }
        }

        // flush last
        for (name, buf) in &buffers {
            if !buf.is_empty() {
                files.get_mut(name).unwrap().write_all(buf).unwrap();
            }
        }

        // check line counts
        eprintln!("[CHECK] verifying line counts...");
        let mut all_ok = true;
        for (name, exp) in &expected {
            let path = out_dir.join(format!("{}.csv", name));
            let content = fs::read(&path).unwrap();
            let actual = content.iter().filter(|&&b| b == b'\n').count() as u64 - 1;
            if actual != *exp {
                eprintln!("[MISMATCH] {}: expected {} rows, got {}", name, exp, actual);
                all_ok = false;
            }
        }
        if all_ok {
            eprintln!("[CHECK] all line counts OK");
        }
    });

    contigs.par_iter().enumerate().for_each_with(tx, |tx, (acc_idx, (_acc, seq))| {
        for m in ac.find_overlapping_iter(seq) {
            let idx = m.pattern().as_usize();
            let cut = m.start() as i64 + enzymes[idx].cut_offset;

            tx.send(Hit {
                enzyme_idx:   idx as u32,
                acc_idx:      acc_idx as u32,
                motif_start:  m.start() as u64,
                cut_position: cut,
            })
            .expect("Writer thread died unexpectedly");
        }
    });

    writer_handle.join().expect("Writer thread panicked");
    eprintln!("[OK] done.");
}