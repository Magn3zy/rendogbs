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
    path::PathBuf,
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

// message through channel = one hit
struct Hit {
    enzyme_name: String,
    row: String, // "accession,motif_start,cut_position"
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

fn main() {
    let cli = Cli::parse();

    // at leats 1 writer + N-1 searchers
    let search_threads = cli.parallel.saturating_sub(1).max(1);

    rayon::ThreadPoolBuilder::new()
        .num_threads(search_threads)
        .build_global()
        .unwrap();

    fs::create_dir_all(&cli.out_dir).unwrap();

    // load enzymes 
    let mut enzymes: Vec<Enzyme> = Vec::new();
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(&cli.enzymes_run)
        .expect("Cannot open enzymes_run.csv");

    for rec in rdr.records() {
        let rec = rec.expect("Invalid CSV row");
        enzymes.push(Enzyme {
            name: rec[0].to_string(),
            sequence: rec[1].to_string(),
            cut_offset: rec[2].parse().expect("Invalid cut_offset"),
        });
    }
    eprintln!("[INFO] {} patterns loaded", enzymes.len());

    // build AhoCorasick
    let patterns: Vec<&str> = enzymes.iter().map(|e| e.sequence.as_str()).collect();
    let ac = AhoCorasick::builder()
        .ascii_case_insensitive(true)
        .build(&patterns)
        .expect("AhoCorasick build failed");

    // load fasta
    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = needletail::parse_fastx_file(&cli.r#ref)
        .expect("Cannot read FASTA");

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
    eprintln!(
        "[INFO] {} contigs | {} search threads + 1 writer thread",
        contigs.len(),
        search_threads
    );

    // search and writter thead
    let (tx, rx) = mpsc::sync_channel::<Hit>(65536); // backpressure buffer

    let out_dir = cli.out_dir.clone();

    // writter thread outside of rayonpool
    let writer_handle = thread::spawn(move || {
        // lazy open after first Hit
        let mut writers: HashMap<String, BufWriter<fs::File>> = HashMap::new();

        for hit in rx {
            let w = writers.entry(hit.enzyme_name.clone()).or_insert_with(|| {
                let path = out_dir.join(format!("{}.csv", hit.enzyme_name));
                let mut bw = BufWriter::with_capacity(
                    1 << 20, // 1MB buffer
                    fs::File::create(&path).unwrap()
                );
                writeln!(bw, "accession,motif_start,cut_position").unwrap();
                eprintln!("  -> {}.csv (created)", hit.enzyme_name);
                bw
            });

            writeln!(w, "{}", hit.row).unwrap();
        }

        // channel closed - flush writter
        for (_, mut w) in writers {
            w.flush().unwrap();
        }
    });

    //  parallel search par_iter between contigs and find_overlapping_iter
    contigs.par_iter().for_each_with(tx, |tx, (acc, seq)| {
        for m in ac.find_overlapping_iter(seq) {
            let idx = m.pattern().as_usize();
            let enzyme = &enzymes[idx];
            let cut = m.start() as i64 + enzyme.cut_offset;

            tx.send(Hit {
                enzyme_name: enzyme.name.clone(),
                row: format!("{},{},{}", acc, m.start(), cut),
            })
            .expect("Writer thread died unexpectedly");
        }
    });
    // tx drop rx loop inside writter ended

    writer_handle.join().expect("Writer thread panicked");
    eprintln!("[OK] done.");
}