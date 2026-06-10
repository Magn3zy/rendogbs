/// Input:  --enzymes-run  enzymes_run.csv
///         --reference    genome.fa / .fa.gz
///         --out-dir      output directory
///         --threads      N contigs processed in parallel
///
/// Output: cuts/EcoRI.csv, cuts/MseI.csv ...
///         columns: accession,strand,motif_start,cut_position

use std::{
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::{Arc, Mutex},
};
/// hashovaci tabulka, filesystem, data do bufferu zapis po blocich, objekt reprez cestu, arc a mutex pro paralelni beh

use aho_corasick::{AhoCorasick, MatchKind};
use clap::Parser;
use csv::ReaderBuilder;
use rayon::prelude::*;
//// aho, csv čtení, cli parser,

#[derive(Parser)]
#[command(name = "rendogbs_finder")]
//// clap generuje parser argumentu, struktura drzici parametry z prikazove radky
struct Cli {
    #[arg(long)]
    enzymes_run: PathBuf,

    #[arg(long, short = 'r')]
    reference: PathBuf,

    #[arg(long, short = 'o')]
    out_dir: PathBuf,

    #[arg(long, short = 't', default_value_t = 1)]
    threads: usize,
}

/// start programu, tvorba kofigurace poolu, tvori globalni pool, vytvoreni adresre vystupu pokud neni, pole pro enzymy, otevreni csv s motivy
fn main() {
    let cli = Cli::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    fs::create_dir_all(&cli.out_dir).unwrap();

    // columns: enzyme_name, expanded_sequence, cut_offset, ...
    let mut names: Vec<String> = Vec::new();
    let mut seqs: Vec<String> = Vec::new();
    let mut cuts: Vec<i64> = Vec::new();

    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(&cli.enzymes_run)
        .expect("Cannot open enzymes_run.csv");

    for rec in rdr.records() {
        let rec = rec.unwrap();

        let seq = rec[1].to_uppercase();

        names.push(rec[0].to_string());
        seqs.push(seq);
        cuts.push(rec[2].parse().unwrap());
    }
    /// pridani po radku motivu, jako string, cut ofset 0 based
    eprintln!("[INFO] {} patterns loaded", seqs.len());
    /// vypsani poctu nactenych motivu
    // Automat for all motifs, sdilene vlastnictvi mezi vlakny, builder, normalni prohledavani, postaveni automatu nad všemi motivy, sdilene nazvy a offsety
    let ac = Arc::new(
        AhoCorasick::builder()
            .match_kind(MatchKind::Standard)
            .build(&seqs)
            .expect("AhoCorasick build failed"),
    );

    let names = Arc::new(names);
    let cuts = Arc::new(cuts);

    // fasta load to ram, otevreni po contigach, smycka jede pres contigy dokud jsou, rozdeleni dle mezer, prevedeni pismen
    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new();

    let mut reader = needletail::parse_fastx_file(&cli.reference)
        .expect("Cannot read FASTA");

    while let Some(rec) = reader.next() {
        let rec = rec.unwrap();

        let acc = std::str::from_utf8(rec.id())
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap_or("?")
            .to_string();
        /// 
        let seq: Vec<u8> = rec
            .seq()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();

        contigs.push((acc, seq));
    }

    eprintln!(
        "[INFO] {} contigs | {} threads",
        contigs.len(),
        cli.threads
    );

    // unique enzyme names
    let mut enzyme_set: Vec<String> = names.iter().cloned().collect();
    enzyme_set.sort();
    enzyme_set.dedup();

    // shared results
    let results: Arc<Mutex<HashMap<String, Vec<String>>>> =
        Arc::new(Mutex::new(
            enzyme_set
                .iter()
                .map(|n| (n.clone(), Vec::new()))
                .collect(),
        ));

    // parallel contigs
    contigs.par_iter().for_each(|(acc, seq)| {
        let mut local: HashMap<String, Vec<String>> = HashMap::new();

        for m in ac.find_iter(seq) {
            let idx = m.pattern().as_usize();

            let motif_start = m.start() as i64;
            let cut_position = motif_start + cuts[idx];

            local
                .entry(names[idx].clone())
                .or_default()
                .push(format!(
                    "{},{},{},{}",
                    acc,
                    '+',
                    motif_start,
                    cut_position
                ));
        }

        let mut g = results.lock().unwrap();

        for (name, mut rows) in local {
            g.entry(name).or_default().append(&mut rows);
        }
    });

    // csv output
    let g = results.lock().unwrap();

    for name in &enzyme_set {
        let safe = name.replace(['/', ':'], "_");

        let path = cli.out_dir.join(format!("{}.csv", safe));

        let mut w = BufWriter::new(
            fs::File::create(&path).unwrap()
        );

        writeln!(
            w,
            "accession,strand,motif_start,cut_position"
        )
        .unwrap();

        if let Some(rows) = g.get(name) {
            for row in rows {
                writeln!(w, "{}", row).unwrap();
            }
        }

        eprintln!("  -> {}.csv", name);
    }

    eprintln!("[OK] done.");
}