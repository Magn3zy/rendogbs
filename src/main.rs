// Input:  --enzymes-run  enzymes_run.csv
//         --reference    genome.fa / .fa.gz
//         --out-dir      output directory
//         --threads      N contigs processed in parallel
//
// Output: cuts/EcoRI.csv, cuts/MseI.csv ...
//        columns: accession,strand,motif_start,cut_position

use std::{
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    path::PathBuf,
    // sync::{Arc, Mutex}, zamykat/nezamykat - bud kazde vlakno ma hashovaci tabulku a slouci se to nakonec nezapisuje si to navzajem, nebo teda zapisuje navzajem a mutex hlida aby se nezapisovalo pres sebe vlastne, asi otestovat co bude rychlejsi ten lock by mohl delat problemy
};
// hashovaci tabulka, filesystem, data do bufferu zapis po blocich, objekt reprez cestu, arc a mutex pro paralelni beh

use aho_corasick::{AhoCorasick};
use clap::Parser;
use csv::ReaderBuilder;
use rayon::prelude::*;
// aho, csv čtení, cli parser, parallel

#[derive(Parser)]
#[command(name = "rendogbs_finder")]
// clap generuje parser argumentu, struktura drzici parametry z prikazove radky
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

// start programu, tvorba kofigurace poolu, tvori globalni pool, vytvoreni adresre vystupu pokud neni, pole pro enzymy, otevreni csv s motivy
fn main() {
    let cli = Cli::parse(); // parsovani argumentu

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global() // nastaveni limitu pro cely program
        .unwrap();

    fs::create_dir_all(&cli.out_dir).unwrap(); // vytvoreni adresare results

    // nacteni dat enzymu z csv
    let mut names: Vec<String> = Vec::new();
    let mut seqs: Vec<String> = Vec::new();
    let mut cuts: Vec<i64> = Vec::new();

    let mut rdr = ReaderBuilder::new()
        .has_headers(true) // preskoceni hlavicky
        .from_path(&cli.enzymes_run)
        .expect("Cannot open enzymes_run.csv");

    for rec in rdr.records() {
        let rec = rec.unwrap();
        names.push(rec[0].to_string()); // nazev to string
        seqs.push(rec[1].to_string()); // string
        cuts.push(rec[2].parse().unwrap()); // cislo i64
    }

    eprintln!("[INFO] {} patterns loaded", seqs.len());

    let ac = AhoCorasick::builder()
            .ascii_case_insensitive(true) // test bez prevodu fasta
            .build(&seqs)
            .expect("AhoCorasick build failed");

    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new(); // sekvence jako bajty
    let mut reader = needletail::parse_fastx_file(&cli.reference)
        .expect("Cannot read FASTA");
// smycka pro precteni fasta
    while let Some(rec) = reader.next() {
        let rec = rec.unwrap(); // contig po contigu do konce souboru
        let acc = std::str::from_utf8(rec.id()) // prevedeni kodovani hlavicky z utf8 na str
            .unwrap()
            .split_whitespace() // ignoruje > zacatek nazvu contigu
            .next()
            .unwrap_or("?")
            .to_string();
        let seq = rec.seq().to_vec();
        contigs.push((acc, seq));
    }

    eprintln!("[INFO] {} contigs | {} threads", contigs.len(), cli.threads);

    // parallel contigs
    let merged: HashMap<String, Vec<String>> = contigs
        .par_iter()
        .map(|(acc, seq)| {
            let mut local: HashMap<String, Vec<String>> = HashMap::new();
            for m in ac.find_overlapping_iter(seq) { // najde vsechny shody!
                let idx = m.pattern().as_usize(); // relativni cut pozice
                let cut = m.start() as i64 + cuts[idx]; // cut position absolutni ve vlakne
                local
                    .entry(names[idx].clone()) //klic nazev enzymu
                    .or_default() // neexistujici klic tzn prazdny vektor
                    .push(format!("{},{},{}", acc, m.start(), cut));
            }
            local
        })
        // slouceni HashMaps z vlaken a presunuti dat
        .reduce(HashMap::new, |mut a, b| {
            for (k, mut v) in b {
                a.entry(k).or_default().append(&mut v);
            }
            a
        });

    // csv output
    let mut enzyme_set: Vec<String> = names.clone();
    enzyme_set.sort(); // abecedni serazeni aby fungovalo odstraneni duplicit
    enzyme_set.dedup(); // na serazenem vektoru

    for name in &enzyme_set {
        let path = cli.out_dir.join(format!("{}.csv", name.replace(['/', ':'], "_"))); // nahrazeni nevhodnych znaku sobouru 
        let mut w = BufWriter::new(fs::File::create(&path).unwrap()); // otevreni souboru pro zapis
        writeln!(w, "accession,motif_start,cut_position").unwrap(); // hlavicka
        for row in merged.get(name).into_iter().flatten() { 
            writeln!(w, "{}", row).unwrap();
        }
        eprintln!("  -> {}.csv", name);
    }

    eprintln!("[OK] done.");
}