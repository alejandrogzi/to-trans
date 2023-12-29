//! to-trans: a high-performance transcriptome builder from fasta + GTF/GFF
//!
//! This is a command line tool that builds a transcriptome
//! from a fasta file and a GTF/GFF file.
//!
//! Usage:
//! to-trans <fasta> <gtf/gff> [<opt>] [<output>] [<threads>]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use rayon::prelude::*;

use colored::Colorize;

use num_cpus;

use clap::{self, Parser};

use to_trans::*;

const COMPLEMENT: [u8; 128] = {
    let mut nt = [0; 128];
    nt[b'A' as usize] = b'T';
    nt[b'T' as usize] = b'A';
    nt[b'C' as usize] = b'G';
    nt[b'G' as usize] = b'C';
    nt[b'a' as usize] = b't';
    nt[b't' as usize] = b'a';
    nt[b'c' as usize] = b'g';
    nt[b'g' as usize] = b'c';
    nt[b'N' as usize] = b'N';
    nt[b'n' as usize] = b'n';
    nt
};

#[derive(Parser, Debug)]
#[command(
    author = "Alejandro Gonzales-Irribarren",
    version = "0.2.0",
    about = "High-performance transcriptome builder from fasta + GTF/GFF"
)]
struct Args {
    /// FASTA file
    #[clap(short = 'f', long, help = "Path to .fa file", value_name = "FASTA")]
    fasta: PathBuf,

    /// GTF file
    #[clap(
        short = 'g',
        long,
        help = "Path to annotation file",
        value_name = "GTF/GFF"
    )]
    gtf: PathBuf,

    // Exon or CDS
    #[clap(
        short = 'm',
        long,
        help = "Feature to extract from GTF/GFF file (exon or CDS)",
        default_value = "exon",
        value_name = "FEATURE"
    )]
    mode: String,

    /// Output file
    #[clap(
        short,
        long,
        default_value = "transcriptome.fa",
        help = "Path to output file"
    )]
    out: PathBuf,

    /// Number of threads
    #[clap(
        short = 't',
        long,
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    threads: usize,
}

fn main() {
    let time = std::time::Instant::now();

    let args = Args::parse();

    if args.threads == 0 {
        println!(
            "{} {}",
            "Error:".bright_red().bold(),
            "Number of threads must be greater than 0!"
        );
        std::process::exit(1);
    }

    if std::fs::metadata(&args.fasta).unwrap().len() == 0 {
        println!(
            "{} {}",
            "Error:".bright_red().bold(),
            "Input file is empty!"
        );
        std::process::exit(1);
    }

    if args.fasta == args.out {
        println!(
            "{} {}",
            "Error:".bright_red().bold(),
            "Input and output files must be different!"
        );
        std::process::exit(1);
    }

    run(args);

    let duration = time.elapsed();
    println!("Elapsed: {:?}", duration);
}

fn run(args: Args) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let gtf = reader(&args.gtf).unwrap_or_else(|e| {
        eprintln!("{} {}", "Error:".bright_red().bold(), e);
        std::process::exit(1);
    });
    let records = parallel_parse(&gtf, &args.mode).unwrap_or_else(|e| {
        eprintln!("{} {}", "Error:".bright_red().bold(), e);
        std::process::exit(1);
    });
    let seqs = Fasta::read(&args.fasta).unwrap().records;
    let mut writer = BufWriter::new(File::create(args.out).unwrap());

    // {transcript : (chr, strand, [starts], [ends]))}
    let transcripts: HashMap<String, (String, String, Vec<u32>, Vec<u32>)> = records
        .into_par_iter()
        .fold(
            || HashMap::new(), // local accumulator [per thread]
            |mut local_transcripts, record| {
                let (chr, strand, start, end, transcript) = (
                    record.chr,
                    record.strand,
                    record.start,
                    record.end,
                    record.transcript,
                );

                let entry = local_transcripts
                    .entry(transcript)
                    .or_insert_with(|| (chr, strand, Vec::new(), Vec::new()));
                entry.2.push(start);
                entry.3.push(end);

                local_transcripts
            },
        )
        .reduce(
            || HashMap::new(),
            |mut combined_transcripts, local_transcripts| {
                // merge local accs
                for (transcript, (chr, strand, starts, ends)) in local_transcripts {
                    let combined_entry = combined_transcripts
                        .entry(transcript)
                        .or_insert_with(|| (chr, strand, Vec::new(), Vec::new()));
                    combined_entry.2.extend(starts);
                    combined_entry.3.extend(ends);
                }
                combined_transcripts
            },
        );

    for (transcript, (chr, strand, mut starts, mut ends)) in transcripts {
        if !transcript.is_empty() {
            let mut seq = String::new();
            let mut process_parts = |starts: &mut Vec<u32>, ends: &mut Vec<u32>| {
                for (i, &start) in starts.iter().enumerate() {
                    let end = ends[i];
                    if let Some(part) = get_sequence(&seqs, &chr, start, end, &strand) {
                        seq.push_str(&String::from_utf8(part).expect("Invalid UTF-8 sequence"));
                    }
                }
            };

            match strand.as_str() {
                "+" => {
                    starts.sort_unstable();
                    ends.sort_unstable();

                    process_parts(&mut starts, &mut ends);
                }
                "-" => {
                    starts.sort_unstable_by(|a, b| b.cmp(a));
                    ends.sort_unstable_by(|a, b| b.cmp(a));

                    process_parts(&mut starts, &mut ends);
                }
                _ => (),
            }
            if let Err(e) = writeln!(writer, ">{}\n{}", transcript, seq) {
                eprintln!("Failed to write sequence {}: {}", transcript, e);
            }
        }
    }
}

fn get_sequence(
    seqs: &[seq_io::fasta::OwnedRecord],
    chr: &str,
    start: u32,
    end: u32,
    strand: &str,
) -> Option<Vec<u8>> {
    let mut seq = vec![];
    let start = start - 1;
    for record in seqs {
        //let head = String::from_utf8(record.head.to_vec()).unwrap();
        let x = record.head.split(|&b| b == b' ').next().unwrap().to_vec();
        let head = String::from_utf8(x).unwrap();
        if head.contains(chr) {
            seq = record.seq[start as usize..end as usize].to_vec();

            if strand.contains("-") {
                seq.reverse();
                seq = b_cpm(seq);
            }

            return Some(seq);
        }
    }
    return Some(seq);
}

fn b_cpm(seq: Vec<u8>) -> Vec<u8> {
    let cpm = seq
        .iter()
        .map(|&byte| COMPLEMENT[byte as usize])
        .collect::<Vec<u8>>();
    cpm
}
