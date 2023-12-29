//! to-trans: a high-performance transcriptome builder from fasta + GTF/GFF
//!
//! This is a command line tool that builds a transcriptome
//! from a fasta file and a GTF/GFF file.
//!
//! Usage:
//! to-trans <fasta> <gtf/gff> <opt> <output>

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use rayon::prelude::*;

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
    version = "0.1.0",
    about = "High-performance transcriptome builder from fasta + GTF/GFF"
)]
struct Args {
    /// FASTA file
    #[clap(short, long, help = "Path to .fa file")]
    fasta: PathBuf,
    /// GTF file
    #[clap(short, long, help = "Path to annotation file")]
    gtf: PathBuf,
    // Exon or CDS
    #[clap(
        short,
        long,
        help = "Feature to extract from GTF/GFF file (exon or CDS)",
        default_value = "exon"
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
}

fn main() {
    let time = std::time::Instant::now();

    let args = Args::parse();

    let gtf = reader(&args.gtf).unwrap();
    let records = parallel_parse(&gtf, &args.mode).unwrap();
    let seqs = Fasta::read(&args.fasta).unwrap().records;
    let mut writer = BufWriter::new(File::create(args.out).unwrap());

    // {transcript : (chr, strand, [starts], [ends]))}
    let transcripts: HashMap<String, (String, String, Vec<u32>, Vec<u32>)> = records
        .into_par_iter()
        .fold(
            || HashMap::new(), // Create an empty local HashMap as the accumulator.
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
                    .or_insert_with(|| (chr.clone(), strand.clone(), Vec::new(), Vec::new()));
                entry.2.push(start);
                entry.3.push(end);

                local_transcripts
            },
        )
        .reduce(
            || HashMap::new(), // Create an empty HashMap for the reducing step.
            |mut combined_transcripts, local_transcripts| {
                // Merge the local accumulator into the combined one.
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
            match strand.as_str() {
                "+" => {
                    starts.sort_unstable();
                    ends.sort_unstable();

                    for (i, &start) in starts.iter().enumerate() {
                        let end = ends[i];
                        if let Some(part) = get_sequence(&seqs, &chr, start, end, &strand) {
                            seq.push_str(&String::from_utf8(part).unwrap());
                        }
                    }
                    writeln!(writer, ">{}\n{}", transcript, seq).unwrap();
                }
                "-" => {
                    starts.sort_unstable_by(|a, b| b.cmp(a));
                    ends.sort_unstable_by(|a, b| b.cmp(a));

                    for (i, &start) in starts.iter().enumerate() {
                        let end = ends[i];
                        if let Some(part) = get_sequence(&seqs, &chr, start, end, &strand) {
                            seq.push_str(&String::from_utf8(part).unwrap());
                        }
                    }
                    writeln!(writer, ">{}\n{}", transcript, seq).unwrap();
                }
                _ => (),
            }
        }
    }

    let duration = time.elapsed();
    println!("Elapsed: {:?}", duration);
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
