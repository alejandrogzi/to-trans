//! to-trans: a high-performance transcriptome builder from fasta + GTF/GFF
//!
//! This is a command line tool that builds a transcriptome
//! from a fasta file and a GTF/GFF file.
//!
//! Usage:
//! to-trans <fasta> <gtf/gff> <opt> <output>

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

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
    #[clap(short, long, help = "Path to your .fa file")]
    fasta: PathBuf,
    /// GTF file
    #[clap(short, long, help = "Path to your .gtf file")]
    gtf: PathBuf,
    // Exon or CDS
    #[clap(
        short,
        long,
        help = "Feature to extract from GTF/GFF file [exon or CDS]. Default to exon",
        default_value = "exon"
    )]
    opt: String,
    /// Output file
    #[clap(
        short,
        long,
        default_value = "transcriptome.fa",
        help = "Path to output file. Default to ./transcriptome.fa"
    )]
    out: PathBuf,
}

fn main() {
    let time = std::time::Instant::now();

    let args = Args::parse();

    let transcripts = GeneModel::parse(&args.gtf, args.opt);
    let mut records = Fasta::read(&args.fasta).unwrap().records;
    let mut output = File::create(&args.out).unwrap();

    for (transcript, attributes) in transcripts {
        let chr = &attributes["chr"][0];
        let strand = &attributes["strand"][0];
        let mut starts = attributes["start"]
            .iter()
            .map(|x| x.parse::<u32>().unwrap())
            .collect::<Vec<u32>>();
        let mut ends = attributes["end"]
            .iter()
            .map(|x| x.parse::<u32>().unwrap())
            .collect::<Vec<u32>>();
        let mut seq = String::new();

        match strand.as_str() {
            "+" => {
                starts.sort_unstable();
                ends.sort_unstable();

                for (i, start) in starts.iter().enumerate() {
                    let end = ends[i];
                    let part = get_sequence(&mut records, chr, *start, end, strand).unwrap();
                    seq.push_str(&String::from_utf8(part).unwrap());
                }
            }
            "-" => {
                starts.sort_unstable_by(|a, b| b.cmp(a));
                ends.sort_unstable_by(|a, b| b.cmp(a));

                for (i, start) in starts.iter().enumerate() {
                    let end = ends[i];
                    let part = get_sequence(&mut records, chr, *start, end, strand).unwrap();
                    seq.push_str(&String::from_utf8(part).unwrap());
                }
            }
            _ => continue,
        }
        output
            .write_all(format!(">{}\n{}\n", transcript, seq).as_bytes())
            .unwrap();
    }

    let duration = time.elapsed();
    println!("Elapsed: {:?}", duration);
}

fn get_sequence(
    records: &Vec<seq_io::fasta::OwnedRecord>,
    chr: &str,
    start: u32,
    end: u32,
    strand: &String,
) -> Option<Vec<u8>> {
    let mut seq = vec![];
    let start = start - 1;
    for record in records {
        let head = String::from_utf8(record.head.to_vec()).unwrap();
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
