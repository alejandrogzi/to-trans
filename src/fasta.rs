use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug)]
pub struct Fasta {
    index: HashMap<Vec<u8>, Vec<u8>>,
}

impl Fasta {
    pub fn new(fa: &str) -> Result<Self, Box<dyn Error>> {
        let reader = BufReader::new(File::open(fa)?);
        let mut index = HashMap::new();
        let mut chr = Vec::new();
        let mut seq = Vec::new();

        for line in reader.lines() {
            let line = line?.as_bytes().to_vec();

            if line.is_empty() {
                return Err("Empty line in FASTA file".into());
            }

            if line[0] == b'>' {
                if !seq.is_empty() && !chr.is_empty() {
                    index.insert(chr.clone(), seq.clone());
                }
                chr = line[1..].to_vec();
                seq.clear();
            } else {
                seq = line.to_vec();
            }
        }

        if !seq.is_empty() && !chr.is_empty() {
            index.insert(chr, seq);
        }

        Ok(Fasta { index })
    }
}