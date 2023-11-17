#![allow(dead_code)]

use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Deserialize)]
pub struct GeneModel {
    seqname: String,
    source: String,
    feature: String,
    start: String,
    end: String,
    score: String,
    strand: String,
    frame: String,
    attribute: String,
}

impl GeneModel {
    pub fn parse(file: &PathBuf) -> HashMap<String, HashMap<String, Vec<String>>> {
        let mut transcripts: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .has_headers(false)
            .from_path(file)
            .unwrap()
            .into_deserialize::<GeneModel>();

        for line in rdr {
            let line = line.unwrap();
            match line.feature.as_str() {
                "exon" => {
                    let id = de_attr(line.attribute);
                    let chr = line.seqname;
                    let s = line.start;
                    let e = line.end;
                    let strand = line.strand;

                    let transcript = transcripts.entry(id).or_insert(HashMap::new());

                    transcript.entry("chr".to_string()).or_insert(vec![chr]);
                    transcript
                        .entry("strand".to_string())
                        .or_insert(vec![strand]);

                    let start = transcript.entry("start".to_string()).or_insert(Vec::new());
                    start.push(s);

                    let end = transcript.entry("end".to_string()).or_insert(Vec::new());
                    end.push(e);
                }
                _ => continue,
            }
        }
        transcripts
    }
}

pub fn de_attr(input: String) -> String {
    let bytes = input.as_bytes().iter().enumerate();
    let mut incr = 0;
    let mut id = String::new();

    if input.contains("=") {
        incr = 1
    } else {
        incr = 2
    }

    let mut start = 0;
    for (i, byte) in bytes {
        if *byte == b';' {
            let word = &input[start..i];
            if !word.is_empty() && word.starts_with("transcript_id") {
                let (_, value) = pair(word).unwrap();
                id = value;
            }
            start = i + incr;
        }
    }
    return id;
}

pub fn pair(line: &str) -> Option<(String, String)> {
    let mut bytes = line.as_bytes().iter();
    let i = if let Some(pos) = bytes.position(|b| *b == b' ' || *b == b'=') {
        pos
    } else {
        line.len()
    };

    let key = &line[..i];
    let value = &line[i + 1..].trim_matches('"');
    // let value = get_transcript(*&line[i + 1..].trim_matches('"'), ".")?;

    Some((key.to_string(), value.to_string()))
}

pub fn get_transcript(transcript: &str, sep: &str) -> Option<String> {
    let bytes = transcript.as_bytes().iter().enumerate();
    let start = 0;
    let mut word = String::new();

    for (i, byte) in bytes {
        if *byte == sep.as_bytes()[0] {
            word.push_str(&transcript[start..i]);
        }
    }
    if !word.is_empty() {
        return Some(word);
    } else {
        return Some(transcript.to_string());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gff_de_attr() {
        let input = "gene_id=ENSG00000223972;transcript_id=ENST00000406473;";
        assert_eq!(de_attr(input.to_string()), "ENST00000406473");
    }

    #[test]
    fn gtf_de_attr() {
        let input = "gene_id \"ENSG00000223972\"; transcript_id \"ENST00000406473\";";
        assert_eq!(de_attr(input.to_string()), "ENST00000406473");
    }

    #[test]
    fn gff_pair() {
        let line = "transcript_id=ENST00000406473";
        assert_eq!(
            pair(line),
            Some(("transcript_id".to_string(), "ENST00000406473".to_string()))
        );
    }

    #[test]
    fn gtf_pair() {
        let line = "transcript_id \"ENST00000406473\"";
        assert_eq!(
            pair(line),
            Some(("transcript_id".to_string(), "ENST00000406473".to_string()))
        );
    }
}
