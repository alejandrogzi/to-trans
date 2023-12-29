#![allow(dead_code)]

use crate::errors::*;
use memchr::memchr;

#[derive(Debug)]
pub struct GeneModel {
    pub chr: String,
    pub feature: String,
    pub start: u32,
    pub end: u32,
    pub strand: String,
    pub transcript: String,
}

impl GeneModel {
    pub fn parse(line: &str, opt: &str) -> Self {
        let fields = line.split('\t').collect::<Vec<&str>>();
        let feature = fields[2];

        if feature == opt {
            let id = de_attr(fields[8]).unwrap();
            // let id = mem_attr(fields[8]);
            GeneModel {
                chr: fields[0].to_string(),
                feature: feature.to_string(),
                start: fields[3].parse::<u32>().unwrap(),
                end: fields[4].parse::<u32>().unwrap(),
                strand: fields[6].to_string(),
                transcript: id,
            }
        } else {
            GeneModel {
                chr: String::new(),
                feature: String::new(),
                start: 0,
                end: 0,
                strand: String::new(),
                transcript: String::new(),
            }
        }
    }
}

pub fn de_attr(input: &str) -> Result<String> {
    let bytes = input.as_bytes().iter().enumerate();
    let mut incr = 0;
    let mut id = String::new();

    if input.contains("=") {
        incr += 1
    } else {
        incr += 2
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

    if id.is_empty() {
        return Err(Error::AttributeError(input.to_string()));
    }

    Ok(id)
}

pub fn pair(line: &str) -> Result<(String, String)> {
    let mut bytes = line.as_bytes().iter();
    let i = if let Some(pos) = bytes.position(|b| *b == b' ' || *b == b'=') {
        pos
    } else {
        line.len()
    };

    if i + 1 > line.len() {
        return Err(Error::ParseAttributeError(line.to_string()));
    }

    let key = &line[..i];
    let value = &line[i + 1..].trim_matches('"');
    // let value = get_transcript(*&line[i + 1..].trim_matches('"'), ".")?;

    if key.is_empty() || value.is_empty() {
        return Err(Error::ParseAttributeError(line.to_string()));
    }

    Ok((key.to_string(), value.to_string()))
}

fn mem_attr(line: &str) -> String {
    let bytes = line.as_bytes();
    let core = &bytes[memchr(b't', bytes).unwrap()..];
    let comma = &core[..memchr(b';', core).unwrap()];
    let quote = &comma[memchr(b'"', comma).unwrap() + 1..];

    let transcript = String::from_utf8(quote.to_vec())
        .unwrap()
        .trim_end_matches('"')
        .to_string();

    transcript
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gff_de_attr() {
        let input = "gene_id=ENSG00000223972;transcript_id=ENST00000406473;";
        assert_eq!(de_attr(input).unwrap(), "ENST00000406473");
    }

    #[test]
    fn gtf_de_attr() {
        let input = "gene_id \"ENSG00000223972\"; transcript_id \"ENST00000406473\";";
        assert_eq!(de_attr(input).unwrap(), "ENST00000406473");
    }

    #[test]
    fn gff_pair() {
        let line = "transcript_id=ENST00000406473";
        assert_eq!(
            pair(line),
            Ok(("transcript_id".to_string(), "ENST00000406473".to_string()))
        );
    }

    #[test]
    fn gtf_pair() {
        let line = "transcript_id \"ENST00000406473\"";
        assert_eq!(
            pair(line),
            Ok(("transcript_id".to_string(), "ENST00000406473".to_string()))
        );
    }

    #[test]
    fn attr_error() {
        let line = "transcript_idENST00000406473";
        assert_eq!(
            pair(line),
            Err(Error::ParseAttributeError(line.to_string()))
        );
    }

    #[test]
    fn empty_attr() {
        let line = "";
        assert_eq!(de_attr(line), Err(Error::AttributeError(line.to_string())));
    }
}
