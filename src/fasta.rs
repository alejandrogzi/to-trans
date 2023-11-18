use seq_io::fasta::{OwnedRecord, Reader};
use std::path::PathBuf;

pub struct Fasta {
    pub records: Vec<OwnedRecord>,
}

impl Fasta {
    pub fn read(file: &PathBuf) -> Result<Fasta, seq_io::fasta::Error> {
        let mut reader = Reader::from_path(file)?;
        let records: Result<Vec<_>, _> = reader.records().collect();
        let records = records?;
        Ok(Fasta { records })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_read() {
        let fasta_data = b">seq1\nATCG\n>seq2\nGCTA\n";

        let tmp = tempdir::TempDir::new("test_fasta_read").unwrap();
        let path = tmp.path().join("test.fasta");
        std::fs::write(&path, fasta_data).unwrap();

        let result = Fasta::read(&path).unwrap().records;

        assert_eq!(result.len(), 2);
        assert_eq!(result[0].head, b"seq1");
        assert_eq!(result[0].seq, b"ATCG");
        assert_eq!(result[1].head, b"seq2");
        assert_eq!(result[1].seq, b"GCTA");
    }
}
