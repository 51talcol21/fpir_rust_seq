use std::fs::File;
use std::io::{self, BufRead};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FASTQRecord {
    pub id: String,
    pub sequence: String,
    pub quality: String
}

pub fn parse_fastq_file(path: &str) -> std::io::Result<Vec<FASTQRecord>> {
    let file = File::open(path).expect("File is unable to be opened.");
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut sequence: String = String::new();

    let mut sequences_total: Vec<FASTQRecord> = Vec::new();

    while let Some(Ok(line)) = lines.next() {
        if line.starts_with("@") {
            let id = line[1..].trim().to_string();
            while let Some(Ok(seq_line)) = lines.next() {
                if seq_line == "+" {
                    break;
                }
                sequence.push_str(&seq_line);
            }
            let mut quality = line[..].trim().to_string();
            let new_struct = FASTQRecord {
                id: id.clone(),
                sequence: sequence.clone(),
                quality: quality.clone()
            };
            sequences_total.push(new_struct);
            sequence.clear();
            quality.clear();
        }
    }

    Ok(sequences_total)
}