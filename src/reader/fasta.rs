use std::cmp::{max, min};
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead};

use super::fastq::SequenceINFO;
use super::gc_content::GENRecord;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FASTARecord {
    pub id: String,
    pub sequence: String,
}

impl std::fmt::Display for FASTARecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Sequence ID:'{}' \n Sequence:'{}", self.id, self.sequence)
    }
}

pub fn parse_fasta_file(path: &str) -> std::io::Result<SequenceINFO> {
    let file = File::open(path).expect("File is unable to be opened.");
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut sequence: String = String::new();
    let mut sequences_total: Vec<FASTARecord> = Vec::new();
    let mut max_sequence_len = i64::MIN;
    let mut min_sequence_len = i64::MAX;
    let mut mean_sequence_len = 0;
    let mut id: String = String::new();

    while let Some(Ok(line)) = lines.next() {
        if line.starts_with(">") {
            if !id.is_empty() {
                let new_struct = FASTARecord {
                    id: id.clone(),
                    sequence: sequence.clone()
                };
                sequences_total.push(new_struct);
                max_sequence_len = max(max_sequence_len, sequence.len() as i64);
                min_sequence_len = min(min_sequence_len, sequence.len() as i64);
                mean_sequence_len = (mean_sequence_len + sequence.len()) / sequences_total.len();
    
                sequence.clear();
            }
            id = line[1..].trim().to_string();
        }
        else {
            sequence.push_str(&line);
        }
    }

    let new_struct = FASTARecord {
        id: id.clone(),
        sequence: sequence.clone()
    };
    sequences_total.push(new_struct);
    max_sequence_len = max(max_sequence_len, sequence.len() as i64);
    min_sequence_len = min(min_sequence_len, sequence.len() as i64);
    mean_sequence_len = (mean_sequence_len + sequence.len()) / sequences_total.len();

    let sequences_genrecord: Vec<GENRecord> = sequences_total
        .into_iter()
        .map(|rec| GENRecord::FASTARecord(rec))
        .collect();

    let sequence_struct = SequenceINFO {
        sequences: sequences_genrecord,
        sequences_min: min_sequence_len,
        sequences_max: max_sequence_len,
        sequences_mean: mean_sequence_len,
    };

    Ok(sequence_struct)
}