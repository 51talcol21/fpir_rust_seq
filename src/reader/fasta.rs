use std::cmp::{max, min};
use std::fs::File;
use std::io::{self, BufRead};

use super::fastq::SequenceINFO;
use super::gc_content::GENRecord;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FASTARecord {
    pub id: String,
    pub sequence: String,
}

pub fn parse_fasta_file(path: &str) -> std::io::Result<SequenceINFO> {
    let file = File::open(path).expect("File is unable to be opened.");
    let reader = io::BufReader::new(file);

    let mut id: String = String::new();
    let mut sequence: String = String::new();
    let mut sequences_total: Vec<FASTARecord> = Vec::new();
    let mut max_sequence_len = i64::MIN;
    let mut min_sequence_len = i64::MAX;
    let mut mean_sequence_len = 0;

    for line in reader.lines() {
        let each_line = line?;

        if each_line.starts_with(">") {
            // Is start of new sequence, but prior sequence in there.
            // Append prior sequence.
            if !id.is_empty() {
                let new_struct = FASTARecord {
                    id,
                    sequence: sequence.clone()
                };
                sequences_total.push(new_struct);
                max_sequence_len = max(max_sequence_len, sequence.len() as i64);
                min_sequence_len = min(min_sequence_len, sequence.len() as i64);
                mean_sequence_len = (mean_sequence_len + sequence.len()) / sequences_total.len();

                id = String::new();
            }

            id = each_line[1..].trim().to_string();
            sequence.clear();

        } else {
            sequence.push_str(&each_line);
        }
    }
    let new_struct = FASTARecord {
        id,
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