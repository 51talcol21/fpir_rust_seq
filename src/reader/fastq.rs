use std::cmp::{max, min};
use std::fs::File;
use std::fmt;
use std::io::{self, BufRead};

use super::fasta::Percent;
use super::gc_content::{GENRecord, calculate_gc_content};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FASTQRecord {
    pub id: String,
    pub sequence: String,
    pub quality: String,
    pub gc_percent: Percent,
    pub gc_count: usize,
}

impl std::fmt::Display for FASTQRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Sequence ID:'{}' \n Sequence:'{}' \n Quality:'{}'\n", self.id, self.sequence, self.quality)
    }
}

pub struct SequenceINFO {
    pub sequences: Vec<GENRecord>,
    pub read_length_statistics: ReadLengthStatistics,
    pub global_gc_count: usize,
    pub total_nucleotides: usize,
}

pub struct ReadLengthStatistics {
    pub sequences_min: i64,
    pub sequences_max: i64,
    pub sequences_mean: usize,
    pub sequences_median: usize,
}

pub fn parse_fastq_file(path: &str) -> std::io::Result<SequenceINFO> {
    let file = File::open(path).expect("File is unable to be opened.");
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut sequence: String = String::new();

    let mut sequences_total: Vec<GENRecord> = Vec::new();
    let mut max_sequence_len = i64::MIN;
    let mut min_sequence_len = i64::MAX;
    let mut mean_sequence_len = 0;

    let mut global_gc_count = 0;
    let mut total_nucleotides = 0;

    while let Some(Ok(line)) = lines.next() {
        if line.starts_with("@") {
            let mut id = line[1..].trim().to_string();
            while let Some(Ok(seq_line)) = lines.next() {
                if seq_line == "+" {
                    if let Some(Ok(seq_line)) = lines.next() {
                        let mut quality = seq_line[..].trim().to_string();
                        let gc_count = sequence.chars().filter(|c| *c == 'g' || *c =='c' || *c == 'G' || *c == 'C').count();

                        let mut new_struct = FASTQRecord {
                            id: id.clone(),
                            sequence: sequence.clone(),
                            quality: quality.clone(),
                            gc_percent: Percent(0),
                            gc_count,
                        };

                        total_nucleotides += sequence.len();
                        global_gc_count += gc_count;
                        new_struct.gc_percent = Percent((calculate_gc_content(&GENRecord::FASTQRecord(new_struct.clone())) * 1000.00) as u32);
                        sequences_total.push(GENRecord::FASTQRecord(new_struct));
                        max_sequence_len = max(max_sequence_len, sequence.len() as i64);
                        min_sequence_len = min(min_sequence_len, sequence.len() as i64);
                        mean_sequence_len = (mean_sequence_len + sequence.len()) / sequences_total.len();
            
                        sequence.clear();
                        quality.clear();
                        id.clear();
                        break
                    }
                }
                sequence.push_str(&seq_line);
            }
        }
    }

    let mut sequences_genrecord: Vec<GENRecord> = sequences_total
        .into_iter()
        .map(|rec| rec)
        .collect();

    let sequences_len = sequences_genrecord.len();
    sequences_genrecord.select_nth_unstable_by(sequences_len / 2, |sec_a, seq_b| sec_a.by_sequence_length().cmp(&seq_b.by_sequence_length()));
    let median = &sequences_genrecord[sequences_genrecord.len() / 2].by_sequence_length();
    
    let sequence_statistic_struct = ReadLengthStatistics {
        sequences_min: min_sequence_len,
        sequences_max: max_sequence_len,
        sequences_mean: mean_sequence_len,
        sequences_median: *median,
    };

    let sequence_struct = SequenceINFO {
        sequences: sequences_genrecord,
        read_length_statistics: sequence_statistic_struct,
        global_gc_count,
        total_nucleotides,
    };

    Ok(sequence_struct)
}