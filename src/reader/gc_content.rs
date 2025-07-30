use super::{fasta::FASTARecord, fastq::{FASTQRecord, SequenceINFO}};
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GENRecord {
    FASTARecord(FASTARecord),
    FASTQRecord(FASTQRecord)
}

impl GENRecord {
    pub fn by_sequence_length(&self) -> usize {
        match self {
            GENRecord::FASTARecord(seq) => seq.sequence.len(),
            GENRecord::FASTQRecord(seq) => seq.sequence.len(),
        }
    }
}

impl std::fmt::Display for GENRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            GENRecord::FASTARecord(s) => write!(f, "{}", s),
            GENRecord::FASTQRecord(s) => write!(f, "{}", s)
        }
    }
}

pub fn calculate_gc_content(record: &GENRecord) -> f64 {
    let sequence = match record {
        GENRecord::FASTARecord(fasta) => &fasta.sequence,
        GENRecord::FASTQRecord(fastq) => &fastq.sequence,
    };

    let gc_count = sequence.chars().filter(|c| *c == 'g' || *c =='c' || *c == 'G' || *c == 'C').count();
    let gc_percent = match gc_count {
        0 => 0.0,
        _ => (gc_count as f64) / (sequence.len() as f64) * 100.0
    };
    gc_percent
}

pub fn calculate_n50_score(sequences: &mut SequenceINFO) -> usize {
    // Sort by sequence length
    sequences.sequences.sort_by(|a, b| b.by_sequence_length().cmp(&a.by_sequence_length()));
    let mut cumulative = 0;
    for each_seq in sequences.sequences.iter() {
        cumulative += each_seq.by_sequence_length();
        if cumulative >= sequences.total_nucleotides / 2 {
            return each_seq.by_sequence_length();
        }
    }
    0
}