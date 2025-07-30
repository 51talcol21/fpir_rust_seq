use std::fmt;
use super::fasta::FASTARecord;
use super::fastq::FASTQRecord;

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

#[derive(Debug, Eq, Clone, PartialEq)]
pub struct Percent(pub u32);

impl std::fmt::Display for Percent {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}%", self.0 as f32 / 1000.0)
    }
}