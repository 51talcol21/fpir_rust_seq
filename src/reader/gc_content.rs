use super::fasta::FASTARecord;

pub enum GENRecord {
    FASTARecord(FASTARecord),
}

pub fn calculate_gc_content(record: &GENRecord) -> f64 {
    let sequence = match record {
        GENRecord::FASTARecord(fasta) => &fasta.sequence,
    };

    let gc_count = sequence.chars().filter(|c| *c == 'g' || *c =='c' || *c == 'G' || *c == 'C').count();
    let gc_percent = match gc_count {
        0 => 0.0,
        _ => (gc_count as f64) / (sequence.len() as f64) * 100.0
    };
    gc_percent
}