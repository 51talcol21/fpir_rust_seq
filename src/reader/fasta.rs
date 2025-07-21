use std::fs::File;
use std::io::{self, BufRead};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FASTARecord {
    pub id: String,
    pub sequence: String,
}

pub fn parse_fasta_file(path: &str) -> std::io::Result<Vec<FASTARecord>> {
    let file = File::open(path).expect("File is unable to be opened.");
    let reader = io::BufReader::new(file);

    let mut id: String = String::new();
    let mut sequence: String = String::new();
    let mut sequences_total: Vec<FASTARecord> = Vec::new();

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

    Ok(sequences_total)
}