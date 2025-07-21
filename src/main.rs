mod reader {
    pub mod fasta;
    pub mod fastq;
    pub mod gc_content;
}

use clap::Parser;
use reader::fastq::SequenceINFO;
use std::str::FromStr;
use std::fmt;

use crate::reader::gc_content::{GENRecord};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

impl FromStr for Format {
    type Err = String;

    fn from_str(input_frmt: &str) -> Result<Format, Self::Err> {
        match input_frmt.to_lowercase().as_str() {
            "fasta" | "fa" => Ok(Format::Fasta),
            "fastq" | "fq" => Ok(Format::Fastq),
            _ => Err(format!("Invalid Format: {}", input_frmt)),
        }
    }
}

impl std::fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = match *self {
            Format::Fasta => "FASTA",
            Format::Fastq => "FASTQ",
        };
        write!(f, "{}", s)
    }
}

impl fmt::Display for SequenceINFO {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "---Sequence Information--")?;
        writeln!(f, "Min Length: {}nt", self.sequences_min)?;
        writeln!(f, "Max Length: {}nt", self.sequences_max)?;
        writeln!(f, "Mean Length: {}nt", self.sequences_mean)?;
        writeln!(f, "--------Sequences--------")?;

        for record in &self.sequences {
            writeln!(f, "{}", record)?;
        }

        Ok(())
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Input file
    #[arg(short, long)]
    input: String,

    /// Choose the format
    #[arg(short, long, value_parser = ["fastq", "fasta"])]
    format: Option<String>,
}

fn detect_format_from_filename(filename: &str) -> Option<&'static str> {
    if filename.ends_with(".fastq") || filename.ends_with(".fastq.gz") {
        Some("FASTQ")
    } else if filename.ends_with(".fasta") || filename.ends_with(".fasta.gz") || filename.ends_with(".fa") {
        Some("FASTA")
    } else {
        None
    }
}

fn main() {
    let args = Args::parse();
    let path = &args.input;

    let format = if let Some(fmt) = args.format {
        Format::from_str(&fmt)
    } else {
        match detect_format_from_filename(path) {
            Some(ext) => Format::from_str(ext),
            None => {
                eprintln!("Could not auto-detect format from file extension.");
                std::process::exit(1);
            }
        }
    };

    match format {
        Ok(result) => {
            println!("The input is {} and the file format is {}.", args.input, result);
            match result {
                Format::Fasta => {
                    if let Ok(result) = reader::fasta::parse_fasta_file(&args.input) {
                        for each_seq in result.sequences.iter() {
                            let our_record = &each_seq;

                            if let GENRecord::FASTARecord(rec) = &our_record {
                                println!("GC percent for record (ID: {}) is {}", rec.id, rec.gc_percent);
                            }
                        }
                        println!("{}", result);
                    }
                }
                Format::Fastq => {
                    if let Ok(result) = reader::fastq::parse_fastq_file(&args.input) {
                        for each_seq in result.sequences.iter() {
                            let our_record = &each_seq;

                            if let GENRecord::FASTQRecord(rec) = &our_record {
                                println!("GC percent for record (ID: {}) is {}", rec.id, rec.gc_percent);
                            }
                        }
                        println!("{}", result);
                    }
                }
            }
        }
        Err(error) => println!("Error: {}", error),
    }

    ()
}