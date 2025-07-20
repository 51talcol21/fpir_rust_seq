mod reader {
    pub mod fasta;
    pub mod gc_content;
}

use clap::Parser;
use std::str::FromStr;
use std::fmt;

use crate::reader::gc_content::{calculate_gc_content, GENRecord};

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

fn main() {
    let args = Args::parse();

    let format = if let Some(fmt) = args.format {
        Format::from_str(&fmt)
    } else {
        println!("We will have to auto detect it!");
        Ok(Format::Fastq)
    };

    match format {
        Ok(result) => {
            println!("The input is {} and the file format is {}!", args.input, result);
            match result {
                Format::Fasta => {
                    if let Ok(_result) = reader::fasta::parse_fasta_file(&args.input) {
                        println!("Successfully ran {:?}", _result);
                        for each_seq in _result.iter() {
                            let ourRec = GENRecord::FASTARecord(each_seq.clone());
                            println!("This is{:?}", calculate_gc_content(&ourRec));
                        }
                    }
                }
                Format::Fastq => {
                    println!("Not implemented yet :(")
                }
            }
        }
        Err(error) => println!("Error: {}", error),
    }

    ()
}