mod parser {
    pub mod fasta;
    pub mod fastq;
    pub mod stats;
    pub mod genrecord;
}

use clap::Parser;
use parser::stats::calculate_n50_score;
use parser::genrecord::{GENRecord, SequenceINFO, ReadLengthStatistics};
use std::io::{Write, BufWriter};
use std::{str::FromStr, fs::File};
use std::fmt;

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
        writeln!(f, "{}", self.read_length_statistics)?;

        for record in &self.sequences {
            writeln!(f, "{}", record)?;
        }

        Ok(())
    }
}

impl fmt::Display for ReadLengthStatistics {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "-------Read Length Statistics-------")?;
        writeln!(f, "Min Length: {}nt", self.sequences_min)?;
        writeln!(f, "Max Length: {}nt", self.sequences_max)?;
        writeln!(f, "Mean Length: {}nt", self.sequences_mean)?;
        writeln!(f, "Median Length: {}nt", self.sequences_median)?;

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

    /// Compute global GC content
    #[arg(long)]
    gc_global: bool,

    /// Compute per-sequence GC content
    #[arg(long)]
    gc_per_sequence: bool,

    /// Show sequences alongside ID, essentially the fasta/fastq format itself.
    #[arg(long)]
    show_sequences: bool,

    /// Show sequences alongside ID, essentially the fasta/fastq format itself.
    #[arg(long)]
    show_read_statistics: bool,

    /// Output file path
    #[arg(long)]
    output: Option<String>,

    #[arg(long)]
    show_n50: bool,

    #[arg(long, default_value_t = 0)]
    min_length: usize,

    #[arg(long, default_value_t = usize::MAX)]
    max_length: usize,
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

pub struct OutputInfo {
    pub input_file_name: String,
    pub output_file_name: String,
    pub show_gc_global: bool,
    pub show_gc_per_sequence: bool,
    pub show_sequences: bool,
    pub show_read_statistics: bool,
    pub show_n50: bool,
}

pub struct ProcessingValues {
    pub min_length: usize,
    pub max_length: usize,
}

fn write_output(records: &mut SequenceINFO, options: &OutputInfo, input_filters: &ProcessingValues) -> std::io::Result<()> {
    let output_file = File::create(options.output_file_name.clone()).expect("File creation failed");
    let mut output = BufWriter::new(output_file);

    writeln!(output, "File Analyzed: {}\n", options.input_file_name.clone())?;

    if input_filters.min_length > 0 {
        writeln!(output, "Minimum Sequence Length : {}nt", input_filters.min_length)?;
    }

    if input_filters.max_length < usize::MAX {
        writeln!(output, "Maximum Sequence Length : {}nt\n", input_filters.max_length)?;
    }

    if options.show_read_statistics {
        writeln!(output, "{}", records.read_length_statistics)?;
    }

    if options.show_gc_per_sequence {
        writeln!(output, "-------GC Sequence Per Line---------")?;
        writeln!(output, "ID : GC_Percent")?;
        for each_seq in &records.sequences {
            match each_seq {
                GENRecord::FASTARecord(record) => {
                    writeln!(output, "{} : {}", record.id, record.gc_percent)?;
                }
                GENRecord::FASTQRecord(record) => {
                    writeln!(output, "{} : {}", record.id, record.gc_percent)?;
                }
            }
        }
        writeln!(output, "")?;
    }

    if options.show_gc_global {
        writeln!(output, "-------GC Global Statistics---------")?;
        let global_gc = (records.global_gc_count as f64 / records.total_nucleotides as f64) * 100.0;
        writeln!(output, "Global-GC Count: {}\nTotal nt Count: {}", records.global_gc_count, records.total_nucleotides)?;
        writeln!(output, "Global-GC Percentage: {}%\n", global_gc)?;
    }

    if options.show_sequences {
        writeln!(output, "-------All Sequences Shown----------")?;
        for each_seq in &records.sequences {
            match each_seq {
                GENRecord::FASTARecord(record) => {
                    writeln!(output, "{} : {}", record.id, record.sequence)?;
                }
                GENRecord::FASTQRecord(record) => {
                    writeln!(output, "{} : {} : {}", record.id, record.sequence, record.quality)?;
                }
            }
        }
    }

    if options.show_n50 {
        writeln!(output, "\n-----------------N50----------------")?;
        writeln!(output, "N50 value: {}", calculate_n50_score(records))?;
    }

    Ok(())
}

fn main() {
    let args = Args::parse();

    let path = &args.input;
    
    let output_path = match &args.output {
        Some(path) => path,
        None => "output.txt"
    };

    let output_flags = OutputInfo {
        input_file_name: path.to_string(),
        output_file_name: output_path.to_string(),
        show_gc_global: args.gc_global,
        show_gc_per_sequence: args.gc_per_sequence,
        show_read_statistics: args.show_read_statistics,
        show_sequences: args.show_sequences,
        show_n50: args.show_n50
    };

    let input_flags: ProcessingValues = ProcessingValues { 
        min_length: args.min_length,
        max_length: args.max_length,
    };

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
            match result {
                Format::Fasta => {
                    if let Ok(mut result) = parser::fasta::parse_fasta_file(&args.input, &input_flags) {
                        write_output(&mut result, &output_flags, &input_flags).expect("Unable to write to file");
                    }
                }
                Format::Fastq => {
                    if let Ok(mut result) = parser::fastq::parse_fastq_file(&args.input, &input_flags) {
                        write_output(&mut result, &output_flags, &input_flags).expect("Unable to write to file");
                    }
                }
            }
        }
        Err(error) => println!("Error: {}", error),
    }

    println!("Successfully parsed file: '{}'! :)", path);
    ()
}