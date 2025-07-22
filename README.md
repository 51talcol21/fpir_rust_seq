# fpir

fpir, or 'FastX Parser In Rust' is a lightweight, fast, and extensible command-line tool written in Rust for parsing and analyzing FASTA and FASTQ files.

- Parse FASTA and FASTQ files
- Compute per-sequence GC content
- Compute global GC content across all reads
- Calculate basic read length statistics:
  - Mean read length
  - Max read length
- Clean, scriptable command-line interface (CLI)

## Current Features (so far)

- Parse FASTA and FASTQ files
- Compute per-sequence GC content
- Compute global GC content across all reads
- Calculate basic read length statistics:
  - Mean read length
  - Max read length
- Clean, scriptable command-line interface (CLI)

## Planned Features

- More error checking for misformed files/cut off sequences.
- Output in different formats (CSV).
- Median, N50 for fasta.
- Filtering on min/max length, gc content, min quality on fastq.
- Quality score stats 
- Base counts for global and per sequence
- Random subsampling
- More UI for loading/progress.

## Installation

Clone and build from source:

```bash
git clone https://github.com/yourname/fastaseq.git
cd fastaseq
cargo build --release
```

## Usage

After building the project and installing, go to the directory or run it whatever it is.  
Base program is this: `./target/release/fpir`.

To parse a file run the flag `--input FILEINPUT.fasta/fastq`  

By default it outputs to output.txt, but can also specify output file `--output FILENAME.txt`

Current flags
- `--gc-global`: Outputs global GC content
- `--gc-per-sequence`: Outputs the GC content for each ID.
- `--show-read-statistics`: Min/max/mean read length in nt.
- `--show-sequences`: Basically just outputs the original file in the same file as `ID : Sequence : Quality (if fastq)`