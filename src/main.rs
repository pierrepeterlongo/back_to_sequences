use clap::Parser;
mod exact_count;
use exact_count::validate_kmers;



///////////////////////// MAIN /////////////////////////

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input fasta file containing the original kmers
    #[arg(long)]
    in_fasta_reads: String,

    /// Input fasta file containing the original kmers
    #[arg(long)]
    in_fasta_kmers: String,

    /// Output fasta file containing the original kmers
    #[arg(long)]
    out_fasta_reads: String,

    /// If provided, output text file containing the kmers that occur in the reads with their number of occurrences
    #[arg(long, default_value_t = String::from(""))]
    out_txt_kmers: String,

    /// Number of times to greet
    #[arg(short, long, default_value_t = 31)]
    kmer_size: usize,

    /// Used original kmer strand (else canonical kmers are considered)
    #[arg(long, default_value_t = false)]
    stranded: bool,
}

fn main() {
    let args = Args::parse();
    let _ = validate_kmers(args.in_fasta_reads, 
        args.in_fasta_kmers, 
        args.out_fasta_reads, 
        args.out_txt_kmers, 
        args.kmer_size, 
        args.stranded);
}
