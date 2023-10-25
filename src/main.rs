use clap::Parser;
mod exact_count;
use exact_count::back_to_sequences;



///////////////////////// MAIN /////////////////////////

/// Extract sequences that contain some kmers
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input fasta or fastq [.gz] file containing the original sequences (eg. reads)
    #[arg(long)]
    in_sequences: String, 

    /// Input fasta file containing the original kmers
    #[arg(long)]
    in_kmers: String, 

    /// Output file containing the filtered original sequences (eg. reads). 
    /// It will be automatically in fasta or fastq format depending on the input file.
    #[arg(long)]
    out_sequences: String, 

    /// If provided, output text file containing the kmers that occur in the reads with their number of occurrences
    #[arg(long, default_value_t = String::from(""))]
    out_kmers: String,

    /// Size of the kmers to index and search
    #[arg(short, long, default_value_t = 31)]
    kmer_size: usize,

    /// Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
    /// Minimal threshold of the ratio  (%) of kmers that must be found in a sequence to keep it (default 0%). 
    /// Thus by default, if no kmer is found in a sequence, it is not output.
    #[arg(short, long, default_value_t = 0.0)]
    min_threshold: f32,

    
    /// Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
    /// Maximal threshold of the ratio (%) of kmers that must be found in a sequence to keep it (default 100%). 
    /// Thus by default, there is no limitation on the maximal number of kmers found in a sequence. 
    #[arg(long, default_value_t = 100.0)]
    max_threshold: f32,

    /// Used original kmer strand (else canonical kmers are considered)
    #[arg(long, default_value_t = false)]
    stranded: bool,

    /// Query the reverse complement of reads. Useless without the --stranded option
    #[arg(long, default_value_t = false)]
    query_reverse: bool,
}

fn main() {
    let args = Args::parse();
    if args.stranded == false && args.query_reverse == true {
        eprintln!("Warning: --query-reverse is useless without --stranded");
    }
    if args.min_threshold > args.max_threshold {
        panic!("Error: --min-threshold must be <= --max-threshold");
    }
    let _ = back_to_sequences(args.in_sequences, 
        args.in_kmers, 
        args.out_sequences, 
        args.out_kmers, 
        args.kmer_size,
        args.min_threshold, 
        args.max_threshold, 
        args.stranded,
        args.query_reverse);
}
