//! Back to sequences: find the origin of kmers

/* std use */

/* crates use */
use clap::Parser as _;

/* project use */
use back_to_sequences::back_to_sequences;
use back_to_sequences::cli::Args;

///////////////////////// MAIN /////////////////////////

fn main() {
    (|| {
        let args = Args::parse();

        if !args.stranded && args.query_reverse {
            eprintln!("Warning: --query-reverse is useless without --stranded");
        }

        if args.min_threshold > args.max_threshold {
            return Err(eprintln!(
                "Error: --min-threshold must be <= --max-threshold"
            ));
        }

        back_to_sequences(
            args.in_sequences,
            args.in_kmers,
            args.out_sequences,
            args.out_kmers,
            args.kmer_size,
            args.min_threshold,
            args.max_threshold,
            args.stranded,
            args.query_reverse,
            args.no_low_complexity,
        )
    })()
    .map_err(|()| std::process::exit(1))
    .ok();
}
