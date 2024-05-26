//! Back to sequences: find the origin of kmers

/* std use */

// TODO : limit the number of threads with an option

use std::env;

/* crates use */
use clap::Parser as _;

/* project use */
use back_to_sequences::back_to_sequences;
use back_to_sequences::back_to_multiple_sequences;
use back_to_sequences::cli::Args;
use back_to_sequences::kmer_counter::KmerCounterWithLog;

///////////////////////// MAIN /////////////////////////

fn main() {
    (|| {
        let args = Args::parse();

        // Set the number of threads for rayon
        // If the number of threads is not set, rayon will use the number of logical cores
        env::set_var("RAYON_NUM_THREADS", args.threads.to_string());


        // If out_sequences and out_kmers are not provided, we do nothing, we can quit
        if args.out_sequences.is_empty() && args.out_kmers.is_empty() {
            return Err(eprintln!("Warning: no output file provided, nothing to do"));
        }

        // If out_kmers is not provided but output_kmer_positions is true, warn that it has no effect
        if args.out_kmers == "" && args.output_kmer_positions {
            eprintln!("Warning: --output_kmer_positions has no effect without --out-kmers");
        }

        // If out_kmers is not provided but counted_kmer_threshold is set, this has no effect
        if args.out_kmers == "" && args.counted_kmer_threshold > 0 {
            eprintln!("Warning: --counted-kmer-threshold has no effect without --out-kmers");
        }

        if !args.stranded && args.query_reverse {
            eprintln!("Warning: --query-reverse is useless without --stranded");
        }

        if args.min_threshold > args.max_threshold {
            return Err(eprintln!(
                "Error: --min-threshold must be <= --max-threshold"
            ));
        }

        if args.in_sequences == "" && args.in_filelist != ""{
            if args.out_filelist == ""{
                return Err(eprintln!(
                    "Error: --in-filelist requires --out-filelist"
                ));
            }

            if args.output_kmer_positions{
                return Err(eprintln!(
                    "Error: --in-filelist and --output-kmer-positions are mutually exclusive (for now)"
                ));
            }
            back_to_multiple_sequences(
                args.in_filelist,
                args.in_kmers,
                args.out_filelist,
                args.out_kmers,
                args.output_mapping_positions,
                args.kmer_size,
                args.counted_kmer_threshold,
                args.min_threshold,
                args.max_threshold,
                args.stranded,
                args.query_reverse,
                args.no_low_complexity,
            )
        }
        else{
            if args.output_kmer_positions{
                // Use KmerCounterWithLog to log the match position of kmers in the reads
                back_to_sequences::<std::sync::Mutex<KmerCounterWithLog>>(
                    args.in_sequences,
                    args.in_kmers,
                    args.out_sequences,
                    args.out_kmers,
                    args.output_mapping_positions,
                    args.kmer_size,
                    args.counted_kmer_threshold,
                    args.min_threshold,
                    args.max_threshold,
                    args.stranded,
                    args.query_reverse,
                    args.no_low_complexity,
                )
            }
            else {
                // Use atomic_counter::RelaxedCounter to only count the number of kmers in the reads
                back_to_sequences::<atomic_counter::RelaxedCounter>(
                    args.in_sequences,
                    args.in_kmers,
                    args.out_sequences,
                    args.out_kmers,
                    args.output_mapping_positions,
                    args.kmer_size,
                    args.counted_kmer_threshold,
                    args.min_threshold,
                    args.max_threshold,
                    args.stranded,
                    args.query_reverse,
                    args.no_low_complexity,
                )
            }
        }
>>>>>>> main
    })()
    .map_err(|()| std::process::exit(1))
    .ok();
}
