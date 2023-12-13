//! Back to sequences: find the origin of kmers

#![warn(missing_docs)]

/* std use */
use std::io::Write as _;

/* crates use */
use atomic_counter::{AtomicCounter, RelaxedCounter};

/* mod declarations */
pub mod cli;
pub mod consts;
pub mod count;
pub mod kmer_hash;
pub mod sequence_normalizer;


/* project use */


/// Extract sequences that contain some kmers and
/// output the kmers that occur in the reads with their number of occurrences
pub fn back_to_sequences(
    in_fasta_reads: String,
    in_fasta_kmers: String,
    out_fasta_reads: String,
    out_txt_kmers: String,
    kmer_size: usize,
    min_threshold: f32,
    max_threshold: f32,
    stranded: bool,
    query_reverse: bool,
    no_low_complexity: bool,
) -> Result<(), ()> {
    // check that in_fasta_reads is a non empty file if it exists:
    if !in_fasta_reads.is_empty() {
        cli::validate_non_empty_file(in_fasta_reads.clone())?;
    }
    cli::validate_non_empty_file(in_fasta_kmers.clone())?;
    // check that in_fasta_kmers is a non empty file:

    let (kmer_set, kmer_size) =
        kmer_hash::index_kmers::<RelaxedCounter>(in_fasta_kmers, kmer_size, stranded, no_low_complexity)
            .map_err(|e| eprintln!("Error indexing kmers: {}", e))?;

    count::kmers_in_fasta_file_par(
        in_fasta_reads,
        &kmer_set,
        kmer_size,
        out_fasta_reads.clone(),
        min_threshold,
        max_threshold,
        stranded,
        query_reverse,
    )?;
    println!(
        "Filtered sequences with exact kmer count are in file {}",
        out_fasta_reads
    );

    // if the out_kmers_file is not empty, we output counted kmers in the out_kmers_file file
    if !out_txt_kmers.is_empty() {
        (|| -> std::io::Result<_> {
            // prints all kmers from kmer_set that have a count > 0
            let mut output = std::fs::File::create(&out_txt_kmers)?;
            for (kmer, count) in kmer_set.iter() {
                if count.get() > 0 {
                    output.write_all(kmer)?;
                    writeln!(output, " {}", count.get())?;
                }
            }
            Ok(())
        })()
        .map_err(|e| eprintln!("Error writing the kmers file: {}", e))?;

        println!(
            "kmers with their number of occurrences in the original sequences are in file {}",
            out_txt_kmers
        );
    }
    Ok(())
}
