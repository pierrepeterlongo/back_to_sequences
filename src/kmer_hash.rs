//! Kmer hash declarations

/* std use */

/* crates use */
use ahash::AHashMap as HashMap;
use fxread::initialize_reader;

/* project use */
use crate::sequence_normalizer::SequenceNormalizer;

/// index all kmers of size kmer_size in the fasta file
/// returns a hashmap with the kmers as keys and their count as values, initialized to 0
pub fn index_kmers<T: Default>(
    file_name: String,
    kmer_size: usize,
    stranded: bool,
) -> anyhow::Result<(HashMap<Vec<u8>, T>, usize)> {
    let mut kmer_set = HashMap::new();
    let reverse_complement = if stranded { Some(false) } else { None };

    let mut reader = initialize_reader(&file_name)?;
    loop {
        let Some(record) = reader.next_record()? else {
            break;
        };
        let acgt_sequence = record.seq();
        // for each kmer of the sequence, insert it in the kmer_set
        for i in 0..(acgt_sequence.len() - kmer_size + 1) {
            let kmer = &acgt_sequence[i..(i + kmer_size)];
            kmer_set.insert(
                SequenceNormalizer::new(kmer, reverse_complement)
                    .iter()
                    .collect(),
                Default::default(), // RelaxedCounter::new(0)
            );
        }
    }
    println!(
        "Indexed {} kmers, each of size {}",
        kmer_set.len(),
        kmer_size
    );

    Ok((kmer_set, kmer_size))
}
