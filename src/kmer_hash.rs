//! Kmer hash declarations

/* std use */

/* crates use */
use ahash::AHashMap as HashMap;
use entropy::shannon_entropy;
use fxread::initialize_reader;

/* project use */
use crate::sequence_normalizer::SequenceNormalizer;

/// given a kmer as a &[u8] return a tuple boolean, position
/// if the kmer contains an non ACGT letter, return false and the position of the first non ACGT letter
/// else return true and 0 as position
pub fn first_non_acgt(kmer: &[u8]) -> (bool, usize) {
    for (i, &byte) in kmer.iter().enumerate() {
        if byte != b'A' && byte != b'C' && byte != b'G' && byte != b'T' {
            return (false, i);
        }
    }
    (true, 0)
}


// /// given a kmer as a &[u8] check that it contains only ACGT letters
// /// return true if it is the case, false otherwise
// fn is_acgt(kmer: &[u8]) -> bool {
//     for &byte in kmer {
//         if byte != b'A' && byte != b'C' && byte != b'G' && byte != b'T' {
//             return false;
//         }
//     }
//     true
// }



/// index all kmers of size kmer_size in the fasta file
/// returns a hashmap with the kmers as keys and their count as values, initialized to 0
pub fn index_kmers<T: Default>(
    file_name: String,
    kmer_size: usize,
    stranded: bool,
    no_low_complexity: bool,
) -> anyhow::Result<(HashMap<Vec<u8>, T>, usize)> {
    let mut kmer_set = HashMap::new();
    let reverse_complement = if stranded { Some(false) } else { None };

    let mut reader = initialize_reader(&file_name)?;
    loop {
        let Some(mut record) = reader.next_record()? else {
            break;
        };
        record.upper();
        let acgt_sequence = record.seq();

        
        // for each kmer of the sequence, insert it in the kmer_set
        if acgt_sequence.len() < kmer_size {
            continue;
        }
        let mut i = 0;
        while i < acgt_sequence.len() - kmer_size + 1 {
        // for mut i in 0..(acgt_sequence.len() - kmer_size + 1) {
            let kmer = &acgt_sequence[i..(i + kmer_size)];
            let first_non_acgt = first_non_acgt(kmer);
            if first_non_acgt.0 == false {
                // If the kmer contains a non acgt letter, we jump to the next possible kmer
                i = i + first_non_acgt.1 + 1;
                continue;
            }
            else { 
                // If the entropy is too low, the kmer is not inserted
                if no_low_complexity && shannon_entropy(kmer) < 1.0 { 
                    i = i + 1;
                    continue;
                }
                kmer_set.insert(
                    SequenceNormalizer::new(kmer, reverse_complement)
                        .iter()
                        .collect(),
                    Default::default(), // RelaxedCounter::new(0)
                );
            }
            i = i + 1;
        }
    }
    println!(
        "Indexed {} kmers, each of size {}",
        kmer_set.len(),
        kmer_size
    );

    Ok((kmer_set, kmer_size))
}
