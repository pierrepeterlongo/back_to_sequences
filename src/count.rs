//! Count declarations

/* std use */
use std::fs::File;
use std::io::{BufWriter, Write};

/* crates use */
use ahash::AHashMap as HashMap;
use anyhow::Context as _;

/* project use */
use crate::chunks::{NO_WRITER, Pipeline, WithoutId};
use crate::kmer_counter::KmerCounter;
use crate::matched_sequences::MatchedSequence;
use crate::sequence_normalizer::SequenceNormalizer;


/// Reverse complement a sequence in place.
pub fn rev_comp(seq: &mut [u8]) {
    // Reverse the sequence
    seq.reverse();

    // Complement the sequence
    seq.iter_mut().for_each(|c| {
        if *c & 2 == 0 {
            *c ^= 21;
        } else {
            *c ^= 4;
        }
    });
}

/// for each sequence of a given fasta file, count the number of indexed kmers it contains
/// and output the sequence if its ratio of indexed kmers is in ]min_threshold, max_threshold]
#[allow(clippy::too_many_arguments)]
pub fn kmers_in_fasta_file_par<T, D>(
    file_name: String,
    kmer_set: &HashMap<Vec<u8>, T>,
    kmer_size: usize,
    out_fasta: String,
    min_threshold: f32,
    max_threshold: f32,
    stranded: bool,
    query_reverse: bool,
    map_both_strands: bool,
) -> anyhow::Result<(usize, usize, usize)>
where
    T: KmerCounter,
    D: MatchedSequence + Send + 'static,
{
    let reader = match file_name.as_str() {
        "" => needletail::parse_fastx_stdin()?,
        f  => needletail::parse_fastx_file(f)?,
    };

    let mut output_file = BufWriter::new(
        File::create(out_fasta).context("Error: failed to open the sequence file for writing")?,
    );

    Pipeline::<Option<D>>::run(
        reader,
        // map
        |record| {
            if query_reverse {
                // we need to reverse complement the sequence first
                rev_comp(record.seq);
            }
            let total_nucleotides = record.seq.len();
            let total_kmer = record.seq.len() - kmer_size + 1;

            let proxy_shared_kmers = shared_kmers_par::<_, D>(
                kmer_set,
                &record.seq,
                record.read_id,
                kmer_size,
                stranded,
                map_both_strands,
                );

            let match_kmer = proxy_shared_kmers.match_count();

            *record.extra = Some(proxy_shared_kmers);

            (total_nucleotides, total_kmer, match_kmer)
        },
        // reduce
        (
            || (0, 0, 0),
            |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
        ),
        // writer
        |record| {
            let matched_sequence = record.extra.as_ref().unwrap();
            let percent_shared_kmers = matched_sequence.percent_shared_kmers();

            if percent_shared_kmers > min_threshold && percent_shared_kmers <= max_threshold {
                // supports the user defined thresholds

                // let record_as_string = seq.1.to_owned();
                let iter = record.seq.split(|&x| x == b'\n');


                // output_file.write_all(iter.next().unwrap())?; // write the original header of the record
                let header = std::str::from_utf8(record.id).unwrap();
                writeln!(output_file, ">{}{}", header, matched_sequence)?; // append metrics
                for line in iter {
                    output_file.write_all(line)?;
                    output_file.write_all(b"\n")?;
                }
            } // end read contains at least one indexed kmer
            Ok(())
        })
}

/// for each sequence of a given fasta file, count the number of indexed kmers it contains
pub fn only_kmers_in_fasta_file_par<T, D>(
    file_name: String,
    kmer_set: &HashMap<Vec<u8>, T>,
    kmer_size: usize,
    stranded: bool,
    query_reverse: bool,
) -> anyhow::Result<(usize, usize, usize)>
where
    T: KmerCounter,
    D: MatchedSequence + Send + 'static,
{

    let reader = match file_name.as_str() {
        "" => needletail::parse_fastx_stdin()?,
        f  => needletail::parse_fastx_file(f)?,
    };

    Pipeline::<Option<D>, WithoutId>::run(
        reader,
        // map
        |record| {
            if query_reverse {
                // we need to reverse complement the sequence first
                rev_comp(record.seq);
            }
            let total_nucleotides = record.seq.len();
            let total_kmer = record.seq.len() - kmer_size + 1;

            let proxy_shared_kmers = shared_kmers_par::<_, D>(
                kmer_set,
                &record.seq,
                record.read_id,
                kmer_size,
                stranded,
                false, // in this case we map only the kmer or its reverse complement not both
            );
            let match_kmer = proxy_shared_kmers.match_count();

            *record.extra = Some(proxy_shared_kmers);

            (total_nucleotides, total_kmer, match_kmer)
        },
        // reduce
        (
            || (0, 0, 0),
            |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2)
        ),
        // writer
        NO_WRITER,
        )
}

/// count the number of indexed kmers in a given read
pub fn shared_kmers_par<C, D>(
    kmer_set: &HashMap<Vec<u8>, C>,
    read: &[u8],
    read_id: usize,
    kmer_size: usize,
    stranded: bool,
    map_both_strands: bool,
) -> D
where
    C: KmerCounter,
    D: MatchedSequence + Sized,
{
    // prints the read in human readable format
    //println!("Read {}: {}", read_id, String::from_utf8_lossy(read));
    if read.len() < kmer_size {
        return D::new(0);
    }
    let mut result = D::new(read.len() - kmer_size + 1);
    let reverse_complement = if stranded { Some(false) } else { None };

    let mut buf = [0].repeat(kmer_size);

    if !map_both_strands {
        // if we do not map both strands, we only map the kmer or its reverse complement
        let canonical_kmer = buf.as_mut_slice();
        // For computing the numbe of positions covered by at least a kmer, we need to keep track of the first uncovered position
        let mut first_uncovered_position = 0;

        for i in 0..(read.len() - kmer_size + 1) {
            let kmer = &read[i..(i + kmer_size)];
            let sequence_normalizer = SequenceNormalizer::new(kmer, reverse_complement);
            sequence_normalizer.copy_to_slice(canonical_kmer);
            if let Some(kmer_counter) = kmer_set.get(canonical_kmer) {
                result.add_match(i, sequence_normalizer.is_raw());
                if first_uncovered_position < i {
                    result.add_covered_base(kmer_size);
                }
                else {
                    result.add_covered_base(kmer_size + i - first_uncovered_position);
                }
                first_uncovered_position = i + kmer_size;


                kmer_counter.add_match(crate::kmer_counter::KmerMatch {
                    id_read: (read_id),
                    position: (i),
                    forward: (sequence_normalizer.is_raw()),
                });
            } 
        }
        result
    } else {
        // if we map both strands, we map the kmer and its reverse complement
        // Note that if --stranded is not set, the mapping is always detected in forward strand
        let normalizer_kmer = buf.as_mut_slice();
        // For computing the numbe of positions covered by at least a kmer, we need to keep track of the first uncovered position
        let mut first_uncovered_position = 0;

        for i in 0..(read.len() - kmer_size + 1) {
            let kmer = &read[i..(i + kmer_size)];
            let sequence_normalizer = SequenceNormalizer::new(kmer, reverse_complement);
            sequence_normalizer.copy_to_slice(normalizer_kmer);


            if let Some(kmer_counter) = kmer_set.get(normalizer_kmer) {
                result.add_match(i, true);
                kmer_counter.add_match(crate::kmer_counter::KmerMatch {
                    id_read: (read_id),
                    position: (i),
                    forward: true,
                });
                if first_uncovered_position < i {
                    result.add_covered_base(kmer_size);
                }
                else {
                    result.add_covered_base(kmer_size + i - first_uncovered_position);
                }
                first_uncovered_position = i + kmer_size;
            } else {
                // forward did not match, we try the reverse one
                rev_comp(normalizer_kmer);
                if let Some(kmer_counter) = kmer_set.get(normalizer_kmer) {
                    result.add_match(i, false);
                    kmer_counter.add_match(crate::kmer_counter::KmerMatch {
                        id_read: (read_id),
                        position: (i),
                        forward: false,
                    });
                    if first_uncovered_position < i {
                        result.add_covered_base(kmer_size);
                    }
                    else {
                        result.add_covered_base(kmer_size + i - first_uncovered_position);
                    }
                    first_uncovered_position = i + kmer_size;
                }   
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    /* crate use */
    use biotest::Format as _;

    /* project use */
    use super::*;
    use crate::matched_sequences;

    #[test]
    fn test_rev_comp() {
        let mut seq = b"AAACC".to_vec();
        rev_comp(&mut seq);
        assert_eq!(seq, b"GGTTT".to_vec());
    }

    #[test]
    fn shared_kmers() -> anyhow::Result<()> {
        let kmer_size = 15;
        let seq_len = 150;

        let mut rng = biotest::rand();
        let s_generator = biotest::Sequence::builder().sequence_len(seq_len).build()?;
        let k_generator = biotest::Fasta::builder().sequence_len(kmer_size).build()?;

        let mut sequence = vec![];
        s_generator.record(&mut sequence, &mut rng)?;

        let mut data = vec![];
        k_generator.records(&mut data, &mut rng, 100)?;
        data.extend(b">sequence\n");
        data.extend(&sequence);

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();
        let kmers_in_path = temp_path.join("kmers_in.fasta");

        std::fs::File::create(&kmers_in_path)?.write_all(&data)?;

        let (kmer_set_cano, _) = crate::kmer_hash::index_kmers::<atomic_counter::RelaxedCounter>(
            kmers_in_path.display().to_string(),
            kmer_size,
            false,
            false,
        )?;

        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_cano,
                &sequence,
                42,
                kmer_size,
                false,
                false
            )
            .count,
            136
        );

        let mut random_sequence = vec![];
        s_generator.record(&mut random_sequence, &mut rng)?;

        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_cano,
                &random_sequence,
                42,
                kmer_size,
                false,
                false
            )
            .count,
            0
        );

        let to_small_sequence = random_sequence[10..20].to_vec();

        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_cano,
                &to_small_sequence,
                42,
                kmer_size,
                false,
                false
            )
            .count,
            0
        );

        let (kmer_set_both, _) = crate::kmer_hash::index_kmers::<atomic_counter::RelaxedCounter>(
            kmers_in_path.display().to_string(),
            kmer_size,
            false,
            true,
        )?;

        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_both,
                &sequence,
                42,
                kmer_size,
                false,
                true
            )
            .count,
            136
        );

        rev_comp(&mut sequence);

        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_both,
                &sequence,
                42,
                kmer_size,
                true,
                true
            )
            .count,
            136
        );

        Ok(())
    }
}
