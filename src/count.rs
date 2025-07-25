//! Count declarations

/* std use */
use std::fs::File;
use std::io::{self};
use std::io::{stdin, BufWriter, Write};
// use std::sync::Mutex;

/* crates use */
use ahash::AHashMap as HashMap;
use anyhow::Context as _;
use fxread::{initialize_reader, initialize_stdin_reader};
use rayon::prelude::*;

/* project use */
use crate::kmer_counter::KmerCounter;
use crate::matched_sequences::MatchedSequence;
use crate::sequence_normalizer::SequenceNormalizer;
use crate::kmer_prefiltration::KmerPrefiltration;

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
    prefilter: &KmerPrefiltration,
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
    const CHUNK_SIZE: usize = 32; // number of records
    const INPUT_CHANNEL_SIZE: usize = 8; // in units of CHUNK_SIZE records
    const OUTPUT_CHANNEL_SIZE: usize = 8; // in units of CHUNK_SIZE records

    struct Chunk<D> {
        id: usize,
        records: Vec<(usize, fxread::Record, Option<D>)>,
        // records: Vec<(usize, fxread::Record, Box<dyn MatchedSequence>)>,
    }

    let (input_tx, input_rx) = std::sync::mpsc::sync_channel::<Chunk<D>>(INPUT_CHANNEL_SIZE);
    let (output_tx, output_rx) = std::sync::mpsc::sync_channel::<Chunk<D>>(OUTPUT_CHANNEL_SIZE);

    let mut output_file = BufWriter::new(
        File::create(out_fasta).context("Error: failed to open the sequence file for writing")?,
    );

    let mut output_record = move |(_read_id, record, matched_sequence): (
        usize,
        fxread::Record,
        Option<D>,
    )|
          -> std::io::Result<()> {
        let percent_shared_kmers = matched_sequence.as_ref().unwrap().percent_shared_kmers();

        if percent_shared_kmers > min_threshold && percent_shared_kmers <= max_threshold {
            // supports the user defined thresholds

            let record_as_string = record.as_str().trim().as_bytes();
            let mut iter = record_as_string.split(|&x| x == b'\n');

            output_file.write_all(iter.next().unwrap())?; // write the original header of the record
            writeln!(output_file, "{}", matched_sequence.as_ref().unwrap())?; // append metrics
            for line in iter {
                output_file.write_all(line)?;
                output_file.write_all(b"\n")?;
            }
        } // end read contains at least one indexed kmer
        Ok(())
    };

    let reader_thread = std::thread::spawn(move || -> anyhow::Result<()> {
        let mut reader = if file_name.is_empty() {
            initialize_stdin_reader(stdin().lock()).unwrap()
        } else {
            initialize_reader(&file_name).unwrap()
        };

        let mut read_id = 0;
        for id in 0.. {
            let mut vec = Vec::with_capacity(CHUNK_SIZE);
            for _ in 0..CHUNK_SIZE {
                match reader.next_record()? {
                    None => break,
                    Some(record) => vec.push((read_id, record, None)),
                }
                read_id += 1;
            }
            if vec.is_empty() || input_tx.send(Chunk { id, records: vec }).is_err() {
                return Ok(());
            }
        }
        unreachable!()
    });

    let writer_thread = std::thread::spawn(move || -> io::Result<_> {
        // buffer for reordering the output (because Rayon::iter::ParallelBridge() does not
        // preserve the original order of the items)
        let mut buffer = HashMap::<usize, Vec<_>>::new();

        for id in 0.. {
            let records: Vec<(usize, fxread::Record, Option<D>)> = match buffer.remove(&id) {
                Some(vec) => vec,
                None => loop {
                    match output_rx.recv() {
                        Err(_) => {
                            assert!(buffer.is_empty());
                            return Ok(());
                        }
                        Ok(chunk) => {
                            if chunk.id == id {
                                break chunk.records;
                            } else {
                                buffer.insert(chunk.id, chunk.records);
                            }
                        }
                    }
                },
            };
            records.into_iter().try_for_each(&mut output_record)?;
        }
        unreachable!()
    });

    let (total_nucleotides, total_kmer, match_kmer) = input_rx
        .into_iter()
        .par_bridge()
        .map(move |mut chunk| {
            let mut tt_kmer = 0;
            let mut match_kmer = 0;
            let mut tt_nt = 0; // total number of nucleotides

            for (read_id, record, nb_shared_kmers) in &mut chunk.records {
                tt_nt += record.seq().len();
                tt_kmer += record.seq().len() - kmer_size + 1;
                record.upper();
                if query_reverse {
                    record.rev_comp(); // reverse the sequence in place
                }
                let proxy_shared_kmers = shared_kmers_par::<_, D>(
                    kmer_set,
                    prefilter,
                    record.seq(),
                    *read_id,
                    kmer_size,
                    stranded,
                    map_both_strands,
                );

                match_kmer += proxy_shared_kmers.match_count();

                *nb_shared_kmers = Some(proxy_shared_kmers);
            }
            output_tx.send(chunk).ok(); // result ignored because an error on output_tx.send() would come together with an io
                                        // error in the writer thread
            (tt_nt, tt_kmer, match_kmer)
        })
        .reduce(|| (0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

    reader_thread
        .join()
        .unwrap()
        .context("Error reading the sequences")?;
    writer_thread
        .join()
        .unwrap()
        .context("Error writing the sequences")?;

    Ok((total_nucleotides, total_kmer, match_kmer))
}

/// for each sequence of a given fasta file, count the number of indexed kmers it contains
pub fn only_kmers_in_fasta_file_par<T, D>(
    file_name: String,
    kmer_set: &HashMap<Vec<u8>, T>,
    prefilter: &KmerPrefiltration,
    kmer_size: usize,
    stranded: bool,
    query_reverse: bool,
) -> anyhow::Result<(usize, usize, usize)>
where
    T: KmerCounter,
    D: MatchedSequence + Send + 'static,
{
    const CHUNK_SIZE: usize = 32; // number of records
    const INPUT_CHANNEL_SIZE: usize = 8; // in units of CHUNK_SIZE records
    const OUTPUT_CHANNEL_SIZE: usize = 8; // in units of CHUNK_SIZE records

    struct Chunk<D> {
        records: Vec<(usize, fxread::Record, Option<D>)>,
        // records: Vec<(usize, fxread::Record, Box<dyn MatchedSequence>)>,
    }

    let (input_tx, input_rx) = std::sync::mpsc::sync_channel::<Chunk<D>>(INPUT_CHANNEL_SIZE);
    let (_output_tx, _output_rx) = std::sync::mpsc::sync_channel::<Chunk<D>>(OUTPUT_CHANNEL_SIZE);

    let reader_thread = std::thread::spawn(move || -> anyhow::Result<()> {
        let mut reader = if file_name.is_empty() {
            initialize_stdin_reader(stdin().lock()).unwrap()
        } else {
            initialize_reader(&file_name).unwrap()
        };

        let mut read_id = 0;
        for _ in 0.. {
            let mut vec = Vec::with_capacity(CHUNK_SIZE);
            for _ in 0..CHUNK_SIZE {
                match reader.next_record()? {
                    None => break,
                    Some(record) => vec.push((read_id, record, None)),
                }
                read_id += 1;
            }

            if vec.is_empty() || input_tx.send(Chunk { records: vec }).is_err() {
                return Ok(());
            }
        }
        unreachable!()
    });

    let (total_nucleotides, total_kmer, match_kmer) = input_rx
        .into_iter()
        .par_bridge()
        .map(move |mut chunk| {
            let mut tt_nt = 0; // total number of nucleotides
            let mut tt_kmer = 0;
            let mut match_kmer = 0;

            for (read_id, record, nb_shared_kmers) in &mut chunk.records {
                tt_nt += record.seq().len();
                tt_kmer += record.seq().len() - kmer_size + 1;
                record.upper();
                if query_reverse {
                    record.rev_comp(); // reverse the sequence in place
                }
                let proxy_shared_kmers = shared_kmers_par::<_, D>(
                    kmer_set,
                    prefilter,
                    record.seq(),
                    *read_id,
                    kmer_size,
                    stranded,
                    false, // in this case we map only the kmer or its reverse complement not both
                );
                match_kmer += proxy_shared_kmers.match_count();

                *nb_shared_kmers = Some(proxy_shared_kmers);
            }

            (tt_nt, tt_kmer, match_kmer)
        })
        .reduce(|| (0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

    reader_thread
        .join()
        .unwrap()
        .context("Error reading the sequences")?;

    Ok((total_nucleotides, total_kmer, match_kmer))
}

/// count the number of indexed kmers in a given read
pub fn shared_kmers_par<C, D>(
    kmer_set: &HashMap<Vec<u8>, C>,
    prefilter: &KmerPrefiltration,
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
        
        let positions_worse_to_test = prefilter.potiential_kmer_positions(read);
        for &i in &positions_worse_to_test { // OPTIMIZED
        // for i in 0..read.len()-kmer_size+1 { // NOT OPTIMIZED
            let kmer = &read[i..(i + kmer_size)];
            let sequence_normalizer = SequenceNormalizer::new(kmer, reverse_complement);
            sequence_normalizer.copy_to_slice(canonical_kmer);
            if let Some(kmer_counter) = kmer_set.get(canonical_kmer) {
                ///// DBUG START ///////
                // if !positions_worse_to_test.contains(&i) {
                //     panic!("kmer {:?}, pos {} (seq len {}) is in the read but not in the positions worse to test", 
                //     str::from_utf8(kmer), 
                //     i, 
                //     read.len()
                // );
                // }
                ///// DBUG ENDS ////////
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
        let positions_worse_to_test = prefilter.potiential_kmer_positions(read);
        for i in positions_worse_to_test {
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
        let prefilter = KmerPrefiltration::from_kmer_set(
            kmer_set_cano.keys().cloned().collect::<Vec<_>>().as_slice(),
            0.1,
            15,
            11,
        );
        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_cano,
                &prefilter,
                &sequence,
                42,
                kmer_size,
                false,
                false
            )
            .count,
            135
        );

        let mut random_sequence = vec![];
        s_generator.record(&mut random_sequence, &mut rng)?;

        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_cano,
                &prefilter,
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
                &prefilter,
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
                &prefilter,
                &sequence,
                42,
                kmer_size,
                false,
                true
            )
            .count,
            135
        );

        rev_comp(&mut sequence);

        assert_eq!(
            shared_kmers_par::<_, matched_sequences::MachedCount>(
                &kmer_set_both,
                &prefilter,
                &sequence,
                42,
                kmer_size,
                true,
                true
            )
            .count,
            135
        );

        Ok(())
    }
}