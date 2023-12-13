//! Count declarations

/* std use */
use std::fs::File;
use std::io::{self};
use std::io::{stdin, BufWriter, Write};

/* crates use */
use ahash::AHashMap as HashMap;
use fxread::{initialize_reader, initialize_stdin_reader};
use rayon::prelude::*;

use atomic_counter::AtomicCounter as _;

/* project use */
use crate::sequence_normalizer::SequenceNormalizer;

/// round a float to a given number of decimals
fn round(x: f32, decimals: u32) -> f32 {
    let y = 10i32.pow(decimals) as f32;
    (x * y).round() / y
}

/// for each sequence of a given fasta file, count the number of indexed kmers it contains
/// and output the sequence if its ratio of indexed kmers is in ]min_threshold, max_threshold]
pub fn kmers_in_fasta_file_par(
    file_name: String,
    kmer_set: &HashMap<Vec<u8>, atomic_counter::RelaxedCounter>,
    kmer_size: usize,
    out_fasta: String,
    min_threshold: f32,
    max_threshold: f32,
    stranded: bool,
    query_reverse: bool,
) -> Result<(), ()> {
    const CHUNK_SIZE: usize = 32; // number of records
    const INPUT_CHANNEL_SIZE: usize = 8; // in units of CHUNK_SIZE records
    const OUTPUT_CHANNEL_SIZE: usize = 8; // in units of CHUNK_SIZE records

    struct Chunk {
        id: usize,
        records: Vec<(fxread::Record, Option<usize>)>,
    }

    let (input_tx, input_rx) = std::sync::mpsc::sync_channel::<Chunk>(INPUT_CHANNEL_SIZE);
    let (output_tx, output_rx) = std::sync::mpsc::sync_channel::<Chunk>(OUTPUT_CHANNEL_SIZE);

    let mut output_file =
        BufWriter::new(File::create(out_fasta).map_err(|e| {
            eprintln!("Error: failed to open the sequence file for writing: {}", e)
        })?);

    let mut output_record =
        move |(record, nb_shared_kmers): (fxread::Record, Option<usize>)| -> std::io::Result<()> {
            // round percent_shared_kmers to 3 decimals and transform to percents
            let percent_shared_kmers = round(
                (nb_shared_kmers.unwrap() as f32 / (record.seq().len() - kmer_size + 1) as f32)
                    * 100.0,
                2,
            );
            if percent_shared_kmers > min_threshold && percent_shared_kmers <= max_threshold {
                // supports the user defined thresholds

                let record_as_string = record.as_str().trim().as_bytes();
                let mut iter = record_as_string.split(|&x| x == b'\n');

                output_file.write_all(iter.next().unwrap())?; // write the original header of the record
                writeln!(
                    output_file,
                    " {} {}",
                    nb_shared_kmers.unwrap(),
                    percent_shared_kmers
                )?; // append metrics
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

        for id in 0.. {
            let mut vec = Vec::with_capacity(CHUNK_SIZE);
            for _ in 0..CHUNK_SIZE {
                match reader.next_record()? {
                    None => break,
                    Some(record) => vec.push((record, None)),
                }
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
            let records = match buffer.remove(&id) {
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

    input_rx
        .into_iter()
        .par_bridge()
        .try_for_each(move |mut chunk| {
            for (record, nb_shared_kmers) in &mut chunk.records {
                record.upper();
                if query_reverse {
                    record.rev_comp(); // reverse the sequence in place
                }
                *nb_shared_kmers = Some(shared_kmers_par(
                    kmer_set,
                    record.seq(),
                    kmer_size,
                    stranded,
                ));
            }
            output_tx.send(chunk)
        })
        .ok(); // result ignored because an error on output_tx.send() would come together with an io
               // error in the writer thread

    reader_thread
        .join()
        .unwrap()
        .map_err(|e| eprintln!("Error reading the sequences: {}", e))?;
    writer_thread
        .join()
        .unwrap()
        .map_err(|e| eprintln!("Error writing the sequences: {}", e))
}

/// count the number of indexed kmers in a given read
pub fn shared_kmers_par(
    kmer_set: &HashMap<Vec<u8>, atomic_counter::RelaxedCounter>,
    read: &[u8],
    kmer_size: usize,
    stranded: bool,
) -> usize {
    let mut shared_kmers_count = 0;
    let reverse_complement = if stranded { Some(false) } else { None };

    let mut buf = [0].repeat(kmer_size);
    let canonical_kmer = buf.as_mut_slice();

    if read.len() < kmer_size {
        return 0;
    }
    for i in 0..(read.len() - kmer_size + 1) {
        let kmer = &read[i..(i + kmer_size)];
        SequenceNormalizer::new(kmer, reverse_complement).copy_to_slice(canonical_kmer);
        if let Some(kmer_counter) = kmer_set.get(canonical_kmer) {
            shared_kmers_count += 1;
            kmer_counter.inc();
        }
    }
    shared_kmers_count
}
