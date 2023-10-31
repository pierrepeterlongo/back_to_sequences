use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::{BufWriter,Write, stdin};
use std::io::{self};
use std::num::NonZeroU8;
use auto_enums::auto_enum;
use fxread::{initialize_reader,initialize_stdin_reader};
use atomic_counter::{RelaxedCounter, AtomicCounter};
use rayon::prelude::*;

fn base_complement(base: u8) -> Option<NonZeroU8>
{
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        c if c < 128 => c,
        // all non-ascii character yield None so that we abord in case we try to reverse complement
        // a multibyte unicode character.
        _ => 0,
    }.try_into().ok()
}
 
/// Zero-copy object for normalizing a sequence
///
/// - converts to uppercase ascii
/// - optionally reverse complements the sequence
///
/// To avoid extra memory allocations, the original slice is kept in place and the resulting
/// sequence is returned through an iterator (`fn iter()`).
/// 
struct SequenceNormalizer<'a>
{
    raw: &'a [u8],
    reverse_complement: bool,
}

/// forward lookup-table for the SequenceNormalizer
#[ctor::ctor]
static FORWARD_MAP : [u8;256] = {
    let mut a = [0;256];
    a.iter_mut().enumerate().for_each(|(i, c)|{ *c = (i as u8).to_ascii_uppercase()});
    a
};
/// reverse lookup-table for the SequenceNormalizer
#[ctor::ctor]
static REVERSE_MAP : [Option<NonZeroU8>;256] = {
    let mut a = [None;256];
    a.iter_mut().enumerate().for_each(|(i, c)|{ *c = base_complement((i as u8).to_ascii_uppercase())});
    a
};
 
impl<'a> SequenceNormalizer<'a>
{
    /// - `raw` is the original sequence (ascii string)
    /// - `reverse_complement` may be:
    ///     - `Some(false)` to get the original sequence
    ///     - `Some(true)` to get the reverse complement
    ///     - `None` to get the canonical sequence
    fn new(raw: &'a [u8], reverse_complement: Option<bool>) -> Self
    {
        Self{raw, reverse_complement: reverse_complement.unwrap_or_else(|| {
            let forward = Self::iter_impl(raw, false);
            let reverse = Self::iter_impl(raw, true);
            reverse.cmp(forward).is_lt()
        })}
    }

    #[auto_enum(Iterator)]
    fn iter_impl(raw: &[u8] , reverse_complement: bool) -> impl Iterator<Item=u8> + '_
    {
        if reverse_complement {
            raw.iter().rev().map(|c| {
                REVERSE_MAP[*c as usize]
                    .expect("cannot complement base (contains non-ascii byte: 0x{:x})").into()})
        } else {
            raw.iter().map(|c| FORWARD_MAP[*c as usize])
        }
    }

    /// Get an iterator on the normalized sequence
    fn iter(&self) -> impl Iterator<Item=u8> + '_
    {
        Self::iter_impl(self.raw, self.reverse_complement)
    }

    /// Copy the normalized sequence into a slice
    ///
    /// panics for the slice has a different length
    fn copy_to_slice(&self, dest: &mut [u8])
    {
        assert_eq!(dest.len(), self.raw.len());
        for (i, c) in self.iter().enumerate() {
            dest[i] = c;
        }
    }
}





/// check that a file name corresponds to a non empty file:
fn validate_non_empty_file(in_file: String) -> Result<(), ()> {
    if let Ok(metadata) = fs::metadata(in_file.clone()) {
        // Check if the file exists
        if ! metadata.is_file() {
            return Err(eprintln!("{:#} exists, but it's not a file.", in_file));
        }
    } else {
        return Err(eprintln!("The {} file does not exist or there was an error checking its existence.", in_file));
    }
    Ok(())
}

/// index all kmers of size kmer_size in the fasta file
/// returns a hashmap with the kmers as keys and their count as values, initialized to 0
fn index_kmers<T:Default>(file_name: String, kmer_size: usize, stranded: bool) -> anyhow::Result<(HashMap<Vec<u8>, T>, usize)> {
    let mut kmer_set = HashMap::new();
    let reverse_complement = if stranded { Some(false) } else { None };

    let mut reader = initialize_reader(&file_name)?;
    loop {
        let Some(record) = reader.next_record()? else { break };
        let acgt_sequence = record.seq();
        // for each kmer of the sequence, insert it in the kmer_set
        for i in 0..(acgt_sequence.len() - kmer_size + 1) {
            let kmer = &acgt_sequence[i..(i + kmer_size)];
            kmer_set.insert(
                SequenceNormalizer::new(kmer, reverse_complement).iter().collect(),
                Default::default(), // RelaxedCounter::new(0)
                );
        }
    }
    println!("Indexed {} kmers, each of size {}", kmer_set.len(), kmer_size);
    
    Ok((kmer_set, kmer_size))
}

/// round a float to a given number of decimals
fn round(x: f32, decimals: u32) -> f32 {
    let y = 10i32.pow(decimals) as f32;
    (x * y).round() / y
}

/// for each sequence of a given fasta file, count the number of indexed kmers it contains
/// and output the sequence if its ratio of indexed kmers is in ]min_threshold, max_threshold]
fn count_kmers_in_fasta_file_par(file_name: String, 
    kmer_set:  &HashMap<Vec<u8>, atomic_counter::RelaxedCounter>, 
    kmer_size: usize, 
    out_fasta: String, 
    min_threshold: f32, 
    max_threshold: f32,
    stranded: bool, 
    query_reverse: bool) -> Result<(), ()>{

    const CHUNK_SIZE: usize = 32; // number of records
    const INPUT_CHANNEL_SIZE:  usize = 8; // in units of CHUNK_SIZE records
    const OUTPUT_CHANNEL_SIZE: usize = 8; // in units of CHUNK_SIZE records

    struct Chunk {
        id: usize,
        records: Vec<(fxread::Record, Option<usize>)>,
    }

    let (input_tx,   input_rx) = std::sync::mpsc::sync_channel::<Chunk>(INPUT_CHANNEL_SIZE);
    let (output_tx, output_rx) = std::sync::mpsc::sync_channel::<Chunk>(OUTPUT_CHANNEL_SIZE);
    
    let mut output_file = BufWriter::new(
        File::create(out_fasta).map_err(
            |e| eprintln!("Error: failed to open the sequence file for writing: {}", e))?);

    let mut output_record = move |(record, nb_shared_kmers): (fxread::Record, Option<usize>)| -> std::io::Result<()> {
        // round percent_shared_kmers to 3 decimals and transform to percents
        let percent_shared_kmers = round((nb_shared_kmers.unwrap() as f32 / (record.seq().len() - kmer_size + 1) as f32) * 100.0, 2) ;
        if percent_shared_kmers > min_threshold && percent_shared_kmers <= max_threshold { // supports the user defined thresholds

            let record_as_string = record.as_str().trim().as_bytes(); 
            let mut iter = record_as_string.split(|&x| x == b'\n');

            output_file.write_all(iter.next().unwrap())?;   // write the original header of the record
            write!(output_file, " {} {}\n", nb_shared_kmers.unwrap(), percent_shared_kmers)?; // append metrics
            for line in iter {
                output_file.write_all(line)?;
                output_file.write_all(b"\n")?;
            }
        } // end read contains at least one indexed kmer
        Ok(())
    };


    let reader_thread = std::thread::spawn(move || -> anyhow::Result<()> {
        let mut reader = 
            if file_name.len() > 0 { initialize_reader(&file_name).unwrap() } 
            else { initialize_stdin_reader(stdin().lock()).unwrap() };

        for id in 0.. {
            let mut vec = Vec::with_capacity(CHUNK_SIZE);
            for _ in 0..CHUNK_SIZE {
                match reader.next_record()? {
                    None => break,
                    Some(record) => vec.push((record, None)),
                }
            }
            if vec.is_empty() || input_tx.send(Chunk{ id, records: vec }).is_err() {
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
                        },
                        Ok(chunk) =>
                            if chunk.id == id {
                                break chunk.records;
                            } else {
                                buffer.insert(chunk.id, chunk.records);
                            }
                    }
                },
            };
            records.into_iter().try_for_each(&mut output_record)?;
        }
        unreachable!()
    });

    input_rx.into_iter().par_bridge().try_for_each(move |mut chunk| {
        for (record, nb_shared_kmers) in &mut chunk.records {
            if query_reverse {
                record.rev_comp();  // reverse the sequence in place
            }
            *nb_shared_kmers = Some(count_shared_kmers_par(kmer_set, record.seq(), kmer_size, stranded));
        }
        output_tx.send(chunk)
    })
    .ok(); // result ignored because an error on output_tx.send() would come together with an io
           // error in the writer thread

    reader_thread.join().unwrap()
        .map_err(|e| eprintln!("Error reading the sequences: {}", e))?;
    writer_thread.join().unwrap()
        .map_err(|e| eprintln!("Error writing the sequences: {}", e))
}

/// count the number of indexed kmers in a given read
fn count_shared_kmers_par(kmer_set:  &HashMap<Vec<u8>, atomic_counter::RelaxedCounter>, read: &[u8], kmer_size: usize, stranded: bool) -> usize {
    let mut shared_kmers_count = 0;
    let reverse_complement = if stranded { Some(false) } else { None };

    let mut buf = [0].repeat(kmer_size);
    let canonical_kmer = buf.as_mut_slice();

    for i in 0..(read.len() - kmer_size + 1) {
        let kmer = &read[i..(i + kmer_size)];
        SequenceNormalizer::new(kmer, reverse_complement).copy_to_slice(canonical_kmer);
        if let Some(kmer_counter) = kmer_set.get(canonical_kmer)
        {
            shared_kmers_count += 1;
            // kmer_set[&canonical_kmer] += 1;
            // kmer_set.insert(canonical_kmer, 1 + kmer_set[&canonical_kmer] );
            
            // *kmer_set.get_mut(&canonical_kmer).unwrap().add(1);
            kmer_counter.inc();
        }
    }
    shared_kmers_count
}


/// Extract sequences that contain some kmers and 
/// output the kmers that occur in the reads with their number of occurrences
pub fn back_to_sequences(in_fasta_reads: String, 
    in_fasta_kmers: String, 
    out_fasta_reads:String, 
    out_txt_kmers: String, 
    kmer_size: usize, 
    min_threshold: f32, 
    max_threshold: f32,
    stranded: bool,
    query_reverse: bool) -> Result<(),()> {
      
    
    // check that in_fasta_reads is a non empty file if it exists:
    if in_fasta_reads.len() > 0 {
        validate_non_empty_file(in_fasta_reads.clone())?;
    }
    validate_non_empty_file(in_fasta_kmers.clone())?;
    // check that in_fasta_kmers is a non empty file:

    let (kmer_set, kmer_size) = index_kmers::<RelaxedCounter>(in_fasta_kmers, kmer_size, stranded)
        .map_err(|e| eprintln!("Error indexing kmers: {}", e))?;

    count_kmers_in_fasta_file_par(in_fasta_reads, &kmer_set, kmer_size, out_fasta_reads.clone(), min_threshold, max_threshold, stranded, query_reverse)?;
    println!("Filtered sequences with exact kmer count are in file {}", out_fasta_reads);


    // if the out_kmers_file is not empty, we output counted kmers in the out_kmers_file file
    if out_txt_kmers.len() > 0 {
        (|| -> io::Result<_> {
            // prints all kmers from kmer_set that have a count > 0
            let mut output = File::create(&out_txt_kmers)?;
            for (kmer, count) in kmer_set.iter() {
                if count.get() > 0 {
                    output.write_all(kmer)?;
                    write!(output, " {}\n", count.get())?;
                }
            }
            Ok(())
        })().map_err(|e| eprintln!("Error writing the kmers file: {}", e))?;

        println!("kmers with their number of occurrences in the original sequences are in file {}", out_txt_kmers);
    }
    Ok(())
}
