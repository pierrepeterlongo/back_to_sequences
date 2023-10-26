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
}





/// check that a file name corresponds to a non empty file:
fn validate_non_empty_file(in_file: String) {
    if let Ok(metadata) = fs::metadata(in_file.clone()) {
        // Check if the file exists
        if ! metadata.is_file() {
            panic!("{:#} exists, but it's not a file.", in_file);
        }
    } else {
        panic!("The {} file does not exist or there was an error checking its existence.", in_file);
    }
}

/// index all kmers of size kmer_size in the fasta file
/// returns a hashmap with the kmers as keys and their count as values, initialized to 0
fn index_kmers<T:Default>(file_name: String, kmer_size: usize, stranded: bool) -> Result<(HashMap<Vec<u8>, T>, usize), io::Error> {
    let mut kmer_set = HashMap::new();
    let reverse_complement = if stranded { Some(false) } else { None };

    let reader = initialize_reader(&file_name).unwrap();
    for record in reader {
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
    query_reverse: bool) -> std::io::Result<()>{

    // structure for buffering and reordering the output
    //
    // Rayon::iter::ParallelBridge() does not preserve the original order of the items. We buffer
    // the result in a hashmap and reorder them before output.
    struct Output {
        file: BufWriter<File>,
        buffer: HashMap<usize, (fxread::Record, usize)>,
        id: usize,
    }
    let output = std::sync::Mutex::new(Output{
        file: BufWriter::new(File::create(out_fasta)?),
        buffer: HashMap::new(),
        id: 0,
    });
    let output_record = |file: &mut BufWriter<_>, record: &fxread::Record, nb_shared_kmers| -> std::io::Result<()>
    {
        // round percent_shared_kmers to 3 decimals and transform to percents
        let percent_shared_kmers = round((nb_shared_kmers as f32 / (record.seq().len() - kmer_size + 1) as f32) * 100.0, 2) ;
        if percent_shared_kmers > min_threshold && percent_shared_kmers <= max_threshold { // supports the user defined thresholds

            let record_as_string = record.as_str().trim().as_bytes(); 
            let mut iter = record_as_string.split(|&x| x == b'\n');

            file.write_all(iter.next().unwrap())?;   // write the original header of the record
            write!(file, " {} {}\n", nb_shared_kmers, percent_shared_kmers)?; // append metrics
            for line in iter {
                file.write_all(line)?;
                file.write_all(b"\n")?;
            }
        } // end read contains at least one indexed kmer
        Ok(())
    };

    let (tx, rx) = std::sync::mpsc::sync_channel(1024);
    let (_, result) = rayon::join(move ||{// lance deux threads 
        if file_name.len() > 0 {
            let reader = initialize_reader(&file_name).unwrap();
            for record in reader {
                tx.send(record).unwrap();
            }
        }
        else {
            let input = stdin().lock();
            let reader = initialize_stdin_reader(input).unwrap();

            for record in reader {
                tx.send(record).unwrap();
            }
        }
    }, ||{
        rx.into_iter().enumerate().par_bridge().try_for_each(|(id, mut record)| -> std::io::Result<()>{
            if query_reverse {
                record.rev_comp();  // reverse the sequence in place
            }
            let nb_shared_kmers = count_shared_kmers_par(kmer_set, record.seq(), kmer_size, stranded);

            let mut out = output.lock().unwrap();
            out.buffer.insert(id, (record, nb_shared_kmers));
            assert!(id >= out.id);
            for i in out.id.. {
                match out.buffer.remove(&i) {
                    None => break,
                    Some((rec, nb)) => {
                        output_record(&mut out.file, &rec, nb)?;
                        out.id += 1;
                    }
                 }
            }
            Ok(())
        }) // end of for each
    }); // end of rayon join
    result?;
    assert!(output.lock().unwrap().buffer.is_empty());
    Ok(())
}

/// count the number of indexed kmers in a given read
fn count_shared_kmers_par(kmer_set:  &HashMap<Vec<u8>, atomic_counter::RelaxedCounter>, read: &[u8], kmer_size: usize, stranded: bool) -> usize {
    let mut shared_kmers_count = 0;
    let mut canonical_kmer = Vec::new();
    let reverse_complement = if stranded { Some(false) } else { None };

    for i in 0..(read.len() - kmer_size + 1) {
        let kmer = &read[i..(i + kmer_size)];
        canonical_kmer.clear();
        canonical_kmer.extend(SequenceNormalizer::new(kmer, reverse_complement).iter());
        if kmer_set.contains_key(&canonical_kmer){
            shared_kmers_count += 1;
            // kmer_set[&canonical_kmer] += 1;
            // kmer_set.insert(canonical_kmer, 1 + kmer_set[&canonical_kmer] );
            
            // *kmer_set.get_mut(&canonical_kmer).unwrap().add(1);
            kmer_set[&canonical_kmer].inc();

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
    query_reverse: bool) -> std::io::Result<()> {
      
    // check that inkmers and reads_file are non empty files:
    if in_fasta_reads.len() > 0 {
        validate_non_empty_file(in_fasta_reads.clone());
    }
    validate_non_empty_file(in_fasta_kmers.clone());

    
    match index_kmers::<RelaxedCounter>(in_fasta_kmers, kmer_size, stranded) {

        Ok((kmer_set, kmer_size)) => {
            let _ = count_kmers_in_fasta_file_par(in_fasta_reads, &kmer_set, kmer_size, out_fasta_reads.clone(), min_threshold, max_threshold, stranded, query_reverse);
            println!("Filtered sequences with exact kmer count are in file {}", out_fasta_reads);
            

            // if the out_kmers_file is not empty, we output counted kmers in the out_kmers_file file
            if out_txt_kmers.len() > 0 {
                
                // prints all kmers from kmer_set that have a count > 0
                let mut output = File::create(out_txt_kmers.clone())?;
                // let mut output = File::create(out_kmers_file);
                for (kmer, count) in kmer_set.iter() {
                    if count.get() > 0 {
                        output.write_all(kmer)?;
                        write!(output, " {}\n", count.get())?;
                    }
                }
            println!("kmers with their number of occurrences in the original sequences are in file {}", out_txt_kmers);
            }
        }

        Err(err) => eprintln!("Error indexing kmers: {}", err),
    }
    Ok(())

}
