//! Kmer hash declarations

/* std use */
use std::io::Write as _;

/* crates use */
use ahash::AHashMap as HashMap;
use entropy::shannon_entropy;
use fxread::{initialize_reader, initialize_stdin_reader};
use rayon::prelude::*;

use atomic_counter::AtomicCounter as _;

/* project use */
use crate::sequence_normalizer::SequenceNormalizer;

/// given a kmer as a &[u8] check that it contains only ACGT letters
/// return true if it is the case, false otherwise
fn is_acgt(kmer: &[u8]) -> bool {
    for &byte in kmer {
        if byte != b'A' && byte != b'C' && byte != b'G' && byte != b'T' {
            return false;
        }
    }
    true
}

/// round a float to a given number of decimals
fn round(x: f32, decimals: u32) -> f32 {
    let y = 10i32.pow(decimals) as f32;
    (x * y).round() / y
}

/// Struct that store KmerCount
pub struct KmerCount<T> {
    inner: ahash::AHashMap<Vec<u8>, T>,
    kmer_size: usize,
}

impl<T> KmerCount<T> {
    /// Intialize a new KmerCounter
    pub fn from_path<P>(filename: P, kmer_size: usize, stranded: bool) -> anyhow::Result<Self>
    where
        P: AsRef<std::path::Path>,
        T: std::default::Default,
    {
        let mut kmer_count = ahash::AHashMap::new();
        let reverse_complement = if stranded { Some(false) } else { None };

        let mut reader = initialize_reader(filename.as_ref().to_str().unwrap())?;
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
            for i in 0..(acgt_sequence.len() - kmer_size + 1) {
                let kmer = &acgt_sequence[i..(i + kmer_size)];
                if is_acgt(kmer) {
                    //TODO if not: get the position of the last non acgt letter and jump to the next potential possible kmer
                    // If the entropy is too low, the kmer is not inserted
                    if shannon_entropy(kmer) < 1.0 {
                        // TODO: make this an option.
                        continue;
                    }
                    kmer_count.insert(
                        SequenceNormalizer::new(kmer, reverse_complement)
                            .iter()
                            .collect(),
                        Default::default(), // RelaxedCounter::new(0)
                    );
                }
            }
        }

        println!(
            "Indexed {} kmers, each of size {}",
            kmer_count.len(),
            kmer_size
        );

        Ok(KmerCount {
            inner: kmer_count,
            kmer_size,
        })
    }
}

impl<T> KmerCount<T>
where
    T: atomic_counter::AtomicCounter,
{
    /// for each sequence of a given fasta file, count the number of indexed kmers it contains
    /// and output the sequence if its ratio of indexed kmers is in ]min_threshold, max_threshold]
    pub fn filter_path<P>(
        &self,
        file_name: Option<P>,
        out_fasta: P,
        min_threshold: f32,
        max_threshold: f32,
        stranded: bool,
        query_reverse: bool,
    ) -> Result<(), ()>
    where
        P: AsRef<std::path::Path> + std::marker::Send + std::clone::Clone,
    {
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
            std::io::BufWriter::new(std::fs::File::create(out_fasta).map_err(|e| {
                eprintln!("Error: failed to open the sequence file for writing: {}", e)
            })?);

        let k = self.kmer_size;
        let mut output_record = move |(record, nb_shared_kmers): (
            fxread::Record,
            Option<usize>,
        )|
              -> std::io::Result<()> {
            // round percent_shared_kmers to 3 decimals and transform to percents
            let percent_shared_kmers = round(
                (nb_shared_kmers.unwrap() as f32 / (record.seq().len() - k + 1) as f32) * 100.0,
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

        let file_name_str = file_name.map(|x| x.as_ref().to_str().unwrap().to_string());
        let reader_thread = std::thread::spawn(move || -> anyhow::Result<()> {
            let mut reader = if let Some(path) = file_name_str {
                initialize_reader(&path).unwrap()
            } else {
                initialize_stdin_reader(std::io::stdin().lock()).unwrap()
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

        let writer_thread = std::thread::spawn(move || -> std::io::Result<_> {
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
                    *nb_shared_kmers = Some(self.count_kmer(record.seq(), stranded));
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
            .map_err(|e| eprintln!("Error writing the sequences: {}", e))?;

        Ok(())
    }

    /// count the number of indexed kmers in a given read
    fn count_kmer(&self, read: &[u8], stranded: bool) -> usize {
        let mut shared_kmers_count = 0;
        let reverse_complement = if stranded { Some(false) } else { None };

        let mut buf = [0].repeat(self.kmer_size);
        let canonical_kmer = buf.as_mut_slice();

        if read.len() < self.kmer_size {
            return 0;
        }

        for i in 0..(read.len() - self.kmer_size + 1) {
            let kmer = &read[i..(i + self.kmer_size)];
            SequenceNormalizer::new(kmer, reverse_complement).copy_to_slice(canonical_kmer);
            if let Some(count) = self.as_ref().get(canonical_kmer) {
                shared_kmers_count += 1;
                // kmer_set[&canonical_kmer] += 1;
                // kmer_set.insert(canonical_kmer, 1 + kmer_set[&canonical_kmer] );

                // *kmer_set.get_mut(&canonical_kmer).unwrap().add(1);
                count.inc();
            }
        }
        shared_kmers_count
    }
}

impl KmerCount<atomic_counter::RelaxedCounter> {
    /// Write in kmer_count in csv format, by default separator is ','
    pub fn to_csv<P>(&self, path: P, separator: Option<u8>) -> std::io::Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        let mut output = std::fs::File::create(path)?;

        for (kmer, count) in self.as_ref().iter() {
            if count.get() > 0 {
                output.write_all(kmer)?;
                output.write_all(&[separator.unwrap_or(b',')])?;
                writeln!(output, "{}", count.get())?;
            }
        }

        Ok(())
    }
}

impl<T> AsRef<ahash::AHashMap<Vec<u8>, T>> for KmerCount<T> {
    fn as_ref(&self) -> &ahash::AHashMap<Vec<u8>, T> {
        &self.inner
    }
}

impl<T> AsMut<ahash::AHashMap<Vec<u8>, T>> for KmerCount<T> {
    fn as_mut(&mut self) -> &mut ahash::AHashMap<Vec<u8>, T> {
        &mut self.inner
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use atomic_counter::AtomicCounter;

    use std::convert::AsRef as _;
    use std::io::Read as _;

    const FASTA_FILE: &[u8] = b">random_seq 0
GGACTGAGGAGGATACACCTATGGA
>random_seq 1
AACGCGCGGCAGCAGGTCCATTACT
>random_seq 2
GGCAGTAACTTTACATGCAGATACT
>random_seq 3
TTTATGGTGAGATTGACTTAACGGT
>random_seq 4
ATTGCGACGCTGCCCCCTCTGCGTA
>random_seq 5
TATCGCTGAAAGGCTCTTTGATCGA
>random_seq 6
GGTGCCCGTCAAGTTTATTTTTTCT
>random_seq 7
CTGCAATTAGAGAATCTATCGTCTC
>random_seq 8
GGTGCATTCATAGTAGCTGCAGAAA
>random_seq 9
ATGTTCGGACATGCCCCTATTGGGG
";

    #[test]
    fn canonical_from_path() -> anyhow::Result<()> {
        let mut input = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        input.write_all(FASTA_FILE)?;

        let kmer_count: KmerCount<u64> = KmerCount::<u64>::from_path(input.path(), 5, false)?;
        let mut kmers: Vec<Vec<u8>> = kmer_count.as_ref().keys().map(|x| x.to_vec()).collect();
        kmers.sort();

        assert_eq!(
            kmers,
            vec![
                b"AAACT", b"AAAGT", b"AACAT", b"AACGC", b"AACGG", b"AACTT", b"AAGGC", b"AAGTC",
                b"AATAG", b"AATCT", b"AATGC", b"AATGG", b"AATTG", b"ACATG", b"ACCAT", b"ACCGT",
                b"ACCTA", b"ACCTG", b"ACGAT", b"ACGCA", b"ACGCG", b"ACGCT", b"ACGGG", b"ACTAT",
                b"ACTGA", b"ACTGC", b"ACTTA", b"ACTTG", b"AGAAT", b"AGACG", b"AGAGC", b"AGATA",
                b"AGCAG", b"AGCCT", b"AGCGA", b"AGCTA", b"AGCTG", b"AGGAT", b"AGGTC", b"AGGTG",
                b"AGTAA", b"AGTAG", b"AGTAT", b"AGTCA", b"AGTCC", b"AGTTA", b"ATACA", b"ATAGA",
                b"ATAGG", b"ATCAA", b"ATCGA", b"ATCGC", b"ATCTA", b"ATCTC", b"ATCTG", b"ATGAA",
                b"ATGCA", b"ATGCC", b"ATGGA", b"ATGTA", b"ATGTC", b"ATTAC", b"ATTAG", b"ATTCA",
                b"ATTGA", b"ATTGC", b"ATTGG", b"CAAAG", b"CAATA", b"CAATC", b"CAGAA", b"CAGAG",
                b"CAGCA", b"CAGCG", b"CAGTA", b"CAGTC", b"CATAA", b"CATAG", b"CATGC", b"CATTA",
                b"CATTC", b"CCATA", b"CCCTA", b"CCGAA", b"CCGTC", b"CCTCA", b"CCTGC", b"CGAAC",
                b"CGACG", b"CGATA", b"CGATC", b"CGCAA", b"CGCAG", b"CGGAC", b"CGGCA", b"CGTCA",
                b"CGTTA", b"CTATC", b"CTCAC", b"CTCAG", b"CTCTA", b"CTGAA", b"CTGCA", b"CTGCC",
                b"CTTAA", b"CTTGA", b"CTTTA", b"GAACA", b"GAATC", b"GACGA", b"GACGC", b"GAGAC",
                b"GAGCC", b"GATAC", b"GATCA", b"GCACC", b"GCAGA", b"GCAGC", b"GCGAC", b"GCGTA",
                b"GCTAC", b"GCTGA", b"GGACA", b"GGACC", b"GGATA", b"GGCAC", b"GGGCA", b"GGTGA",
                b"GTAAA", b"GTAAC", b"GTCAA", b"GTCCA", b"GTGCA", b"GTGTA", b"GTTAA", b"GTTTA",
                b"TACTA", b"TATGA", b"TCAAA", b"TCCGA", b"TCGCA", b"TCTAA", b"TCTCA", b"TGAAA",
                b"TGCAA", b"TGTAA"
            ]
        );

        let counts: Vec<u64> = kmer_count.as_ref().values().cloned().collect();
        assert_eq!(counts, vec![0; counts.len()]);
        Ok(())
    }

    #[test]
    fn stranded_from_path() -> anyhow::Result<()> {
        let mut input = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        input.write_all(FASTA_FILE)?;

        let kmer_count: KmerCount<u64> = KmerCount::<u64>::from_path(input.path(), 5, true)?;
        let mut kmers: Vec<Vec<u8>> = kmer_count.as_ref().keys().map(|x| x.to_vec()).collect();
        kmers.sort();

        assert_eq!(
            kmers,
            vec![
                b"AACGC", b"AACGG", b"AACTT", b"AAGGC", b"AAGTT", b"AATCT", b"ACATG", b"ACCTA",
                b"ACGCG", b"ACGCT", b"ACGGT", b"ACTGA", b"ACTTA", b"ACTTT", b"AGAAT", b"AGATA",
                b"AGATT", b"AGCAG", b"AGCTG", b"AGGAT", b"AGGCT", b"AGGTC", b"AGTAA", b"AGTAG",
                b"AGTTT", b"ATACA", b"ATACT", b"ATAGT", b"ATCGA", b"ATCGC", b"ATCGT", b"ATCTA",
                b"ATGCA", b"ATGCC", b"ATGGA", b"ATGGT", b"ATGTT", b"ATTAC", b"ATTAG", b"ATTCA",
                b"ATTGA", b"ATTGC", b"ATTGG", b"CAAGT", b"CAATT", b"CACCT", b"CAGAA", b"CAGAT",
                b"CAGCA", b"CAGGT", b"CAGTA", b"CATAG", b"CATGC", b"CATTA", b"CATTC", b"CCATT",
                b"CCCGT", b"CCCTA", b"CCGTC", b"CCTAT", b"CGACG", b"CGCTG", b"CGGAC", b"CGGCA",
                b"CGTCA", b"CGTCT", b"CTATC", b"CTATG", b"CTATT", b"CTCTG", b"CTGAA", b"CTGAG",
                b"CTGCA", b"CTGCC", b"CTGCG", b"CTTAA", b"CTTTA", b"CTTTG", b"GAATC", b"GACAT",
                b"GACGC", b"GACTG", b"GACTT", b"GAGAT", b"GATAC", b"GATCG", b"GATTG", b"GCAAT",
                b"GCAGA", b"GCAGC", b"GCAGG", b"GCAGT", b"GCATT", b"GCGAC", b"GCGTA", b"GCTCT",
                b"GCTGA", b"GCTGC", b"GGACA", b"GGACT", b"GGATA", b"GGCAG", b"GGCTC", b"GGTCC",
                b"GGTGA", b"GGTGC", b"GTAAC", b"GTAGC", b"GTCAA", b"GTCCA", b"GTCTC", b"GTGAG",
                b"GTGCA", b"GTGCC", b"GTTCG", b"GTTTA", b"TAACG", b"TAACT", b"TACAC", b"TACAT",
                b"TAGAG", b"TAGCT", b"TAGTA", b"TATCG", b"TATGG", b"TATTG", b"TCAAG", b"TCATA",
                b"TCCAT", b"TCGCT", b"TCGGA", b"TCGTC", b"TCTAT", b"TCTGC", b"TGAAA", b"TGACT",
                b"TGAGA", b"TGAGG", b"TGATC", b"TGCAA", b"TGCAG", b"TGCAT", b"TGCCC", b"TGCGA",
                b"TGCGT", b"TGTTC", b"TTAAC", b"TTACA", b"TTACT", b"TTAGA", b"TTATG", b"TTCAT",
                b"TTCGG", b"TTGAC", b"TTGAT", b"TTGCG", b"TTTAC", b"TTTGA"
            ]
        );

        let counts: Vec<u64> = kmer_count.as_ref().values().cloned().collect();
        assert_eq!(counts, vec![0; counts.len()]);
        Ok(())
    }

    #[test]
    fn canonical_filter_path() -> anyhow::Result<()> {
        let mut kmers_in = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        kmers_in.write_all(FASTA_FILE)?;

        let mut reads_in = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        reads_in.write_all(FASTA_FILE)?;
        reads_in.write_all(FASTA_FILE)?;

        let output = tempfile::Builder::new().suffix(".fasta").tempfile()?;

        let kmer_count: KmerCount<atomic_counter::RelaxedCounter> =
            KmerCount::from_path(kmers_in.path(), 5, false)?;
        let _ = kmer_count.filter_path(
            Some(reads_in.path()),
            output.path(),
            0.0,
            100.0,
            false,
            false,
        );
        let mut counts: Vec<usize> = kmer_count.as_ref().values().map(|x| x.get()).collect();
        counts.sort();

        assert_eq!(
            counts,
            vec![
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                4, 6, 6, 6, 6, 8
            ]
        );

        Ok(())
    }

    const COUNT_OUTPUT: &[u8] = b"AAACT 2
AAAGT 2
AACAT 2
AACGC 2
AACGG 2
AACTT 4
AAGGC 2
AAGTC 2
AATAG 2
AATCT 4
AATGC 2
AATGG 2
AATTG 2
ACATG 4
ACCAT 2
ACCGT 2
ACCTA 2
ACCTG 2
ACGAT 2
ACGCA 2
ACGCG 2
ACGCT 2
ACGGG 2
ACTAT 2
ACTGA 2
ACTGC 2
ACTTA 2
ACTTG 2
AGAAT 2
AGACG 2
AGAGC 2
AGATA 2
AGCAG 2
AGCCT 2
AGCGA 2
AGCTA 2
AGCTG 2
AGGAT 2
AGGTC 2
AGGTG 2
AGTAA 4
AGTAG 2
AGTAT 2
AGTCA 2
AGTCC 2
AGTTA 2
ATACA 2
ATAGA 2
ATAGG 4
ATCAA 2
ATCGA 2
ATCGC 2
ATCTA 2
ATCTC 2
ATCTG 2
ATGAA 2
ATGCA 4
ATGCC 2
ATGGA 4
ATGTA 2
ATGTC 2
ATTAC 2
ATTAG 2
ATTCA 2
ATTGA 2
ATTGC 4
ATTGG 2
CAAAG 2
CAATA 2
CAATC 2
CAGAA 2
CAGAG 2
CAGCA 2
CAGCG 4
CAGTA 2
CAGTC 2
CATAA 2
CATAG 4
CATGC 4
CATTA 2
CATTC 2
CCATA 4
CCCTA 2
CCGAA 2
CCGTC 2
CCTCA 2
CCTGC 2
CGAAC 2
CGACG 2
CGATA 4
CGATC 2
CGCAA 2
CGCAG 2
CGGAC 2
CGGCA 2
CGTCA 2
CGTTA 2
CTATC 2
CTCAC 2
CTCAG 2
CTCTA 2
CTGAA 2
CTGCA 8
CTGCC 6
CTTAA 2
CTTGA 2
CTTTA 2
GAACA 2
GAATC 2
GACGA 2
GACGC 2
GAGAC 2
GAGCC 2
GATAC 4
GATCA 2
GCACC 4
GCAGA 6
GCAGC 6
GCGAC 2
GCGTA 2
GCTAC 2
GCTGA 2
GGACA 2
GGACC 2
GGATA 2
GGCAC 2
GGGCA 6
GGTGA 2
GTAAA 2
GTAAC 2
GTCAA 4
GTCCA 2
GTGCA 2
GTGTA 2
GTTAA 2
GTTTA 2
TACTA 2
TATGA 2
TCAAA 2
TCCGA 2
TCGCA 2
TCTAA 2
TCTCA 2
TGAAA 2
TGCAA 2
TGTAA 2
";

    #[test]
    fn to_csv() -> anyhow::Result<()> {
        let mut kmers_in = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        kmers_in.write_all(FASTA_FILE)?;

        let mut reads_in = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        reads_in.write_all(FASTA_FILE)?;
        reads_in.write_all(FASTA_FILE)?;

        let output = tempfile::Builder::new().suffix(".fasta").tempfile()?;

        let kmer_count: KmerCount<atomic_counter::RelaxedCounter> =
            KmerCount::from_path(kmers_in.path(), 5, false)?;
        let _ = kmer_count.filter_path(
            Some(reads_in.path()),
            output.path(),
            0.0,
            100.0,
            false,
            false,
        );

        let mut kmers_out = tempfile::Builder::new().suffix(".txt").tempfile()?;

        kmer_count.to_csv(kmers_out.path(), Some(b' '))?;
        let mut content = Vec::new();
        kmers_out.read_to_end(&mut content)?;

        let mut truth: Vec<Vec<u8>> = COUNT_OUTPUT
            .split(|x| *x == b'\n')
            .map(|x| x.to_vec())
            .collect();
        truth.sort();

        let mut result: Vec<Vec<u8>> = content.split(|x| *x == b'\n').map(|x| x.to_vec()).collect();
        result.sort();

        assert_eq!(result, truth);

        Ok(())
    }
}
