//! Kmer hash declarations

/* std use */

/* crates use */
use ahash::AHashMap as HashMap;
use entropy::shannon_entropy;
use needletail::{Sequence, parse_fastx_file};

/* project use */
use crate::{kmer_counter::KmerCounter, sequence_normalizer::SequenceNormalizer};

/// given a kmer as a &[u8] return a tuple boolean, position
/// if the kmer contains an non ACGT letter, return false and the position of the first non ACGT letter
/// As we used needletail to normalize the sequences, we only check for 'N'
/// else return true and 0 as position
pub fn first_non_acgt(kmer: &[u8]) -> (bool, usize) {
    for (i, &byte) in kmer.iter().enumerate() {
        if byte == b'N' {
            return (false, i);
        }
        // if byte != b'A' && byte != b'C' && byte != b'G' && byte != b'T' {
            // return (false, i);
        // }
    }
    (true, 0)
}

/// index all kmers of size kmer_size in the fasta file
/// returns a hashmap with the kmers as keys and their count as values, initialized to 0
pub fn index_kmers<T: KmerCounter>(
    file_name: String,
    kmer_size: usize,
    stranded: bool,
    no_low_complexity: bool,
) -> anyhow::Result<(HashMap<Vec<u8>, T>, usize)> {
    // ) -> anyhow::Result<(Box<dyn HashMap<Vec<u8>, T>>, usize)> {
    let mut kmer_set = HashMap::new();
    let reverse_complement = if stranded { Some(false) } else { None };

    let mut reader = parse_fastx_file(&file_name).unwrap();
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        // transform norm_seq into a string on the ACGTN alphabet
        let acgt_sequence = norm_seq.iter().cloned().collect::<Vec<u8>>(); // TODO optimize

    // loop {
    //     let Some(mut record) = reader.next()? else {
    //         break;
    //     };
    //     record.upper();
    //     let acgt_sequence = record.seq();

        // for each kmer of the sequence, insert it in the kmer_set
        if acgt_sequence.len() < kmer_size {
            continue;
        }
        let mut i = 0;
        while i < acgt_sequence.len() - kmer_size + 1 {
            // for mut i in 0..(acgt_sequence.len() - kmer_size + 1) {
            let kmer = &acgt_sequence[i..(i + kmer_size)];
            let first_non_acgt = first_non_acgt(kmer);
            if !first_non_acgt.0 {
                // If the kmer contains a non acgt letter, we jump to the next possible kmer
                i = i + first_non_acgt.1 + 1;
                continue;
            } else {
                // If the entropy is too low, the kmer is not inserted
                if no_low_complexity && shannon_entropy(kmer) < 1.0 {
                    i += 1;

                    continue;
                }
                kmer_set.insert(
                    SequenceNormalizer::new(kmer, reverse_complement)
                        .iter()
                        .collect(),
                    Default::default(), // RelaxedCounter::new(0) // TODO call default from kmer_counter (anthony)
                );
            }
            i += 1;
        }
    }
    eprintln!(
        "Indexed {} kmers, each of size {}",
        kmer_set.len(),
        kmer_size
    );

    Ok((kmer_set, kmer_size))
}

#[cfg(test)]
mod tests {
    /* std use */
    use std::io::Write as _;

    /* crate use */
    use atomic_counter::AtomicCounter as _;
    use biotest::values::Generate as _;
    use biotest::Format as _;

    /* project use */
    use super::*;

    #[test]
    fn non_acgt() -> anyhow::Result<()> {
        let mut rng = biotest::rand();

        let sequence = biotest::values::Nucleotides::DnaUpper.generate(&mut rng, 50)?;
        println!("{}", String::from_utf8(sequence.clone()).unwrap());
        assert_eq!((true, 0), first_non_acgt(&sequence));

        let mut sequence = biotest::values::Nucleotides::DnaUpper.generate(&mut rng, 10)?;
        sequence.push(b'@');
        assert_eq!((false, 10), first_non_acgt(&sequence));

        let sequence = biotest::values::Nucleotides::DnaUpper.generate(&mut rng, 0)?;
        assert_eq!((true, 0), first_non_acgt(&sequence));

        Ok(())
    }

    #[test]
    fn build_index_kmers() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let k_generate = biotest::Fasta::builder().sequence_len(16).build()?;

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fasta");

        k_generate.create(&kmers_in_path, &mut rng, 5)?;

        let (index, kmer_size) = index_kmers::<atomic_counter::RelaxedCounter>(
            kmers_in_path.into_os_string().into_string().unwrap(),
            15,
            false,
            false,
        )?;

        assert_eq!(kmer_size, 15);

        let mut keys = index.keys().cloned().collect::<Vec<Vec<u8>>>();
        keys.sort_unstable();
        assert_eq!(
            keys,
            vec![
                b"AACTTGCGACAAGAA".to_vec(),
                b"AAGCATTACCGTGGC".to_vec(),
                b"ACGTTAAGAAGGTTC".to_vec(),
                b"AGCAACCATCTATAA".to_vec(),
                b"ATTGCAGCATGTCTC".to_vec(),
                b"CAAGCATTACCGTGG".to_vec(),
                b"CGAACCTTCTTAACG".to_vec(),
                b"GAACTTGCGACAAGA".to_vec(),
                b"GAGCAACCATCTATA".to_vec(),
                b"GGAGACATGCTGCAA".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();

        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }

    #[test]
    fn build_index_kmers_stranded() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let k_generate = biotest::Fasta::builder().sequence_len(16).build()?;

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fasta");

        k_generate.create(&kmers_in_path, &mut rng, 5)?;

        let (index, kmer_size) = index_kmers::<atomic_counter::RelaxedCounter>(
            kmers_in_path.into_os_string().into_string().unwrap(),
            15,
            true,
            false,
        )?;

        assert_eq!(kmer_size, 15);

        let mut keys = index.keys().cloned().collect::<Vec<Vec<u8>>>();
        keys.sort_unstable();
        assert_eq!(
            keys,
            vec![
                b"AACTTGCGACAAGAA".to_vec(),
                b"CCACGGTAATGCTTG".to_vec(),
                b"CGAACCTTCTTAACG".to_vec(),
                b"GAACCTTCTTAACGT".to_vec(),
                b"GAACTTGCGACAAGA".to_vec(),
                b"GAGACATGCTGCAAT".to_vec(),
                b"GCCACGGTAATGCTT".to_vec(),
                b"GGAGACATGCTGCAA".to_vec(),
                b"TATAGATGGTTGCTC".to_vec(),
                b"TTATAGATGGTTGCT".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();

        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }

    #[test]
    fn build_index_kmers_no_low_complexity() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let k_generate = biotest::Fasta::builder().sequence_len(16).build()?;

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fasta");
        let mut fasta_data = vec![];

        k_generate.records(&mut fasta_data, &mut rng, 5)?;
        fasta_data.extend(b">read_homopolymer\n");
        fasta_data.extend([b'A'; 16]);
        fasta_data.push(b'\n');

        std::fs::File::create(&kmers_in_path)?.write_all(&fasta_data)?;

        let (index, kmer_size) = index_kmers::<atomic_counter::RelaxedCounter>(
            kmers_in_path.into_os_string().into_string().unwrap(),
            15,
            false,
            true,
        )?;

        assert_eq!(kmer_size, 15);

        let mut keys = index.keys().cloned().collect::<Vec<Vec<u8>>>();
        keys.sort_unstable();
        assert_eq!(
            keys,
            vec![
                b"AACTTGCGACAAGAA".to_vec(),
                b"AAGCATTACCGTGGC".to_vec(),
                b"ACGTTAAGAAGGTTC".to_vec(),
                b"AGCAACCATCTATAA".to_vec(),
                b"ATTGCAGCATGTCTC".to_vec(),
                b"CAAGCATTACCGTGG".to_vec(),
                b"CGAACCTTCTTAACG".to_vec(),
                b"GAACTTGCGACAAGA".to_vec(),
                b"GAGCAACCATCTATA".to_vec(),
                b"GGAGACATGCTGCAA".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();
        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }

    #[test]
    fn not_nucleotide() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let k_generate = biotest::Fasta::builder().sequence_len(16).build()?;

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fasta");
        let mut fasta_data = vec![];

        k_generate.records(&mut fasta_data, &mut rng, 5)?;
        k_generate.record(&mut fasta_data, &mut rng)?;
        fasta_data.extend(b"@ACTG\n");
        k_generate.records(&mut fasta_data, &mut rng, 5)?;

        std::fs::File::create(&kmers_in_path)?.write_all(&fasta_data)?;

        let (index, kmer_size) = index_kmers::<atomic_counter::RelaxedCounter>(
            kmers_in_path.into_os_string().into_string().unwrap(),
            15,
            false,
            false,
        )?;

        assert_eq!(kmer_size, 15);

        let mut keys = index.keys().cloned().collect::<Vec<Vec<u8>>>();
        keys.sort_unstable();
        assert_eq!(
            keys,
            vec![
                b"AACTTGCGACAAGAA".to_vec(),
                b"AAGCATTACCGTGGC".to_vec(),
                b"ACGTTAAGAAGGTTC".to_vec(),
                b"AGATTTGTGCTTAAG".to_vec(),
                b"AGCAACCATCTATAA".to_vec(),
                b"ATGTCAGGCTAGTTC".to_vec(),
                b"ATTGCAGCATGTCTC".to_vec(),
                b"ATTTAAAACGGCTTG".to_vec(),
                b"CAAGCATTACCGTGG".to_vec(),
                b"CATTGCGCTTGAAGG".to_vec(),
                b"CCAAGCCGTTTTAAA".to_vec(),
                b"CCTATGCTCACTCAA".to_vec(),
                b"CCTTAAGCACAAATC".to_vec(),
                b"CGAACCTTCTTAACG".to_vec(),
                b"CTATTGTAGGCGCAC".to_vec(),
                b"CTTCAAGCGCAATGA".to_vec(),
                b"GAACTTGCGACAAGA".to_vec(),
                b"GAGCAACCATCTATA".to_vec(),
                b"GGAACTAGCCTGACA".to_vec(),
                b"GGAGACATGCTGCAA".to_vec(),
                b"TCCTATGCTCACTCA".to_vec(),
                b"TCTATTGTAGGCGCA".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();

        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }
}
