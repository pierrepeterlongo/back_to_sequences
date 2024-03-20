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

/// index all kmers of size kmer_size in the fasta file
/// returns a hashmap with the kmers as keys and their count as values, initialized to 0
pub fn index_kmers<T: Default>(
    file_name: String,
    kmer_size: usize,
    stranded: bool,
    no_low_complexity: bool,
) -> anyhow::Result<(HashMap<Vec<u8>, T>, usize)> {
// ) -> anyhow::Result<(Box<dyn HashMap<Vec<u8>, T>>, usize)> {
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
            } else {
                // If the entropy is too low, the kmer is not inserted
                if no_low_complexity && shannon_entropy(kmer) < 1.0 {
                    i = i + 1;
                    continue;
                }
                kmer_set.insert(
                    SequenceNormalizer::new(kmer, reverse_complement)
                        .iter()
                        .collect(),
                    Default::default(), // RelaxedCounter::new(0) // TODO call default from kmer_counter (anthony)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_acgt() {
        let mut rng = crate::tests::rand();

        let sequence = crate::tests::sequence(&mut rng, 50).to_ascii_uppercase();
        assert_eq!((true, 0), first_non_acgt(&sequence));

        let mut sequence = crate::tests::sequence(&mut rng, 10).to_ascii_uppercase();
        sequence.push(b'@');
        assert_eq!((false, 10), first_non_acgt(&sequence));

        let sequence = crate::tests::sequence(&mut rng, 0).to_ascii_uppercase();
        assert_eq!((true, 0), first_non_acgt(&sequence));
    }

    #[test]
    fn build_index_kmers() -> anyhow::Result<()> {
        let mut rng = crate::tests::rand();

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fasta");

        crate::tests::io::write_fasta(&mut rng, 16, 5, 5, &kmers_in_path)?;

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
                b"AAGCATTACCGTGGC".to_vec(),
                b"ACTGTGCAAAAGGGG".to_vec(),
                b"AGACATGAGCAACCA".to_vec(),
                b"AGGATATCGAATTAT".to_vec(),
                b"ATGAATCGCGTGTTA".to_vec(),
                b"CAAGCATTACCGTGG".to_vec(),
                b"CAGACATGAGCAACC".to_vec(),
                b"CAGGATATCGAATTA".to_vec(),
                b"CCCTTTTGCACAGTA".to_vec(),
                b"CTAACACGCGATTCA".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();
        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }

    // anthony : pq "failed to resolve: unresolved import"
    #[test]
    fn build_index_kmers_stranded() -> anyhow::Result<()> {
        let mut rng = crate::tests::rand();

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fastq");

        crate::tests::io::write_fastq(&mut rng, 16, 5, 5, &kmers_in_path)?;

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
                b"ACATGCTGCAATTAC".to_vec(),
                b"ATCCTCTGGAACTTG".to_vec(),
                b"ATGAATCGCGTGTTA".to_vec(),
                b"CATCCTCTGGAACTT".to_vec(),
                b"GACATGCTGCAATTA".to_vec(),
                b"TGAATCGCGTGTTAG".to_vec(),
                b"TGCTCATGTCTGCTG".to_vec(),
                b"TGTACGCAGGATATC".to_vec(),
                b"TTGCTCATGTCTGCT".to_vec(),
                b"TTGTACGCAGGATAT".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();
        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }

    #[test]
    fn build_index_kmers_no_low_complexity() -> anyhow::Result<()> {
        let mut rng = crate::tests::rand();

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fasta");

        let mut fasta_data = crate::tests::fasta::records(&mut rng, 5, 16, 2);
        fasta_data.extend(b">read_homopolymer\n");
        fasta_data.extend([b'A'; 16]);
        fasta_data.push(b'\n');
        fasta_data.extend(crate::tests::fasta::records(&mut rng, 5, 16, 2));

        crate::tests::io::write_buffer(&fasta_data, &kmers_in_path)?;

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
                b"AAGCATTACCGTGGC".to_vec(),
                b"AGACATGAGCAACCA".to_vec(),
                b"AGGATATCGAATTAT".to_vec(),
                b"ATGAATCGCGTGTTA".to_vec(),
                b"CAAGCATTACCGTGG".to_vec(),
                b"CAGACATGAGCAACC".to_vec(),
                b"CAGGATATCGAATTA".to_vec(),
                b"CTAACACGCGATTCA".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();
        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }

    #[test]
    fn not_nucleotide() -> anyhow::Result<()> {
        let mut rng = crate::tests::rand();

        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();

        let kmers_in_path = temp_path.join("kmers_in.fasta");

        let mut fasta_data = crate::tests::fasta::records(&mut rng, 5, 16, 2);
        fasta_data.extend(crate::tests::fasta::record(&mut rng, 10, 16));
        *fasta_data.last_mut().unwrap() = b'@';
        fasta_data.extend(b"ACTG\n");
        fasta_data.extend(crate::tests::fasta::records(&mut rng, 5, 16, 2));

        crate::tests::io::write_buffer(&fasta_data, &kmers_in_path)?;

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
                b"AAGCATTACCGTGGC".to_vec(),
                b"AGCAGACATGAGCAA".to_vec(),
                b"ATGAATCGCGTGTTA".to_vec(),
                b"CAAGCATTACCGTGG".to_vec(),
                b"CAGCAGACATGAGCA".to_vec(),
                b"CTAACACGCGATTCA".to_vec(),
                b"CTATAATTCGATATC".to_vec(),
                b"CTCCCCTTTTGCACA".to_vec(),
                b"CTGTGCAAAAGGGGA".to_vec(),
                b"GGATATCGAATTATA".to_vec(),
            ]
        );

        let mut values = index.values().map(|x| x.get()).collect::<Vec<usize>>();
        values.sort_unstable();
        assert_eq!(values, vec![0usize; keys.len()]);

        Ok(())
    }
}
