use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::{self};
use fxread::initialize_reader;

fn validate_non_empty_file(in_file: String) {
    // check that inkmers is a non empty file:
    // Attempt to get metadata for the file
    if let Ok(metadata) = fs::metadata(in_file.clone()) {
        // Check if the file exists
        if ! metadata.is_file() {
            panic!("{:#} exists, but it's not a file.", in_file);
        }
    } else {
        panic!("The {} file does not exist or there was an error checking its existence.", in_file);
    }
}

pub fn validate_kmers(in_fasta_reads: String, in_fasta_kmers: String, out_fasta_reads:String, out_txt_kmers: String, kmer_size: usize, stranded: bool) -> std::io::Result<()> {
  
    // check that inkmers and reads_file are non empty files:
    validate_non_empty_file(in_fasta_reads.clone());
    validate_non_empty_file(in_fasta_kmers.clone());


    match index_kmers(in_fasta_kmers, kmer_size, stranded) {
        Ok((mut kmer_set, kmer_size)) => {
            let _ = count_kmers_in_fasta_file(in_fasta_reads, &mut kmer_set, kmer_size, out_fasta_reads.clone(), stranded);
            println!("Filtered reads with exact kmer count are in file {}", out_fasta_reads);
            

            // if the out_kmers_file is not empty, we output counted kmers in the out_kmers_file file
            if out_txt_kmers.len() > 0 {
                
                // prints all kmers from kmer_set that have a count > 0
                let mut output = File::create(out_txt_kmers.clone())?;
                // let mut output = File::create(out_kmers_file);
                for (kmer, count) in kmer_set.iter() {
                    if *count > 0 {
                        output.write_all(kmer.as_bytes())?;
                        output.write_all(b" ")?;
                        output.write_all(count.to_string().as_bytes())?;
                        output.write_all(b"\n")?;
                    }
                }
            println!("kmers with their number of occurrences in the original reads are in file {}", out_txt_kmers);
            }
        }

        Err(err) => eprintln!("Error indexing kmers: {}", err),
    }
    Ok(())

}

fn reverse_complement(kmer: &str) -> String {
    kmer.chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => base,
        })
        .collect()
}

fn canonical(kmer: &str, stranded: bool) -> String {
    if stranded  {return kmer.to_string();}
    let rev_comp = reverse_complement(kmer);
    if kmer < &rev_comp {
        kmer.to_string()
    } else {
        rev_comp
    }
}

fn index_kmers(file_name: String, kmer_size: usize, stranded: bool) -> Result<(HashMap<String, i32>, usize), io::Error> {
    let mut kmer_set: HashMap<String, i32> = HashMap::new();

    let reader = initialize_reader(&file_name).unwrap();
    for record in reader {
        let record_as_string = record.as_str_checked().unwrap().trim().as_bytes();
        let mut iter = record_as_string.split(|&x| x == b'\n');
        let _ = iter.next().unwrap();
        let acgt_sequence = iter.next().unwrap().to_owned();
        let string_acgt_sequence = String::from_utf8(acgt_sequence).expect("Found invalid UTF-8");
        // for each kmer of the sequence, insert it in the kmer_set
        for i in 0..(string_acgt_sequence.len() - kmer_size + 1) {
            let kmer = &string_acgt_sequence[i..(i + kmer_size)];
            kmer_set.insert(canonical(&&kmer.to_ascii_uppercase(), stranded), 0);
        }
        // kmer_set.insert(canonical(&&string_acgt_sequence.to_ascii_uppercase()), 0);
    }
    println!("Indexed {} kmers, each of size {}", kmer_set.len(), kmer_size);
    
    Ok((kmer_set, kmer_size))
}

fn count_kmers_in_fasta_file(file_name: String, kmer_set:  &mut HashMap<String, i32>, kmer_size: usize, out_fasta: String, stranded: bool) -> std::io::Result<()>{
    
    let mut output = File::create(out_fasta)?;
    let reader = initialize_reader(&file_name).unwrap();
    for record in reader {
        let record_as_string = record.as_str_checked().unwrap().trim().as_bytes();
        let mut iter = record_as_string.split(|&x| x == b'\n');
        let stringheader = iter.next().unwrap();
        let acgt_sequence = iter.next().unwrap().to_owned();
        let string_acgt_sequence = String::from_utf8(acgt_sequence).expect("Found invalid UTF-8");
        let nb_shared_kmers = count_shared_kmers(kmer_set, &string_acgt_sequence, kmer_size, stranded);
        let ratio_kmers = nb_shared_kmers as f32 / (string_acgt_sequence.len() - kmer_size + 1) as f32;
        if 
        output.write_all(stringheader)?;
        output.write_all(b" ")?;
        output.write_all(nb_shared_kmers.to_string().as_bytes())?; 
        output.write_all(b" ")?;
        output.write_all(ratio_kmers.to_string().as_bytes())?;
        output.write_all(b"\n")?;
        output.write_all(string_acgt_sequence.as_bytes())?;
        output.write_all(b"\n")?;
        for line in iter {
            output.write_all(line)?;
            output.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn count_shared_kmers(kmer_set: &mut HashMap<String, i32>, read: &str, kmer_size: usize, stranded: bool) -> usize {
    let mut shared_kmers_count = 0;

    for i in 0..(read.len() - kmer_size + 1) {
        let kmer = &read[i..(i + kmer_size)];
        let canonical_kmer = canonical(&kmer.to_ascii_uppercase(), stranded);
        if kmer_set.contains_key(&canonical_kmer){
            shared_kmers_count += 1;
            *kmer_set.get_mut(&canonical_kmer).unwrap() += 1;
        }
    }
    shared_kmers_count
}
