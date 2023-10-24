use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::{self};
use fxread::initialize_reader;
use atomic_counter::{RelaxedCounter, AtomicCounter};
use rayon::prelude::*;



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
fn index_kmers<T:Default>(file_name: String, kmer_size: usize, stranded: bool) -> Result<(HashMap<String, T>, usize), io::Error> {
    let mut kmer_set = HashMap::new();

    let reader = initialize_reader(&file_name).unwrap();
    for record in reader {
        let acgt_sequence = record.seq_str_checked().expect("Found invalid UTF-8");
        // for each kmer of the sequence, insert it in the kmer_set
        for i in 0..(acgt_sequence.len() - kmer_size + 1) {
            let kmer = &acgt_sequence[i..(i + kmer_size)];
        }
        // kmer_set.insert(canonical(&&string_acgt_sequence.to_ascii_uppercase()), 0);
    }
    println!("Indexed {} kmers, each of size {}", kmer_set.len(), kmer_size);
    
    Ok((kmer_set, kmer_size))
}

fn round(x: f32, decimals: u32) -> f32 {
    let y = 10i32.pow(decimals) as f32;
    (x * y).round() / y
}

fn count_kmers_in_fasta_file_par(file_name: String, 
    kmer_set:  &HashMap<String, atomic_counter::RelaxedCounter>, 
    kmer_size: usize, 
    out_fasta: String, 
    threshold: f32, 
    stranded: bool, 
    query_reverse: bool) -> std::io::Result<()>{
    let output = File::create(out_fasta)?;
    let write_lock = std::sync::Mutex::new(output);
    let (tx, rx) = std::sync::mpsc::sync_channel(1024);
    let (_, result) = rayon::join(move ||{// lance deux threads 
        let reader = initialize_reader(&file_name).unwrap();
        for record in reader {
            tx.send(record).unwrap();
        }
    }, ||{
        rx.into_iter().par_bridge().try_for_each(|mut record| -> std::io::Result<()>{
            if query_reverse {
                record.rev_comp();  // reverse the sequence in place
            }
            let record_as_string = record.as_str().trim().as_bytes(); 
            let mut iter = record_as_string.split(|&x| x == b'\n');
            let stringheader = iter.next().unwrap();
            let acgt_sequence = record.seq_str_checked().expect("Found invalid UTF-8");
            let nb_shared_kmers = count_shared_kmers_par(kmer_set, acgt_sequence, kmer_size, stranded);
            let ratio_shared_kmers = nb_shared_kmers as f32 / (acgt_sequence.len() - kmer_size + 1) as f32;
            
            if ratio_shared_kmers > threshold{ // supports the user defined threshold
                // round ratio_shared_kmers to 3 decimals and transform to percents
                let percent_shared_kmers = round(ratio_shared_kmers*100.0, 2);
                let mut out = write_lock.lock().unwrap();
                out.write_all(stringheader)?;
                write!(out, " {} {}\n", nb_shared_kmers, percent_shared_kmers)?;
                for line in iter {
                    out.write_all(line)?;
                    out.write_all(b"\n")?;
                }
            } // end read contains at least one indexed kmer
        Ok(())
        }) // end of for each
    }); // end of rayon join
    result
}

fn count_shared_kmers_par(kmer_set:  &HashMap<String, atomic_counter::RelaxedCounter>, read: &str, kmer_size: usize, stranded: bool) -> usize {
    let mut shared_kmers_count = 0;

    for i in 0..(read.len() - kmer_size + 1) {
        let kmer = &read[i..(i + kmer_size)];
        let canonical_kmer = canonical(&kmer.to_ascii_uppercase(), stranded);
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



pub fn validate_kmers(in_fasta_reads: String, 
    in_fasta_kmers: String, 
    out_fasta_reads:String, 
    out_txt_kmers: String, 
    kmer_size: usize, 
    threshold: f32, 
    stranded: bool,
    query_reverse: bool) -> std::io::Result<()> {
      
    // check that inkmers and reads_file are non empty files:
    validate_non_empty_file(in_fasta_reads.clone());
    validate_non_empty_file(in_fasta_kmers.clone());

    
    match index_kmers::<RelaxedCounter>(in_fasta_kmers, kmer_size, stranded) {

        Ok((kmer_set, kmer_size)) => {
            let _ = count_kmers_in_fasta_file_par(in_fasta_reads, &kmer_set, kmer_size, out_fasta_reads.clone(), threshold, stranded, query_reverse);
            println!("Filtered sequences with exact kmer count are in file {}", out_fasta_reads);
            

            // if the out_kmers_file is not empty, we output counted kmers in the out_kmers_file file
            if out_txt_kmers.len() > 0 {
                
                // prints all kmers from kmer_set that have a count > 0
                let mut output = File::create(out_txt_kmers.clone())?;
                // let mut output = File::create(out_kmers_file);
                for (kmer, count) in kmer_set.iter() {
                    if count.get() > 0 {
                        write!(output, "{} {}\n", kmer, count.get())?;
                    }
                }
            println!("kmers with their number of occurrences in the original sequences are in file {}", out_txt_kmers);
            }
        }

        Err(err) => eprintln!("Error indexing kmers: {}", err),
    }
    Ok(())

}
