use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::Write;
use std::io::{self};
use fxread::initialize_reader;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <kmers.fasta> <reads.fasta> <output.fasta>", args[0]);
        return;
    }

    let kmer_file = &args[1];
    let reads_file = &args[2];
    let out_reads_file = &args[3];

    match index_kmers(kmer_file) {
        Ok((kmer_set, kmer_size)) => {
            let _ = read_fasta_file(reads_file, &kmer_set, kmer_size, out_reads_file);
            println!("Done, results are in file {}", out_reads_file);
        }
        Err(err) => eprintln!("Error indexing kmers: {}", err),
    }
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

fn canonical(kmer: &str) -> String {
    let rev_comp = reverse_complement(kmer);
    if kmer < &rev_comp {
        kmer.to_string()
    } else {
        rev_comp
    }
}

fn index_kmers(file_name: &str) -> Result<(HashSet<String>, usize), io::Error> {
    let mut kmer_set = HashSet::new();
    let mut kmer_size = 0;


    let reader = initialize_reader(&file_name).unwrap();
    for record in reader {
        // let header = record.id_str_checked().unwrap().to_string();
        let record_as_string = record.as_str_checked().unwrap().trim().as_bytes();
        let mut iter = record_as_string.split(|&x| x == b'\n');
        let _ = iter.next().unwrap();
        let acgt_sequence = iter.next().unwrap().to_owned();
        let string_acgt_sequence = String::from_utf8(acgt_sequence).expect("Found invalid UTF-8");
        if kmer_size == 0 {
            kmer_size = string_acgt_sequence.len();
        }
        kmer_set.insert(canonical(&&string_acgt_sequence.to_ascii_uppercase()));
    }
    println!("Indexed {} kmers, each of size {}", kmer_set.len(), kmer_size);
    // println!("Kmer size: {}", kmer_size);
    // // prints all kmers in the set
    // for kmers in &kmer_set {
    //     println!("{}", kmers);
    // }
    Ok((kmer_set, kmer_size))
}

fn read_fasta_file(file_name: &str, kmer_set: &HashSet<String>, kmer_size: usize, out_fasta: &str) -> std::io::Result<()>{
    
    let mut output = File::create(out_fasta)?;
    let reader = initialize_reader(&file_name).unwrap();
    for record in reader {
        let record_as_string = record.as_str_checked().unwrap().trim().as_bytes();
        let mut iter = record_as_string.split(|&x| x == b'\n');
        let stringheader = iter.next().unwrap();
        let acgt_sequence = iter.next().unwrap().to_owned();
        let string_acgt_sequence = String::from_utf8(acgt_sequence).expect("Found invalid UTF-8");
        let nb_shared_kmers = count_shared_kmers(kmer_set, &string_acgt_sequence, kmer_size);
        output.write_all(stringheader)?;
        output.write_all(b" ")?;
        output.write_all(nb_shared_kmers.to_string().as_bytes())?; 
        output.write_all(b"\n")?;
        output.write_all(string_acgt_sequence.as_bytes())?;
        output.write_all(b"\n")?;
        for line in iter {
            output.write_all(line)?;
            output.write_all(b"\n")?;
        }

        // println!("{} {}\n{}", header, count_shared_kmers(kmer_set, &record_as_string, kmer_size), record_as_string);
    }
    

    // let file = File::open(file_name)?;
    // let reader = BufReader::new(file);


    // for line in reader.lines() {
    //     let line = line?;   
    //     if !line.is_empty() && !line.starts_with('>') {
    //         // prints the number of shared kmers 
    //         println!("For read '{}', number of shared kmers: {}", line, count_shared_kmers(kmer_set, &line, kmer_size));
    //     }
    // }
    Ok(())
}

fn count_shared_kmers(kmer_set: &HashSet<String>, read: &str, kmer_size: usize) -> usize {
    let mut shared_kmers_count = 0;

    for i in 0..(read.len() - kmer_size + 1) {
        let kmer = &read[i..(i + kmer_size)];
        let canonical_kmer = canonical(&kmer.to_ascii_uppercase());
        if kmer_set.contains(&canonical_kmer){
            shared_kmers_count += 1;
        }
    }
    shared_kmers_count
}
