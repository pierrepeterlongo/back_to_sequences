extern crate clap;
 // also add dependency to Cargo.toml

use std::collections::HashMap;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use fasta::read::FastaReader;
use std::io::Write;



use clap::Parser;

#[derive(Parser,Default,Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(short, long)]
    tsv_file: Option<String>,
    #[arg(short, long)]
    in_fasta: Option<String>,
    #[arg(short, long)]
    out_fasta: Option<String>,
    #[arg(long, default_value_t = 0.0)]
    threshold: f32,
}


// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn get_non_empty_headers (file_name: String, threshold: f32) -> HashMap<String, f32> {
    // declare the map
    let mut map: HashMap<String, f32> = HashMap::new();

    if let Ok(mut lines) = read_lines(file_name) {
        // Consumes the first line, must be 
        // kmers   D1
        lines.next();

        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(ip) = line {
                // a line is either as kmers:id        0
                // or as kmers:id        value bigger than 0
                // in the first case, we do not want to add it to the map
                // in the second case, we want to add it to the map with its value
                let mut iter = ip.split(":");
                let _kmers = iter.next().unwrap(); // don't need this

                let rest = iter.next().unwrap();
                let mut iter2 = rest.split("\t"); // should be as "sequence2	0"
                let sequence_id = String::from(iter2.next().unwrap());
                let ratio_kmers: f32 = iter2.next().unwrap().parse().unwrap();
                if ratio_kmers > threshold {
                    map.insert(sequence_id, ratio_kmers);
                }
            }
        }
    }
    map
}


fn output_reads (map: HashMap<String, f32>, in_fasta: String, out_fasta: String) -> std::io::Result<()> {
    
    // read the input fasta file
    // for each header, check if it is in the map
    // if it is, write the header and the sequence to the output fasta file
    // if it is not, do nothing
    // if the header is not in the map, do nothing
    // close the output fasta file
    // close the input fasta file
    
    let mut cnt = 0;
    let infile = Path::new(&in_fasta);  
    let mut output = File::create(out_fasta.clone())?;
    for [description, seq] in FastaReader::new(infile) {
        // if the header (description, removing the first ">") is in the map, write it to the output fasta file
        let mut header = description.clone();
        header.remove(0);
        if map.contains_key(&header) {
            cnt += 1;
            // write the two lines in the output file
            output.write_all(description.as_bytes())?;
            output.write_all(b" ")?;
            output.write_all(map.get(&header).unwrap().to_string().as_bytes())?;
            output.write_all(b"\n")?;
            output.write_all(seq.as_bytes())?;
            output.write_all(b"\n")?;
        }
    }
    println!("Number of sequences in the output fasta file {} : {}", out_fasta, cnt);
    Ok(())
}


fn main() {
    let args= Cli::parse();
    let map = get_non_empty_headers(args.tsv_file.unwrap(), args.threshold);
    let _ = output_reads (map, args.in_fasta.unwrap(), args.out_fasta.unwrap());
}
