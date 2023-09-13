
 // also add dependency to Cargo.toml

use std::collections::HashMap;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;


// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn parse_input_file (file_name: String) -> HashMap<String, f32> {
    // declare the map
    let mut map: HashMap<String, f32> = HashMap::new();

    if let Ok(mut lines) = read_lines(file_name) {
        // Consumes the line, must be 
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
                if ratio_kmers > 0.0 {
                    map.insert(sequence_id, ratio_kmers);
                }

            }
        }
    }
    map
}


fn main() {
    let map = parse_input_file(String::from("/Users/ppeterlo/workspace/kmer2sequences/output/kmers.tsv"));
    for (key, val) in map.iter() {
        println!("key: {key} val: {val}");
    }
}
