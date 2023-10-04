use clap::{ArgMatches};

use crate::get_km_paths;
use std::fs;
use std::process::Command as RunCommand;

///////////////////////// SUBCOMMAND INDEX KMERS /////////////////////////
/// This subcommand is used to index the kmers of a fasta file
/// It returns a tsv file with the following format:
/// kmers   D1
/// kmers:sequence5	0.009433962264150943

pub fn index_kmers(sub_matches: &ArgMatches) {
    let (kmtricks_path, kmindex_path) = get_km_paths(sub_matches);
    let inkmers = sub_matches.get_one::<String>("IN_KMERS").map(|s| s.clone()).unwrap();
    let outkmers = sub_matches.get_one::<String>("OUT_INDEX").map(|s| s.clone()).unwrap();
    
    let k = sub_matches.get_one::<u32>("K").map(|s| s.clone()).unwrap();
    let s = k - 3; // usees the findere approach 
    let t = sub_matches.get_one::<u32>("T").map(|s| s.clone()).unwrap();

    let bloom_size = sub_matches.get_one::<u64>("BLOOMSIZE").map(|s| s.clone()).unwrap();

    // check that inkmers is a non empty file:
    // Attempt to get metadata for the file
    if let Ok(metadata) = fs::metadata(inkmers.clone()) {
        // Check if the file exists
        if ! metadata.is_file() {
            panic!("{:#} exists, but it's not a file.", inkmers);
        }
    } else {
        panic!("The {} file does not exist or there was an error checking its existence.", inkmers);
    }

    // clear potential previous runs:
    let dir_name = "local_kmer_index";

    let mut cmd = RunCommand::new("rm");
    cmd.arg("-rf");
    cmd.arg(&outkmers);
    cmd.arg(&dir_name);
    let _ = cmd.output().expect("failed to execute process");

    // run kmindex
    let mut cmd = RunCommand::new(&kmindex_path);
    cmd.arg("build");
    cmd.arg("-f");
    cmd.arg(&inkmers);
    cmd.arg("--run-dir");
    cmd.arg(&dir_name);
    cmd.arg("--index");
    cmd.arg(&outkmers);
    cmd.arg("--register-as");
    cmd.arg("kmers");
    cmd.arg("-k");
    cmd.arg(&s.to_string());
    cmd.arg("-t");
    cmd.arg(&t.to_string());
    cmd.arg("--bloom-size");
    cmd.arg(&bloom_size.to_string());
    cmd.arg("--hard-min");
    cmd.arg("1");
    cmd.arg("--km-path");
    cmd.arg(&kmtricks_path);

    println!("Executing {:?}", cmd);
    let _ = cmd.output().expect("failed to execute process");
    // show the command stdout
    println!("kmindex done, the output directory is: {:?}", outkmers);
}