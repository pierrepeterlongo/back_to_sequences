use clap::ArgMatches;
use crate::{get_km_paths, Z};
use std::fs;
use std::process::Command as RunCommand;

///////////////////////// SUBCOMMAND GET HEADERS /////////////////////////
/// This subcommand is used to find headers of a fasta file that 
/// contain kmers that are indexed
/// 
pub fn get_headers (sub_matches: &ArgMatches) {
    let (_, kmindex_path) = get_km_paths(sub_matches);
    
    let index_kmers = sub_matches.get_one::<String>("INKMERS").map(|s| s.clone()).unwrap();
    let infasta = sub_matches.get_one::<String>("INFASTA").map(|s| s.clone()).unwrap();
    let outheaders = sub_matches.get_one::<String>("HEADERS").map(|s| s.clone()).unwrap();
    let t = sub_matches.get_one::<usize>("T").map(|s| s.clone()).unwrap();

    // check that index_kmers is a non empty directory:
    // Attempt to get metadata for the file
    if let Ok(metadata) = fs::metadata(index_kmers.clone()) {
        // Check if the file exists
        if ! metadata.is_dir() {
            panic!("{:#} exists, but it's not a directory.", index_kmers);
        }
    } else {
        panic!("The {} directory does not exist or there was an error checking its existence.", index_kmers);
    }

    // check that infasta is a non empty file:
    // Attempt to get metadata for the file
    if let Ok(metadata) = fs::metadata(infasta.clone()) {
        // Check if the file exists
        if ! metadata.is_file() {
            panic!("{:#} exists, but it's not a file.", infasta);
        }
    } else {
        panic!("The {} file does not exist or there was an error checking its existence.", infasta);
    }

    // clear potential previous runs:
    let mut cmd = RunCommand::new("rm");
    cmd.arg("-rf");
    cmd.arg(&outheaders);

    let printablecmd = format!("{:?}", cmd).replace("\" \"", " ").replace("\"", "");
    println!("Executing {:?}", printablecmd);
    let _ = cmd.output().expect("failed to execute process");

    let mut cmd = RunCommand::new(&kmindex_path);
    cmd.arg("query");
    cmd.arg("--index");
    cmd.arg(&index_kmers);
    cmd.arg("-q");
    cmd.arg(&infasta);
    cmd.arg("-o");
    cmd.arg(&outheaders);
    cmd.arg("-n");
    cmd.arg("kmers");
    cmd.arg("-z");
    cmd.arg(&Z.to_string()); // use the findere approach
    cmd.arg("-r");
    cmd.arg("0.0001");
    cmd.arg("--format");
    cmd.arg("matrix");
    cmd.arg("-t");
    cmd.arg(&t.to_string());

    let printablecmd = format!("{:?}", cmd).replace("\" \"", " ").replace("\"", "");
    println!("Executing {:?}", printablecmd);
    let _ = cmd.output().expect("failed to execute process");

    println!("kmindex done, the output directory is: {:?}", outheaders);


}