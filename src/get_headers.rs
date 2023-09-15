use clap::ArgMatches;
use crate::get_km_paths;
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
    let t = sub_matches.get_one::<u32>("T").map(|s| s.clone()).unwrap();

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
    let _ = cmd.output().expect("failed to execute process");

    // run kmindex
    // # create the kmindex command: 
    // kmindex_cmd="${current_dir_path}/bin/kmindex query --index ${indexed_kmers_file_path} -q ${queried_sequences_file_path} -o ${output_file_path} -n kmers -r 0.0001 --format matrix"

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
    cmd.arg("-r");
    cmd.arg("0.0001");
    cmd.arg("--format");
    cmd.arg("matrix");
    cmd.arg("-t");
    cmd.arg(&t.to_string());

    println!("kmindex command: {:?}", cmd);
    let _ = cmd.output().expect("failed to execute process");

    println!("kmindex done, the output directory is: {:?}", outheaders);


}