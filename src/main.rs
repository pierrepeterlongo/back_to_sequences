use clap::{arg, command, Command, value_parser, ArgMatches};
// extern crate clap;
 // also add dependency to Cargo.toml

use std::collections::HashMap;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use fasta::read::FastaReader;
use std::io::Write;

use std::fs;
use std::process::Command as RunCommand;
use which::which;

// use walkdir::WalkDir;
use glob::glob;


// get the paths of kmtricks and kmindex from arguments and environment variables
fn get_km_paths(sub_matches: &ArgMatches) -> (String, String) {
    let kmindex_path = sub_matches.get_one::<String>("KMINDEX_PATH"); // can be done in one line (see below)
    let kmindex_path = match kmindex_path {
        Some(kmindex_path) => {
            println!("kmindex is installed in {:?}",kmindex_path );
            kmindex_path.clone()
        },
        None => {
            let result = which("kmindex").expect("Path of kmindex undeclared and not found in the $PATH variable"); 
            println!("kmindex is installed in {:?}", result);
            result.to_str().unwrap().to_string()
        }
    };

    let kmtricks_path = match sub_matches.get_one::<String>("KMTRICKS_PATH") {
        Some(kmtricks_path) => {
            println!("kmtricks is installed in {:?}",kmtricks_path );
            kmtricks_path.clone()
        },
        None => {
            let result = which("kmtricks").expect("Path of kmtricks undeclared and not found in the $PATH variable"); 
            println!("kmtricks is installed in {:?}", result);
            result.to_str().unwrap().to_string()
        }
    };

    (kmtricks_path, kmindex_path)
}

///////////////////////// SUBCOMMAND INDEX KMERS /////////////////////////
/// This subcommand is used to index the kmers of a fasta file
/// It returns a tsv file with the following format:
/// kmers   D1
/// kmers:sequence5	0.009433962264150943

fn index_kmers(sub_matches: &ArgMatches) {
    let (kmtricks_path, kmindex_path) = get_km_paths(sub_matches);
    let inkmers = sub_matches.get_one::<String>("IN_KMERS").map(|s| s.clone()).unwrap();
    let outkmers = sub_matches.get_one::<String>("OUT_INDEX").map(|s| s.clone()).unwrap();
    
    let k = sub_matches.get_one::<u32>("K").map(|s| s.clone()).unwrap();
    let t = sub_matches.get_one::<u32>("T").map(|s| s.clone()).unwrap();

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
    cmd.arg(&k.to_string());
    cmd.arg("-t");
    cmd.arg(&t.to_string());
    cmd.arg("--bloom-size");
    cmd.arg("30000000");
    cmd.arg("--hard-min");
    cmd.arg("1");
    cmd.arg("--km-path");
    cmd.arg(&kmtricks_path);

    println!("Executing {:?}", cmd);
    let _ = cmd.output().expect("failed to execute process");
    // show the command stdout
    println!("kmindex done, the output directory is: {:?}", outkmers);
}

///////////////////////// SUBCOMMAND GET HEADERS /////////////////////////
/// This subcommand is used to find headers of a fasta file that 
/// contain kmers that are indexed
/// 
fn find_headers (sub_matches: &ArgMatches) {
    let (_, kmindex_path) = get_km_paths(sub_matches);
    
    let index_kmers = sub_matches.get_one::<String>("INKMERS").map(|s| s.clone()).unwrap();
    let infasta = sub_matches.get_one::<String>("INFASTA").map(|s| s.clone()).unwrap();
    let outheaders = sub_matches.get_one::<String>("OUTHEADERS").map(|s| s.clone()).unwrap();
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
    
///////////////////////// SUBCOMMAND TO_READS /////////////////////////

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn get_non_empty_headers (in_tsv_dir: String, threshold: f32) -> HashMap<String, f32> {
    // declare the map
    let mut map: HashMap<String, f32> = HashMap::new();
    let regex_pattern = format!("{}/**/*.tsv", in_tsv_dir);
    // read the tsv files in the in_tsv_dir directory
    for entry in glob(regex_pattern.as_str()).expect("Failed to read glob pattern") {     
        if let Ok(mut lines) = read_lines(entry.unwrap()) {
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

fn to_reads(sub_matches: &ArgMatches) {
    let intsvdir = sub_matches.get_one::<String>("INTSVDIR").map(|s| s.clone()).unwrap();
    let infasta = sub_matches.get_one::<String>("INFASTA").map(|s| s.clone()).unwrap();
    let outfasta = sub_matches.get_one::<String>("OUTFASTA").map(|s| s.clone()).unwrap();
    let threshold= sub_matches.get_one::<f32>("THRESHOLD").map(|s| s.clone()).unwrap();
    println!(
            "'toreads' was used, intsv is: {:?}, inasta is {:?} out is: {:?}, threshold is: {:?}",
            intsvdir, infasta, outfasta, threshold);
    let map = get_non_empty_headers(intsvdir, threshold);
    let _ = output_reads (map, infasta, outfasta);
}

///////////////////////// MAIN /////////////////////////

fn main() {
    let matches = command!() // requires `cargo` feature
        .propagate_version(true)
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("index_kmers")
            .arg(
                arg!([KMINDEX_PATH])
                    .value_parser(value_parser!(String))
                    .long("kmindex_path")
                    .help("Path to kmindex binary if not in the PATH")
            )
            .arg(
                arg!([KMTRICKS_PATH])
                    .value_parser(value_parser!(String))
                    .long("kmtricks_path")
                    .help("Path to kmtricks binary if not in the PATH")
            )
            .arg(
                arg!([IN_KMERS])
                    .value_parser(value_parser!(String))
                    .long("in_kmers")
                    .required(true)
                    .help("Input fasta file containing the kmers")
            )
            .arg(
                arg!([OUT_INDEX])
                    .value_parser(value_parser!(String))
                    .long("out_index")
                    .required(true)
                    .help("Indexed kmers, using kmindex")
            )
            .arg(
                arg!([K])
                    .value_parser(value_parser!(u32))
                    .short('k')
                    .default_value("31")
                    .help("kmer size")
            )
            .arg(
                arg!([T])
                    .value_parser(value_parser!(u32))
                    .short('t')
                    .default_value("8")
                    .help("kmindex number of threads")
            )
        )

        .subcommand(
            Command::new("get_headers")
            .arg(
                arg!([KMINDEX_PATH])
                    .value_parser(value_parser!(String))
                    .long("kmindex_path")
                    .help("Path to kmindex binary if not in the PATH")
            )
            .arg(
                arg!([KMTRICKS_PATH])
                    .value_parser(value_parser!(String))
                    .long("kmtricks_path")
                    .help("Path to kmtricks binary if not in the PATH")
            )
            .arg(
                arg!([INFASTA])
                    .value_parser(value_parser!(String))
                    .long("in_sequences")
                    .required(true)
                    .help("Input fasta file containing the original sequences in which we search the kmers")
            )
            .arg(
                arg!([INKMERS])
                    .value_parser(value_parser!(String))
                    .long("in_kmer_index")
                    .required(true)
                    .help("Directeroy name containing the indexed kmers, using kmindex (obtained with index_kmers --out_index)")
            )
            .arg(
                arg!([OUTHEADERS])
                    .value_parser(value_parser!(String))
                    .short('o').long("out_headers")
                    .required(true)
                    .help("The output directory that will contain the headers of the sequences that contain kmers in the indexed kmers")
            )
            .arg(
                arg!([T])
                    .value_parser(value_parser!(u32))
                    .short('t')
                    .default_value("8")
                    .help("kmindex number of threads")
            )
        )

        .subcommand(
            Command::new("to_reads")
            .arg(
                arg!([INTSVDIR])
                    .value_parser(value_parser!(String))
                    .long("in_tsv_dir")
                    .required(true)
                    .help("File containing the filtered headers, generated by kmindex")
            )
            .arg(
                arg!([INFASTA])
                    .value_parser(value_parser!(String))
                    .long("in_fasta")
                    .required(true)
                    .help("Input fasta file containing the original sequences")
            )
            .arg(
                arg!([OUTFASTA])
                    .value_parser(value_parser!(String))
                    .long("out_fasta")
                    .required(true)
                    .help("Output fasta file containing the retreived sequences")
            )
            .arg(arg!([THRESHOLD])
                .value_parser(value_parser!(f32))
                .default_value("0.0")
                .long("threshold"))
            .about("Do not retreive sequences whose ratio of shared kmers is below or equal to this value")
        )
        .get_matches();

    match matches.subcommand() {
        Some(("index_kmers", sub_matches)) => {
            index_kmers(&sub_matches)
        }
        ,
        Some(("to_reads", sub_matches)) => {
            to_reads(&sub_matches)
        }
        ,
        Some(("get_headers", sub_matches)) => {
            find_headers(&sub_matches)
        }
        _ => unreachable!("Exhausted list of subcommands and subcommand_required prevents `None`"),
    }

}
