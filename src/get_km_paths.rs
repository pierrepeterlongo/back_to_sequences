use which::which;
use clap::{ArgMatches};

// get the paths of kmtricks and kmindex from arguments and environment variables
pub fn get_km_paths(sub_matches: &ArgMatches) -> (String, String) {
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