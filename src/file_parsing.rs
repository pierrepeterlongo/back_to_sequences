//! Back to sequences: find the origin of kmers

/* std use */
use std::fs::File;
use std::io::{BufRead, BufReader, Error};
use std::path::Path;

/// Parses a file and returns a vector of Strings
/// each line in the file is a String
/// This is used to parse input / output file lists
pub fn read_file_lines(file_path: &str) -> Result<Vec<String>, Error> {
    let path = Path::new(file_path);
    let file = File::open(path)?;

    let reader = BufReader::new(file);
    let mut lines = Vec::new();

    for line in reader.lines() {
        lines.push(line?);
    }

    Ok(lines)
}

#[cfg(test)]
mod tests {
    /* project use */
    use super::*;

    use std::io::Write as _;

    const INPUTS: &[u8] = b"input_1.fasta
input_2.fasta
input_3.fasta
";

    #[test]
    fn inputs_file() -> std::io::Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let temp_path = temp_dir.path();
        let path = temp_path.join("input.lst");

        std::fs::File::create(format!("{}", path.display()))?.write_all(INPUTS)?;

        let inputs = read_file_lines(&format!("{}", path.display()))?;

        assert_eq!(
            inputs,
            vec![
                "input_1.fasta".to_string(),
                "input_2.fasta".to_string(),
                "input_3.fasta".to_string(),
            ]
        );

        Ok(())
    }
}
