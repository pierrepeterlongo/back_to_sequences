//! Define Command Line Interface

/* std use */

/* crates use */
use clap::Parser;

/* project use */

/// Extract sequences that contain some kmers
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Input fasta file containing the original kmers
    ///     Note: back_to_sequences considers the content as a set of kmers
    ///     This means that a kmer is considered only once,
    ///     even if it occurs multiple times in the file.
    ///     If the stranded option is not used (default), a kmer
    ///     and its reverse complement are considered as the same kmer.
    #[arg(long, verbatim_doc_comment)]
    pub in_kmers: String,

    /// Input fasta or fastq [.gz|zst] file containing the original sequences (eg. reads).
    ///     The stdin is used if not provided
    ///     (and if `--in_filelist` is not provided neither)
    #[arg(long, default_value_t = String::from(""), verbatim_doc_comment)]
    pub in_sequences: String,

    /// Input txt file containing in each line a path to a fasta or fastq [.gz|zst] file
    /// containing the original sequences (eg. reads).
    ///     Note1: if this option is used, the `--out_filelist` option must be used.
    ///            The number of lines in out_filelist must be the same as in_filelist
    ///     Note2: Incompatible with `--in_sequences`
    #[arg(long, default_value_t = String::from(""), verbatim_doc_comment)]
    pub in_filelist: String,

    /// Output file containing the filtered original sequences (eg. reads).
    /// It will be automatically in fasta or fastq format depending on the input file.
    /// If not provided, only the in_kmers with their count is output
    #[arg(long, default_value_t = String::from(""), verbatim_doc_comment)]
    pub out_sequences: String,

    /// Output txt file containing in each line a path to a fasta or fastq [.gz] file
    /// that will contain the related output file from the input files list
    #[arg(long, default_value_t = String::from(""), verbatim_doc_comment)]
    pub out_filelist: String,

    /// If provided, output a text file containing the kmers that occur in the reads
    /// with their
    ///  * number of occurrences
    ///     or
    ///  * their occurrence positions if the --output_kmer_positions option is used
    ///     Note: if `--in_filelist` is used the output counted kmers are
    ///     those occurring the last input file of that list
    #[arg(long, default_value_t = String::from(""), verbatim_doc_comment)]
    pub out_kmers: String,

    /// If out_kmers is provided, output only reference kmers whose number of occurrences
    /// is at least equal to this value.
    /// If out_kmers is not provided, this option is ignored
    #[arg(long, default_value_t = 0, verbatim_doc_comment)]
    pub counted_kmer_threshold: usize,

    /// If out_kmers is provided, either only count their number of occurrences (default)
    /// or output their occurrence positions (read_id, position, strand) if this option is used
    #[arg(long, default_value_t = false, verbatim_doc_comment)]
    pub output_kmer_positions: bool,

    /// If provided, output matching positions on sequences in the
    /// out_sequence file(s)
    #[arg(long, default_value_t = false, verbatim_doc_comment)]
    pub output_mapping_positions: bool,

    /// Size of the kmers to index and search
    #[arg(short, long, default_value_t = 31)]
    pub kmer_size: usize,

    /// Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
    /// Minimal threshold of the ratio  (%) of kmers that must be found in a sequence to keep it (default 0%).
    /// Thus by default, if no kmer is found in a sequence, it is not output.
    #[arg(short, long, default_value_t = 0.0, verbatim_doc_comment)]
    pub min_threshold: f32,

    /// Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
    /// Maximal threshold of the ratio (%) of kmers that must be found in a sequence to keep it (default 100%).
    /// Thus by default, there is no limitation on the maximal number of kmers found in a sequence.
    #[arg(long, default_value_t = 100.0, verbatim_doc_comment)]
    pub max_threshold: f32,

    /// Used original kmer strand (else canonical kmers are considered)
    #[arg(long, default_value_t = false)]
    pub stranded: bool,

    /// Query the reverse complement of reads. Useless without the --stranded option
    #[arg(long, default_value_t = false)]
    pub query_reverse: bool,

    /// Do not index low complexity kmers (ie. with a Shannon entropy < 1.0)
    #[arg(long, default_value_t = false)]
    pub no_low_complexity: bool,

    /// Number of threads
    ///    Note: if not provided, the number of threads is set to the number of logical cores
    #[arg(short, long, default_value_t = 0, verbatim_doc_comment)]
    pub threads: usize,
}

/// check that a file name corresponds to a non empty file:
pub fn validate_non_empty_file(in_file: String) -> anyhow::Result<()> {
    if let Ok(metadata) = std::fs::metadata(in_file.clone()) {
        // Check if the file exists
        if !metadata.is_file() {
            anyhow::bail!("{:#} exists, but it's not a file.", in_file)
        }
    } else {
        anyhow::bail!("{:#} is not a file.", in_file)
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_empty_file_test() -> anyhow::Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let directory = temp_dir.into_path();
        let file = directory.join("empty.fasta");

        // test directory
        assert!(
            validate_non_empty_file(directory.clone().into_os_string().into_string().unwrap())
                .is_err()
        );

        // test not exist
        assert!(
            validate_non_empty_file(file.clone().into_os_string().into_string().unwrap()).is_err()
        );

        // test work
        std::fs::File::create(&file)?;
        assert!(
            validate_non_empty_file(file.clone().into_os_string().into_string().unwrap()).is_ok()
        );

        Ok(())
    }
}
