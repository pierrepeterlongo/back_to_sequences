//! Functional test

/* std use */

/* crate use */

/* project use */

use std::io::{Read, Write};

use biotest::Format;

#[test]
fn help() -> std::result::Result<(), anyhow::Error> {
    let mut cmd = assert_cmd::Command::cargo_bin("back_to_sequences")?;

    cmd.args(["--help"]);

    let truth: &[u8] = b"Back to sequences: find the origin of kmers

Usage: back_to_sequences [OPTIONS] --in-kmers <IN_KMERS>

Options:
      --in-sequences <IN_SEQUENCES>    Input fasta or fastq [.gz] file containing the original sequences (eg. reads). THe stdin is used if not provided [default: ]
      --in-kmers <IN_KMERS>            Input fasta file containing the original kmers
      --out-sequences <OUT_SEQUENCES>  Output file containing the filtered original sequences (eg. reads). It will be automatically in fasta or fastq format depending on the input file. If not provided, only the in_kmers with their count is output [default: ]
      --out-kmers <OUT_KMERS>          If provided, output text file containing the kmers that occur in the reads with their number of occurrences [default: ]
  -k, --kmer-size <KMER_SIZE>          Size of the kmers to index and search [default: 31]
  -m, --min-threshold <MIN_THRESHOLD>  Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold] Minimal threshold of the ratio  (%) of kmers that must be found in a sequence to keep it (default 0%). Thus by default, if no kmer is found in a sequence, it is not output [default: 0]
      --max-threshold <MAX_THRESHOLD>  Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold] Maximal threshold of the ratio (%) of kmers that must be found in a sequence to keep it (default 100%). Thus by default, there is no limitation on the maximal number of kmers found in a sequence [default: 100]
      --stranded                       Used original kmer strand (else canonical kmers are considered)
      --query-reverse                  Query the reverse complement of reads. Useless without the --stranded option
      --no-low-complexity              Do not index low complexity kmers (ie. with a Shannon entropy < 1.0)
  -h, --help                           Print help
  -V, --version                        Print version
";

    let assert = cmd.assert();

    assert.success().stdout(truth);

    Ok(())
}

#[test]
fn argument_trouble() -> std::result::Result<(), anyhow::Error> {
    let mut rng = biotest::rand();
    let s_generate = biotest::Fasta::builder().build()?;
    let k_generate = biotest::Fasta::builder().sequence_len(31).build()?;

    let temp_dir = tempfile::tempdir()?;
    let temp_path = temp_dir.path();
    let kmers_in_path = temp_path.join("kmers_in.fasta");
    let kmers_out_path = temp_path.join("kmers_out.fasta");

    let mut reads = vec![];
    s_generate.records(&mut reads, &mut rng, 1)?;
    k_generate.create(&kmers_in_path, &mut rng, 500)?;

    let mut cmd = assert_cmd::Command::cargo_bin("back_to_sequences")?;
    cmd.args([
        "--in-kmers",
        &kmers_in_path
            .clone()
            .into_os_string()
            .into_string()
            .unwrap(),
    ]);

    let assert = cmd
        .assert()
        .stderr("Warning: no output file provided, nothing to do\n");

    assert.failure();

    let mut cmd = assert_cmd::Command::cargo_bin("back_to_sequences")?;
    cmd.args([
        "--in-kmers",
        &kmers_in_path
            .clone()
            .into_os_string()
            .into_string()
            .unwrap(),
        "--out-kmers",
        &kmers_out_path
            .clone()
            .into_os_string()
            .into_string()
            .unwrap(),
        "--query-reverse",
    ])
    .write_stdin(reads.clone());

    let assert = cmd
        .assert()
        .stderr("Warning: --query-reverse is useless without --stranded\n");

    assert.success();

    let mut cmd = assert_cmd::Command::cargo_bin("back_to_sequences")?;
    cmd.args([
        "--in-kmers",
        &kmers_in_path.into_os_string().into_string().unwrap(),
        "--out-kmers",
        &kmers_out_path.into_os_string().into_string().unwrap(),
        "--max-threshold",
        "1",
        "--min-threshold",
        "2",
    ])
    .write_stdin(reads);

    let assert = cmd
        .assert()
        .stderr("Error: --min-threshold must be <= --max-threshold\n");

    assert.failure();

    Ok(())
}

#[test]
fn default_fasta() -> std::result::Result<(), anyhow::Error> {
    let mut cmd = assert_cmd::Command::cargo_bin("back_to_sequences")?;
    let mut rng = biotest::rand();
    let s_generate = biotest::Fasta::builder().build()?;
    let k_generate = biotest::Fasta::builder().sequence_len(10).build()?;

    let temp_dir = tempfile::tempdir()?;
    let temp_path = temp_dir.path();
    let kmers_in_path = temp_path.join("kmers_in.fasta");
    let kmers_out_path = temp_path.join("kmers_out.fasta");
    let reads_out_path = temp_path.join("reads_out.fasta");

    let mut reads = vec![];
    s_generate.records(&mut reads, &mut rng, 100)?;
    k_generate.create(&kmers_in_path, &mut rng, 500)?;

    cmd.args([
        "-k",
        "10",
        "--in-kmers",
        &kmers_in_path.into_os_string().into_string().unwrap(),
        "--out-kmers",
        &kmers_out_path
            .clone()
            .into_os_string()
            .into_string()
            .unwrap(),
        "--out-sequences",
        &reads_out_path
            .clone()
            .into_os_string()
            .into_string()
            .unwrap(),
    ])
    .write_stdin(reads);

    let out = format!(
        "Indexed 499 kmers, each of size 10
Filtered sequences with exact kmer count are in file {}
kmers with their number of occurrences in the original sequences are in file {}
",
        reads_out_path
            .clone()
            .into_os_string()
            .into_string()
            .unwrap(),
        kmers_out_path
            .clone()
            .into_os_string()
            .into_string()
            .unwrap()
    );

    let assert = cmd.assert().stdout(out);

    assert.success();

    // check reads output
    let mut reads_out_content = vec![];
    std::fs::File::open(&reads_out_path)?.read_to_end(&mut reads_out_content)?;

    let mut reads_out_truth = vec![];
    std::fs::File::open("tests/data/reads_out.fasta")?.read_to_end(&mut reads_out_truth)?;

    assert_eq!(reads_out_content, reads_out_truth);

    // check kmers output
    let mut kmers_out_content = vec![];
    std::fs::File::open(&kmers_out_path)?.read_to_end(&mut kmers_out_content)?;
    kmers_out_content.sort_unstable();

    let mut kmers_out_truth = vec![];
    std::fs::File::open("tests/data/kmers_out.csv")?.read_to_end(&mut kmers_out_truth)?;
    kmers_out_truth.sort_unstable();

    assert_eq!(kmers_out_content, kmers_out_truth);

    Ok(())
}
