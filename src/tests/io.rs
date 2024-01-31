//! Manage input/output

/* std use */

/* crates use */
use std::io::Write;

/* project use */
use super::fasta;
use super::fastq;

pub fn write_buffer<P>(buffer: &[u8], path: P) -> std::io::Result<()>
where
    P: AsRef<std::path::Path>,
{
    let mut file = std::io::BufWriter::new(std::fs::File::create(&path)?);

    file.write_all(buffer)?;

    Ok(())
}

/// Write a random fasta in path
pub fn write_fasta<P>(
    rng: &mut rand::rngs::StdRng,
    seq_length: usize,
    comment_length: usize,
    seq_number: usize,
    path: P,
) -> std::io::Result<()>
where
    P: AsRef<std::path::Path>,
{
    let data = fasta::records(rng, comment_length, seq_length, seq_number);
    println!("{:?}", String::from_utf8(data.clone()));

    write_buffer(&data, path)
}

/// Write a random fastq in path
pub fn write_fastq<P>(
    rng: &mut rand::rngs::StdRng,
    seq_length: usize,
    comment_length: usize,
    seq_number: usize,
    path: P,
) -> std::io::Result<()>
where
    P: AsRef<std::path::Path>,
{
    write_buffer(
        &fastq::records(rng, comment_length, seq_length, seq_number),
        path,
    )
}
