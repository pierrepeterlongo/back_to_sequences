//! Fasta generator

/* std use */

/* crates use */

/* project use */

pub fn comment(rng: &mut rand::rngs::StdRng, length: usize) -> Vec<u8> {
    let mut result = Vec::with_capacity(length + 8);

    result.extend(b">read_");
    result.extend(super::name(rng, length));
    result.push(b'\n');

    result
}

pub fn record(rng: &mut rand::rngs::StdRng, comment_length: usize, read_length: usize) -> Vec<u8> {
    let mut result = Vec::with_capacity(comment_length + 7 + read_length + 2);

    result.extend(comment(rng, comment_length));
    result.extend(super::sequence(rng, read_length));
    result.push(b'\n');

    result
}

pub fn records(
    rng: &mut rand::rngs::StdRng,
    comment_length: usize,
    read_length: usize,
    n: usize,
) -> Vec<u8> {
    let mut result = Vec::with_capacity((comment_length + 7 + read_length + 2) * n);

    for _ in 0..n {
        result.extend(record(rng, comment_length, read_length));
    }

    result
}
