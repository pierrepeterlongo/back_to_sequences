//! Declarations of test utils

/* std use */

/* crates use */
use rand::seq::SliceRandom;
use rand::SeedableRng;

/* mod declarations */
pub mod constants;
pub mod fasta;
pub mod fastq;
pub mod io;

/* project use */

pub fn rand() -> rand::rngs::StdRng {
    rand::rngs::StdRng::from_seed(constants::SEED)
}

pub fn name(rng: &mut rand::rngs::StdRng, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            *constants::ALPHABETS
                .choose(rng)
                .unwrap_or_else(|| unreachable!())
        })
        .collect()
}

pub fn sequence(rng: &mut rand::rngs::StdRng, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            *constants::NUCLEOTIDES
                .choose(rng)
                .unwrap_or_else(|| unreachable!())
        })
        .collect()
}

pub fn quality(rng: &mut rand::rngs::StdRng, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            *constants::PHRED33
                .choose(rng)
                .unwrap_or_else(|| unreachable!())
        })
        .collect()
}
