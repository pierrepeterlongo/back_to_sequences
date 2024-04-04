//! Declarations of constants value

/* std use */

/* crates use */

/* project use */

const fn gen_array<const N: usize, const B: usize>() -> [u8; N] {
    let mut array = [0; N];

    let mut i = 0;
    while i < N {
        array[i] = (B + i) as u8;
        i += 1;
    }

    array
}

pub const SEED: [u8; 32] = [42; 32];

pub const NUCLEOTIDES: [u8; 8] = *b"ACTGactg";
pub const PHRED33: [u8; 41] = gen_array::<41, 33>();
pub const _PHRED64: [u8; 40] = gen_array::<40, 64>();
pub const ALPHABETS: [u8; 26] = gen_array::<26, 97>();
