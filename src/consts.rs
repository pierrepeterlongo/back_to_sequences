//! Constant declarations

/* std use */

/* crates use */

/* project use */

fn base_complement(base: u8) -> Option<std::num::NonZeroU8> {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        c if c < 128 => c,
        // all non-ascii character yield None so that we abord in case we try to reverse complement
        // a multibyte unicode character.
        _ => 0,
    }
    .try_into()
    .ok()
}

/// forward lookup-table for the SequenceNormalizer
#[ctor::ctor]
pub static FORWARD_MAP: [u8; 256] = {
    let mut a = [0; 256];
    a.iter_mut()
        .enumerate()
        .for_each(|(i, c)| *c = (i as u8).to_ascii_uppercase());
    a
};
/// reverse lookup-table for the SequenceNormalizer
#[ctor::ctor]
pub static REVERSE_MAP: [Option<std::num::NonZeroU8>; 256] = {
    let mut a = [None; 256];
    a.iter_mut()
        .enumerate()
        .for_each(|(i, c)| *c = base_complement((i as u8).to_ascii_uppercase()));
    a
};
