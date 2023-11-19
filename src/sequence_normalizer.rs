//! Sequence normalization

/* std use */

/* crates use */

/* project use */
use crate::consts::FORWARD_MAP;
use crate::consts::REVERSE_MAP;

/// Zero-copy object for normalizing a sequence
///
/// - optionally reverse complements the sequence
///
/// To avoid extra memory allocations, the original slice is kept in place and the resulting
/// sequence is returned through an iterator (`fn iter()`).
///
pub struct SequenceNormalizer<'a> {
    raw: &'a [u8],
    reverse_complement: bool,
}

impl<'a> SequenceNormalizer<'a> {
    /// - `raw` is the original sequence (ascii string)
    /// - `reverse_complement` may be:
    ///     - `Some(false)` to get the original sequence
    ///     - `Some(true)` to get the reverse complement
    ///     - `None` to get the canonical sequence
    pub fn new(raw: &'a [u8], reverse_complement: Option<bool>) -> Self {
        Self {
            raw,
            reverse_complement: reverse_complement.unwrap_or_else(|| {
                let forward = Self::iter_impl(raw, false);
                let reverse = Self::iter_impl(raw, true);
                reverse.cmp(forward).is_lt()
            }),
        }
    }

    #[auto_enums::auto_enum(Iterator)]
    fn iter_impl(raw: &[u8], reverse_complement: bool) -> impl Iterator<Item = u8> + '_ {
        if reverse_complement {
            raw.iter().rev().map(|c| {
                REVERSE_MAP[*c as usize]
                    .expect("cannot complement base (contains non-ascii byte: 0x{:x})")
                    .into()
            })
        } else {
            raw.iter().map(|c| FORWARD_MAP[*c as usize])
        }
    }

    /// Get an iterator on the normalized sequence
    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        Self::iter_impl(self.raw, self.reverse_complement)
    }

    /// Copy the normalized sequence into a slice
    ///
    /// panics for the slice has a different length
    pub fn copy_to_slice(&self, dest: &mut [u8]) {
        assert_eq!(dest.len(), self.raw.len());
        for (i, c) in self.iter().enumerate() {
            dest[i] = c;
        }
    }
}
