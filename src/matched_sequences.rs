//! matched_sequences declarations

use std::fmt;

/// round a float to a given number of decimals
fn round(x: f32, decimals: u32) -> f32 {
    let y = 10i32.pow(decimals) as f32;
    (x * y).round() / y
}

/// a read matched by a kmer
pub trait MatchedSequence
where
    Self: Sized + std::fmt::Display,
{
    /// Initialize the MatchedSequence
    fn new(mapped_position_size: usize) -> Self;

    /// add a match to the read
    fn add_match(&mut self, position: usize, forward: bool);

    // /// prints the matched read
    // fn to_string(&self) -> String;

    /// returns the percentage of the read that was matched
    fn percent_shared_kmers(&self) -> f32;
}

/// a read matched by a kmer, only counting hte number of matched kmers
pub struct MachedCount {
    /// size of the positions where a match can occur (size_seq - k +1)
    pub mapped_position_size: usize,
    /// number of matched kmers
    pub count: usize,
}

impl fmt::Display for MachedCount {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            " {} {}",
            self.count,
            round(self.percent_shared_kmers(), 5)
        )
    }
}

/// implement MatchedSequence for the very simple Matched count read
impl MatchedSequence for MachedCount {
    fn new(mapped_position_size: usize) -> Self {
        MachedCount {
            mapped_position_size,
            count: 0,
        }
    }

    fn add_match(&mut self, _position: usize, _forward: bool) {
        self.count += 1;
    }

    fn percent_shared_kmers(&self) -> f32 {
        100.0 * self.count as f32 / (self.mapped_position_size as f32)
    }
}

/// a read matched by a kmer, store kmer count and position of matched kmers
pub struct MatchedSequencePositional {
    /// size of the positions where a match can occur (size_seq - k +1)
    pub mapped_position_size: usize,

    /// number of matched kmers
    pub count: usize,

    /// Position of the matched kmers in the read
    /// the boolean indicates if the kmer was mapped in the forward or reverse strand
    pub matched_positions: Vec<(usize, bool)>,
}

impl MatchedSequence for MatchedSequencePositional {
    fn new(mapped_position_size: usize) -> Self {
        MatchedSequencePositional {
            mapped_position_size,
            count: 0,
            matched_positions: Vec::new(),
        }
    }

    fn add_match(&mut self, position: usize, forward: bool) {
        self.count += 1;
        self.matched_positions.push((position, forward));
    }

    // TODO how to avoid to duplicate this code
    fn percent_shared_kmers(&self) -> f32 {
        100.0 * self.count as f32 / (self.mapped_position_size as f32)
    }
}

impl fmt::Display for MatchedSequencePositional {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Pre-allocate vector with estimated capacity
        let capacity = self.matched_positions.len() * 8 + 20; // rough estimate for numbers
        let mut parts = Vec::with_capacity(capacity);
        
        // Add initial count and percent
        parts.push(self.count.to_string());
        parts.push(round(self.percent_shared_kmers(), 5).to_string());
        
        // Store the number of positions covered by the matched kmers
        let mut first_uncovered_position = 0;
        let mut number_covered_positions = 0;

        // Add positions 
        for (position, forward) in &self.matched_positions {
            if first_uncovered_position < *position {
                number_covered_positions += unsafe { crate::KMER_SIZE };
            }
            else {
                number_covered_positions += position + unsafe { crate::KMER_SIZE } - first_uncovered_position;
            }
            first_uncovered_position = position + unsafe { crate::KMER_SIZE };

            if *forward {
                parts.push(position.to_string());
            } else {
                // Avoid allocating a new string for the minus sign
                parts.push(format!("-{}", position));
            }
        }
        
        // Join all parts with spaces
        write!(f, " {} ({})", parts.join(" "), number_covered_positions)

    }
}

#[cfg(test)]
mod tests {
    use crate::KMER_SIZE;

    use super::*;
    #[test]
    fn simple_match() {
        unsafe { KMER_SIZE = 5 };
        let sequence = b"ACGTGACTACGGCATAGCATCGTAGCTGATGTGTCAGCTGTCAGTCA";
        let mut mc = MachedCount::new(sequence.len() - unsafe { KMER_SIZE } + 1);
        mc.add_match(4, true);
        mc.add_match(5, false);
        mc.add_match(6, false);

        assert_eq!(mc.to_string(), " 3 6.97674");
    }

    #[test]
    fn positional_match() {
        unsafe { KMER_SIZE = 5 };
        let sequence = b"ACGTGACTACGGCATAGCATCGTAGCTGATGTGTCAGCTGTCAGTCA";
        let mut mc = MatchedSequencePositional::new(sequence.len() - unsafe { KMER_SIZE } + 1);
        mc.add_match(4, true);
        mc.add_match(5, false);
        mc.add_match(6, false);

        assert_eq!(mc.to_string(), " 3 6.97674 4 -5 -6 (7)");
    }
}
