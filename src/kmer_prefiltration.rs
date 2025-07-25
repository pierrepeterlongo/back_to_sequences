use fastbloom::BloomFilter;
use std::f64::consts::LN_2;
use simd_minimizers;
use packed_seq::{PackedSeqVec, SeqVec};

fn replace_non_acgt(input: &[u8]) -> Vec<u8> {
    input
        .iter()
        .map(|&c| {
            match c {
                b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' => c,
                _ => b'A',
            }
        })
        .collect()
}


/// kmer prefiltration. 
/// Key ideas: 
///  Building: in a bloom filter, instert firt $w$ nucleotides of each indexed kmer
///  Querying: returns true if first $w$ nucleorides are in the bloom filter
/// Todo : use a minimizer instead of first values
/// 
/// 
/// 
///
pub struct KmerPrefiltration {
    /// kmer size
    k: usize,
    /// minimizer size
    msize: usize,
    /// Minimizer window size
    wsize: usize, // wsize = k - msize + 1
    /// The id of the read.
    bloom_filter_size: usize,
    /// The position of the match in the read.
    false_positive_rate: f32,
    /// Whether the match is forward.
    bloom_filter: BloomFilter,
}

impl KmerPrefiltration {
    /// Creates a new KmerPrefiltration with the given expected number of items and false positive rate.
    /// We consider a unique hash function.
    /// Number of bits in the bloom filter is calculated as: m= -n * ln(p) / ln(2)^2
    /// where n is the expected number of items and p is the false positive rate.
    pub fn new(expected_num_items: u32, false_positive_rate: f32, k: usize, msize: usize) -> Self {
        let bloom_filter_size = (-(expected_num_items as f64 * false_positive_rate.ln() as f64) / LN_2.powi(2)) as usize +1usize; // +1 to avoid zero size


        let bloom_filter = BloomFilter::with_num_bits(bloom_filter_size).hashes(1);

        // wsize cannot be < 1:
        let wsize = if k < msize {
            panic!("k cannot be less than msize");
        } else {
            k - msize + 1
        };
        KmerPrefiltration {
            k,
            msize,
            wsize,
            bloom_filter_size,
            false_positive_rate,
            bloom_filter,
        }
    }


    /// Creates a new KmerPrefiltration from keys of a hashmap 
    pub fn from_kmer_set(kmer_set: &[Vec<u8>], false_positive_rate: f32, k: usize, msize: usize) -> Self {
        
        // DEBUG
        // let seq = b"GAGCCGGAAGCGTGGCTGCGTAACGATCACC";
        // let m = 19;
        // let w = 11;

        // let packed_seq = PackedSeqVec::from_ascii(seq);
        // let mut minimizer_positions = Vec::new();
        // simd_minimizers::canonical_minimizer_positions(packed_seq.as_slice(), m, w, &mut minimizer_positions);
        // println!("Minimizer positions {:#?}", minimizer_positions);
        // std::process::exit(0);
        // let _minimizer_values: Vec<_> = simd_minimizers::iter_canonical_minimizer_values(packed_seq.as_slice(), k, &minimizer_positions).collect();

        // FIN DEBUG

        let mut prefilter = KmerPrefiltration::new(kmer_set.len() as u32, false_positive_rate, k, msize);
        prefilter.insert_kmer_set(kmer_set);
        prefilter
    }

    

    fn kmer_to_minimizer(&self, kmer: &[u8]) -> u64 {
        // simd_minimizers::one_canonical_minimizer(kmer, self.msize)
        let mut minimizer_positions = Vec::new();
        let mut superkmer_pos_vec = Vec::new();
        let packed_seq = PackedSeqVec::from_ascii(&replace_non_acgt(kmer));
        simd_minimizers::scalar::canonical_minimizer_and_superkmer_positions_scalar(packed_seq.as_slice(), self.msize, self.wsize, &mut minimizer_positions, &mut superkmer_pos_vec);
        // println!("number of minimizers {} for kmer {:?}",  minimizer_positions.len(), str::from_utf8(kmer));        
        assert!(minimizer_positions.len() == 1);

        // println!("kmer {:?} has minimizer positions {:?}", String::from_utf8(kmer.to_vec()).unwrap(), minimizer_positions);
        for canonical_minimizer_value in simd_minimizers::iter_canonical_minimizer_values(packed_seq.as_slice(), self.msize, &minimizer_positions){
            // println!("kmer {:?} has minimizer value {}", String::from_utf8(kmer.to_vec()).unwrap(), canonical_minimizer_value);
            return canonical_minimizer_value;
        }
        0
    }

    /// Returns true if the kmer is in the bloom filter.
    pub fn contains(&self, kmer: &[u8]) -> bool {
        self.bloom_filter.contains(&self.kmer_to_minimizer(kmer))
    }

    /// Inserts a kmer into the bloom filter.
    fn insert(&mut self, kmer: &[u8]) {
        self.bloom_filter.insert(&self.kmer_to_minimizer(kmer));
        // println!("inserted minimizer {} of kmer {} in the bf", &self.kmer_to_minimizer(kmer), String::from_utf8(kmer.to_vec()).unwrap());
    }

    /// Insert all minimizers from a kmer set into the bloom filter.
    pub fn insert_kmer_set(&mut self, kmer_set: &[Vec<u8>]) {
        for kmer in kmer_set {
            // transform kmer to [u8]
            self.insert(&kmer.as_slice());  
        }
    }

    /// Given a sequence, return a vector containing kmer positions whose
    /// minimizer is in the bloom filter.
    pub fn potiential_kmer_positions(&self, sequence: &[u8]) -> Vec<usize> {
        let mut minimizer_positions = Vec::new();
        let mut superkmer_pos_vec = Vec::new();
        let packed_seq = PackedSeqVec::from_ascii(&replace_non_acgt(sequence));
        simd_minimizers::scalar::canonical_minimizer_and_superkmer_positions_scalar(packed_seq.as_slice(), self.msize, self.wsize, &mut minimizer_positions, &mut superkmer_pos_vec);
        superkmer_pos_vec.push(sequence.len() as u32); // add a dummy value to avoid out of bounds access
        // println!("Minimizer positions: {:?} for seq {:#?}", minimizer_positions, str::from_utf8(&replace_non_acgt(sequence)));
        // println!("Superkmer positions: {:?}", superkmer_pos_vec);
        let last_position = (sequence.len() - self.k + 1) as u32;

        let mut returned_positions = Vec::new();
        for (i, canonical_minimizer_value) in simd_minimizers::iter_canonical_minimizer_values(packed_seq.as_slice(), self.msize, &minimizer_positions).enumerate() {
            // println!("Tested minimizer value: {}", canonical_minimizer_value);
            // kmers from superkmer_pos_vec[i] to superkmer_pos_vec[i+1] are assigned 
            // to minimizer canonical_minimizer_value
            // println!("Testing {}th  minimizer value: {}", i, canonical_minimizer_value);
            if  self.bloom_filter.contains(&canonical_minimizer_value) {
                
                let start = superkmer_pos_vec[i];
                let end = std::cmp::min(
                    std::cmp::min(
                    superkmer_pos_vec[i] + (2 * self.k - self.msize) as u32,
                    superkmer_pos_vec[i + 1]),
                    last_position as u32
                );

                for p in start..end {
                    returned_positions.push(p as usize);
                }
                // println!("Found {}th minimizer value: {} in the bloom filter, positions: {} to {}", i, canonical_minimizer_value, start, end);
            }
            // else {
            //     println!("Did not find {}th minimizer value: {} in the bloom filter", i, canonical_minimizer_value);
            // }
        }
        returned_positions
    }


    /// Display
    pub fn display(&self) {
        println!("KmerPrefiltration: k: {}, msize: {}, wsize: {}, bloom_filter_size: {}, false_positive_rate: {}", 
            self.k, self.msize, self.wsize, self.bloom_filter_size, self.false_positive_rate);
    }
}
