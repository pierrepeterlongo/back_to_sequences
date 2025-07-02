use fastbloom::BloomFilter;
use std::f64::consts::LN_2;
use minimizer_iter::MinimizerBuilder;

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
        // BloomFilter::with_size_and_hashers(bloom_filter_size as usize,1, hasher_1, hasher_2);
        KmerPrefiltration {
            k,
            msize,
            wsize: k - msize + 1,
            bloom_filter_size,
            false_positive_rate,
            bloom_filter,
        }
    }


    /// Creates a new KmerPrefiltration from keys of a hashmap 
    pub fn from_kmer_set(kmer_set: &[Vec<u8>], false_positive_rate: f32, k: usize, msize: usize) -> Self {
        let mut prefilter = KmerPrefiltration::new(kmer_set.len() as u32, false_positive_rate, k, msize);
        prefilter.insert_kmer_set(kmer_set);
        prefilter
    }

    

    fn kmer_to_minimizer(&self, kmer: &[u8]) -> u64 {
        let mut min_iter = MinimizerBuilder::<u64>::new()
            .canonical()
            .minimizer_size(self.msize.into())
            .width(self.wsize as u16)
            .iter(kmer);


        // check that we have exactly one minimizer
        if let Some((minimizer, _position, _forward)) = min_iter.next() {
            if let Some((_minimizer, _position, _forward)) = min_iter.next() {
                // debug
                panic!("Panic: kmer {} has more than one minimizer", String::from_utf8(kmer.to_vec()).unwrap());
            }
            return minimizer;
        } else {
            panic!("Panic: kmer {} has no minimizer", String::from_utf8(kmer.to_vec()).unwrap());
        }
    }

    /// Returns true if the kmer is in the bloom filter.
    pub fn contains(&self, kmer: &[u8]) -> bool {
        self.bloom_filter.contains(&self.kmer_to_minimizer(kmer))
    }

    /// Inserts a kmer into the bloom filter.
    fn insert(&mut self, kmer: &[u8]) {
        self.bloom_filter.insert(&self.kmer_to_minimizer(kmer));
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
        let mut min_iter = MinimizerBuilder::<u64, _>::new()
            .canonical()
            .minimizer_size(self.msize.into())
            .width(self.wsize as u16)
            .iter(&sequence);

        let mut returned_positions = Vec::new();

        // We have a list of position of minimizers
        // We need to find which kmers are assigne to each minimizer. 
        // Isolated minimizer: Given a minimizer position p, with no other minimizer around, kmers from position 
        // p-k+m to p are assigned to this minimizer.
        // General case: Given a minimizer position p1 and next minimizer position p2, kmers from position
        // p1-k+m to p2-k+m-1 are assigned to minimizer position p1. 
        // We report only kmers whose minimizers are "valid" (in the bf). 
        //      Here is the algo:
        // (current_minimizer, current_minimizer_position) = first couple 
        // is_current_minimizer_valid = is current_minimizer in bf
        // (next_minimizer, next_minimizer_position) = second couple 
        // for i from 0 to len sequence - k (included) :
        //     if i >= next_minimizer_position - k + m: // we change the minimizer: 
        //           (current_minimizer, current_minimizer_position) = (next_minimizer, next_minimizer_position)
        //           is_current_minimizer_valid = is current_minimizer in bf
        //           (next_minimizer, next_minimizer_position) = next couple 
        //    
        //     if is_current_minimizer_valid:
        //           add i to return values
        
        let shift_shortcut = self.k - self.msize; // shortcut
        let  (mut current_minimizer, mut current_minimizer_position, _forward) = 
            min_iter.next().unwrap_or((0 as u64, sequence.len(), false));
        let mut is_current_minimizer_valid = self.bloom_filter.contains(&current_minimizer);
        let (mut next_minimizer, mut next_minimizer_position, mut _forward) = min_iter.next().unwrap_or((0 as u64, sequence.len(), false)); // in case we have no more minimizer, we consider the last minimizer as occurring at sequence length. This fills all kmers. 
        for i in 0..sequence.len() - self.k + 1usize {
                        // TODO: hugly : I'd like to make only one comparison.
            // if next_minimizer_position - shift_shortcut is negative,  sustraction creates huge positive result 
            if i as isize >= (next_minimizer_position - shift_shortcut) as isize {
                current_minimizer = next_minimizer; // moved value
                is_current_minimizer_valid = self.bloom_filter.contains(&current_minimizer);
                current_minimizer_position = next_minimizer_position;
                (next_minimizer, next_minimizer_position, _forward) = min_iter.next().unwrap_or((0 as u64, sequence.len(), false));
                // new minimizer should not be at a lower position than previous one.
                assert!(current_minimizer_position < next_minimizer_position); 
            }
            if is_current_minimizer_valid {
                returned_positions.push(i);
            }
        }


        // // no need to store this vector, we can simply iterate over minimizer positions
        // let  (mut current_minimizer, mut current_minimizer_position, _forward) = min_iter.next().unwrap();
        // let mut current_minimizer_in_bf = self.bloom_filter.contains(&current_minimizer);
        // let mut next_iterator = min_iter.next();
        // let (mut next_minimizer, mut next_minimizer_position) = if let Some((min, pos, _forward)) = next_iterator {
        //     (Some(min), Some(pos))
        // } else {
        //     (None, None)
        // };
        // // TODO: l'algo est faux...
        // let mut current_position = 0;
        // loop {
        //     //DEBUG
        //     println!("Current minimizer: {}, position: {}, in bloom filter: {}", current_minimizer, current_minimizer_position, current_minimizer_in_bf);
        //     if current_minimizer_in_bf {
        //         // fill from current_position to min(current_minimizer_position +1, len(sequence) - k )
        //         let end_position = current_minimizer_position.min(sequence.len() - self.k);
        //         // if current_position is already past end_position, we can skip
        //         for  position in current_position..end_position + 1 {
        //             returned_positions.push(position);
        //         }
        //     }
        //     // if next minimizer is not defined, we can stop
        //     if next_minimizer_position == None{
        //         break;
        //     }


        //     // update current and next minimizer for next loop 
        //     current_minimizer = next_minimizer.unwrap();
        //     current_minimizer_position = next_minimizer_position.unwrap();
        //     current_minimizer_in_bf = self.bloom_filter.contains(&current_minimizer);
        //     current_position = current_minimizer_position - self.msize as usize + self.k;

        //     next_iterator = min_iter.next();
        //     if let Some((min, pos, _forward)) = next_iterator {
        //         next_minimizer = Some(min);
        //         next_minimizer_position = Some(pos);
        //     } else {
        //         next_minimizer = None;
        //         next_minimizer_position = None;
        //     }
        // }
        returned_positions
    }


    /// Display
    pub fn display(&self) {
        println!("KmerPrefiltration: k: {}, msize: {}, wsize: {}, bloom_filter_size: {}, false_positive_rate: {}", 
            self.k, self.msize, self.wsize, self.bloom_filter_size, self.false_positive_rate);
    }
}
