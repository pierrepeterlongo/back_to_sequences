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
    msize: u16,
    /// Minimizer window size
    wsize: u16, // wsize = k - msize + 1
    /// The id of the read.
    pub bloom_filter_size: usize,
    /// The position of the match in the read.
    pub false_positive_rate: f32,
    /// Whether the match is forward.
    pub bloom_filter: BloomFilter,
}

impl KmerPrefiltration {
    /// Creates a new KmerPrefiltration with the given expected number of items and false positive rate.
    /// We consider a unique hash function.
    /// Number of bits in the bloom filter is calculated as: m= -n * ln(p) / ln(2)^2
    /// where n is the expected number of items and p is the false positive rate.
    pub fn new(expected_num_items: u32, false_positive_rate: f32, k: usize, msize: u16) -> Self {
        let bloom_filter_size = (-(expected_num_items as f64 * false_positive_rate.ln() as f64) / LN_2.powi(2)) as usize;

        let bloom_filter = BloomFilter::with_num_bits(bloom_filter_size).hashes(1);
        // BloomFilter::with_size_and_hashers(bloom_filter_size as usize,1, hasher_1, hasher_2);
        KmerPrefiltration {
            k,
            msize,
            wsize: k as u16 - msize + 1,
            bloom_filter_size,
            false_positive_rate,
            bloom_filter,
        }
    }

    /// Creates a new KmerPrefiltration from keys of a hashmap 
    pub fn from_kmer_set(kmer_set: &[Vec<u8>], false_positive_rate: f32, k: usize, msize: u16) -> Self {
        let mut prefilter = KmerPrefiltration::new(kmer_set.len() as u32, false_positive_rate, k, msize);
        prefilter.insert_kmer_set(kmer_set);
        prefilter
    }

    fn kmer_to_minimizer(&self, kmer: &[u8]) -> u64 {


        // Build an iterator over minimizers
        // of size 21 with a window of size 11
        // for the sequence "TGATTGCACAATC"
        let min_iter = MinimizerBuilder::<u64>::new()
            .minimizer_size(self.msize.into())
            .width(self.wsize)
            .iter(kmer);
        // todo check that we have exactly one minimizer
        // if min_iter.count() != 1 {
        //     panic!("Panic: kmer {} has no minimizer", String::from_utf8(kmer.to_vec()).unwrap());
        // }
        for (minimizer, _position) in min_iter {
            return minimizer
        }
        // check that this code is unreachable
        unreachable!()
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
            self.insert(kmer);  
        }
    }

    /// Given a sequence, return an iterator over kmer positions whose
    /// minimizer is in the bloom filter.
    pub fn potiential_kmer_positions(&mut self, sequence: &[u8]) -> Vec<usize> {
        let min_iter = MinimizerBuilder::<u64, _>::new_mod()
            .minimizer_size(self.msize.into())
            .width(self.wsize)
            .iter(sequence);
        
        let mut returned_positions = Vec::new();

        // assign each position i to a minimizer position p
        // such that the kmer starting position i as minimizer position p
        // easy: 
        //  if a minimizer is "isolated" (no other minimizers in [p - m +k..p+k])
        //  then all kmers in [p - m + k, p] are assigned to p
        //  else, a new minimizer enters the window position p'>p.
        //        as a minimizer, its value is smaller than the previous minimizer
        //        so all kmers in [p' - m + k, p'] are
        //        assigned to p'
        // thus we can assign each kmer position to a minimizer position by storing: 
        // a vector get_minimizer_position containing minimizer positions - m + k (or zero if <0)
        // then reading this vector from left to right indicates when we change minimizer

        // no need to store this vector, we can simply iterate over minimizer positions
        let  (mut current_minimizer, mut current_minimizer_position) = min_iter.next().unwrap();
        let mut current_minimizer_in_bf = self.bloom_filter.contains(&current_minimizer);
        let mut next_iterator = min_iter.next();
        let mut next_minimizer;
        let mut next_minimizer_position;
        if next_iterator == None {
            next_minimizer = None;
            next_minimizer_position = None;
        }
        else {
            (next_minimizer, next_minimizer_position) = next_iterator.unwrap();
        }
        let mut current_position = 0;
        loop {
            // if the current minimizer is not in the bloom filter, we can skip it
            if current_minimizer_in_bf {
                let mut stop_range = sequence.len() + 1 as u16; 
                if {next_minimizer_position != None}{
                    stop_range = next_minimizer_position;
                }
                for  position in current_position..stop_range {
                    returned_positions.push(position);
                }
            }
            // if next minimizer is not defined, we can stop
            if next_minimizer_position == None{
                break;
            }


            // update current and next minimizer
            current_minimizer = next_minimizer;
            current_minimizer_position = next_minimizer_position;
            current_minimizer_in_bf = self.bloom_filter.contains(&current_minimizer);
            current_position = next_minimizer_position - self.msize as u16 + self.k;

            next_iterator = min_iter.next();
            if next_iterator == None {
                next_minimizer = None;
                next_minimizer_position = None;
            }
            (next_minimizer, next_minimizer_position) = next_iterator.unwrap();
               
        }
        returned_positions   
    }
}
