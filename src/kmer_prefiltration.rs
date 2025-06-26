use bloom::{ASMS,BloomFilter};
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
    /// w size
    wsize: usize,
    /// The id of the read.
    pub expected_num_items: u32,
    /// The position of the match in the read.
    pub false_positive_rate: f32,
    /// Whether the match is forward.
    pub bloom_filter: BloomFilter,
}

impl KmerPrefiltration {
    /// Creates a new KmerPrefiltration with the given expected number of items and false positive rate.
    pub fn new(expected_num_items: u32, false_positive_rate: f32, wsize: usize) -> Self {
        let bloom_filter = BloomFilter::with_rate(false_positive_rate,expected_num_items);
        KmerPrefiltration {
            wsize,
            expected_num_items,
            false_positive_rate,
            bloom_filter,
        }
    }

    /// Creates a new KmerPrefiltration from keys of a hashmap 
    pub fn from_kmer_set(kmer_set: &[Vec<u8>], false_positive_rate: f32, wsize: usize) -> Self {
        let mut prefilter = KmerPrefiltration::new(kmer_set.len() as u32, false_positive_rate, wsize);
        prefilter.insert_kmer_set(kmer_set);
        prefilter
    }

    fn kmer_to_minimizer(&self, kmer: &[u8]) -> Vec<u8> {
        // For simplicity, we just return the first wsize nucleotides as a minimizer.
        // In practice, you might want to implement a more sophisticated minimizer selection.
        // return a vector of the first wsize nucleotides of the kmer
        kmer.iter().take(self.wsize).cloned().collect()
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
}
