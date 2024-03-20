//! various counters, only the counts or loactions in queried sequences

use integer_encoding::*;
use atomic_counter::AtomicCounter as _;

use std::sync::Mutex;

/* project use */
// use crate::kmer_counter::KmerCounter;

/// Trait representing a KmerCounter.
pub trait KmerCounter {
    /// Adds a match to the counter.
    fn add_match(&self, id_read: usize, position: usize, stranded: bool) -> ();

    /// Initializes the counter. (TODO, validate with anthony)
    fn default(&self) -> ();

    /// Returns a string representation of the counter. Anthony ? right way to do this ? 
    fn to_string(&self) -> String;
}
impl KmerCounter for atomic_counter::RelaxedCounter {
    /// Adds a match to the counter.
    fn add_match(&self, _id_read: usize, _position: usize, _stranded: bool) -> () {
        self.inc();
    }

    fn default(&self) -> () { // TODO : validate with anthony
        self.reset();
    }
    
    fn to_string(&self) -> String {
        self.get().to_string() // TODO : validate with anthony
    }
}

#[derive(Default)]
pub(crate) struct KmerCounterWithLog {
    count: usize,
    log: Vec<u8>,
}

impl KmerCounter for Mutex<KmerCounterWithLog>
{
    /// Adds a match to the counter, storing the id_read and position in the log.
    fn add_match(&self, id_read: usize, position: usize, stranded: bool) -> () {
        let mut counter = self.lock().unwrap();
        counter.count += 1;
        counter.log.write_varint(id_read).unwrap();
        counter.log.write_varint(if stranded {-(position as i64)} else {position as i64}).unwrap();
    }

    /// Initializes the counter. TODO: validate with anthony
    fn default(&self) -> () {
        let mut counter = self.lock().unwrap();
        counter.count = 0;
        counter.log.clear();
    }

    fn to_string(&self) -> String {
        let counter = self.lock().unwrap();
        let mut cursor = std::io::Cursor::new(&counter.log);
        let mut result = String::new(); // Create an empty string to store the result
        loop {
            let id_read = match cursor.read_varint::<usize>() {
                Ok(v) => v,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                _ => unreachable!(),
            };
            let (position, stranded) = match cursor.read_varint::<usize>() {
                Ok(v) if v < 0 => (-(v as i64), true),
                Ok(v) if v >= 0 => ((v as i64), false),
                _ => unreachable!(),
            };
            // Append the id_read, position, and stranded information to the result string. Change the format when debugged
            result.push_str(&format!("id_read: {}, position: {}, stranded: {}\n", id_read, position, stranded));
        }
        result // Return the result string
    }
}
