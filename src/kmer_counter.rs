//! various counters, only the counts or loactions in queried sequences

use atomic_counter::AtomicCounter as _;
use integer_encoding::*;

use std::sync::Mutex;

/// Information needed to represent a match between a kmer and a read.
pub struct KmerMatch {
    /// The id of the read.
    pub id_read: usize,
    /// The position of the match in the read.
    pub position: usize,
    /// Whether the match is forward.
    pub forward: bool,
}

/// Trait representing a KmerCounter.
pub trait KmerCounter: Default + Sync + Send {
    /// Adds a match to the counter.
    fn add_match(&self, m: KmerMatch);

    /// Returns a string representation of the counter. Anthony ? right way to do this ?
    fn to_string(&self) -> String;

    /// Returns the count of the counter.
    fn get_count(&self) -> usize;
}

impl KmerCounter for atomic_counter::RelaxedCounter {
    /// Adds a match to the counter.
    fn add_match(&self, _m: KmerMatch) {
        self.inc();
    }

    fn to_string(&self) -> String {
        self.get().to_string()
    }

    fn get_count(&self) -> usize {
        self.get()
    }
}

#[derive(Default)]
/// A KmerCounter that stores the id_read and position of each match in a log.
pub struct KmerCounterWithLog {
    count: usize,
    log: Vec<u8>,
}

impl KmerCounterWithLog {
    fn iter_matches(&self) -> impl std::iter::Iterator<Item = (usize, usize, bool)> + '_ {
        let mut cursor = std::io::Cursor::new(&self.log);
        std::iter::from_fn(move || {
            let id_read = match cursor.read_varint::<usize>() {
                Ok(v) => v,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return None,
                _ => unreachable!(),
            };
            let (position, stranded) = match cursor.read_varint::<i64>() {
                Ok(v) if v < 0 => (-v, true),
                Ok(v) if v >= 0 => (v, false),
                _ => unreachable!(),
            };
            Some((id_read, position as usize, stranded))
        })
    }
}

impl KmerCounter for Mutex<KmerCounterWithLog> {
    /// Adds a match to the counter, storing the id_read and position in the log.
    fn add_match(&self, m: KmerMatch) {
        let mut counter = self.lock().unwrap();
        counter.count += 1;
        counter.log.write_varint(m.id_read).unwrap();
        counter
            .log
            .write_varint(if m.forward {
                -(m.position as i64)
            } else {
                m.position as i64
            })
            .unwrap();
    }

    fn to_string(&self) -> String {
        let counter = self.lock().unwrap();
        let mut result = String::new(); // Create an empty string to store the result
        for (id_read, position, stranded) in counter.iter_matches() {
            // Append the id_read, position, and stranded information to the result string. Change the format when debugged
            result.push_str(&format!("({},{},{}) ", id_read, position, stranded));
        }
        result // Return the result string
    }

    fn get_count(&self) -> usize {
        self.lock().unwrap().count
    }
}

// TODO Pierre: impl Display

#[cfg(test)]
mod tests {
    /* project use */
    use super::*;

    #[test]
    fn relaxed_counter() {
        let counter = atomic_counter::RelaxedCounter::new(0);

        assert_eq!(counter.get_count(), 0);
        assert_eq!(counter.to_string(), "0".to_string());

        for i in 1..11 {
            counter.add_match(KmerMatch {
                id_read: 0,
                position: 0,
                forward: true,
            });

            assert_eq!(counter.get_count(), i);
            assert_eq!(counter.to_string(), format!("{}", i));
        }

        assert_eq!(counter.get_count(), 10);
        assert_eq!(counter.to_string(), "10".to_string());
    }

    #[test]
    fn conter_with_log() {
        let object = KmerCounterWithLog {
            count: 0,
            log: vec![],
        };
        let counter = Mutex::new(object);

        assert_eq!(counter.get_count(), 0);
        assert_eq!(counter.to_string(), "".to_string());

        for i in 1..11 {
            counter.add_match(KmerMatch {
                id_read: i,
                position: i,
                forward: (i % 2) != 0,
            });

            assert_eq!(counter.get_count(), i);
        }

        assert_eq!(counter.get_count(), 10);
        assert_eq!(counter.to_string(), "(1,1,true) (2,2,false) (3,3,true) (4,4,false) (5,5,true) (6,6,false) (7,7,true) (8,8,false) (9,9,true) (10,10,false) ".to_string());
    }
}
