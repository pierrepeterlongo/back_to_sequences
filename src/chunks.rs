#![allow(missing_docs, private_bounds)]

use std::marker::PhantomData;
use std::ops::Range;
use std::sync::mpsc::{SendError, sync_channel};

use ahash::AHashMap as HashMap;
use needletail::{FastxReader, Sequence};
use rayon::prelude::*;

const CHUNK_BUF_SIZE: usize = 65536; // total buffer size
const CHUNK_RECORDS_SIZE: usize = 512; // initial size of the record vector
const INPUT_CHANNEL_SIZE: usize = 8;
const OUTPUT_CHANNEL_SIZE: usize = 8;

enum Never{}

/// Chunk of fastx records
pub struct Chunk<X, O = WithId>
where X: Send,
      O: Send,
{
    chunk_id: usize,
    buf: Vec<u8>,
    records: Vec<InnerRecord<X, O>>,
}


impl<X, O> Chunk<X, O>
where X: Send,
      O: Send + OptionalId,
{
    pub fn try_for_each<F, E>(&mut self, func: &mut F) -> Result<(), E>
        where F: FnMut(Record<'_, X>) -> Result<(), E>
    {
        for rec in &mut self.records {
            let [id, seq] = self.buf.get_disjoint_mut([rec.id_start.get(&rec.seq), rec.seq.clone()]).unwrap();
            func(Record{
                read_id: rec.read_id,
                id, seq,
                extra: &mut rec.extra,
            })?;
        }
        Ok(())
    }

    pub fn for_each<F>(&mut self, func: &mut F)
        where F: FnMut(Record<'_, X>)
    {
        let _ = self.try_for_each(&mut |r| -> Result<_, Never> { func(r); Ok(()) });
    }
}

/// structure representing a fastx record
///
/// This is the type presented in the API.
pub struct Record<'a, X> {
    /// index id of this record (as it appears in the input)
    pub read_id: usize,

    /// id of the fastx record
    pub id: &'a [u8],

    /// normalized fastx sequence
    pub seq: &'a mut [u8],

    /// extra field containing user data
    pub extra: &'a mut X,
}

/// internal structure representing a fastx record
struct InnerRecord<X, O> {
    /// index id of this record (as it appears in the input)
    read_id: usize,

    /// position of the record id string in the chunk buffer (optional)
    id_start: O,

    /// position of the sequence in the chunk buffer
    seq: std::ops::Range<usize>,

    /// extra field containing user data
    extra: X
}

/// Trait for including or excluding the "id" field in the fastx record
trait OptionalId: Sized {
    fn len(id: &[u8]) -> usize;
    fn set(id: &[u8], buf: &mut Vec<u8>) -> Self;
    fn get(&self, seq: &Range<usize>) -> Range<usize>;
}

/// Marker to skip the fastx record id when building the chunks
///
/// As a consequence [Record::id] will always be empty
pub struct WithoutId;

/// Marker to include a copy of the fastx record id when building the chunks
pub struct WithId(usize);

impl OptionalId for WithoutId
{
    fn len(_id: &[u8]) -> usize { 0 }
    fn set(_id: &[u8], _buf: &mut Vec<u8>) -> Self { Self }
    fn get(&self, _seq: &Range<usize>) -> Range<usize> { 0..0 }
}

impl OptionalId for WithId
{
    fn len(id: &[u8]) -> usize {
        id.len()
    }
    fn set(id: &[u8], buf: &mut Vec<u8>) -> Self {
        let start = buf.len();
        buf.extend_from_slice(id);
        Self(start)
    }
    fn get(&self, seq: &Range<usize>) -> Range<usize> {
        self.0..seq.start
    }
}

/// Iterator for converting a fastx input into a sequence of [Chunk]
pub struct ChunksReader<X, O>
where X: Send,
      O: Send,
{
    reader: Box<dyn FastxReader>,
    chunk: Chunk<X, O>,
    chunk_id: usize,
    read_id: usize,
}

impl<X, O> ChunksReader<X, O>
where X: Send,
      O: Send,
{
    /// Create a new chunk
    fn new_chunk(chunk_id: usize) -> Chunk<X, O> {
        Chunk{
            chunk_id,
            buf:     Vec::with_capacity(CHUNK_BUF_SIZE),
            records: Vec::with_capacity(CHUNK_RECORDS_SIZE),
        }
    }
}

impl<X, O> Iterator for ChunksReader<X, O>
where X: Send + Default,
      O: Send + OptionalId,
{
    type Item = anyhow::Result<Chunk<X, O>>;

    fn next(&mut self) -> Option<Self::Item> {
        // println!("chunk state: #{} ({} records, {} bytes)",
        //          self.chunk.chunk_id, self.chunk.records.len(), self.chunk.buf.len());

        let new_chunk = loop {
            let read_id = self.read_id;
            self.read_id += 1;
            let seq_record = match self.reader.next() {
                None if self.chunk.records.is_empty() => return None,
                None         => break Self::new_chunk(self.chunk_id),
                Some(Err(e)) => return Some(Err(e.into())),
                Some(Ok(r))  => r,
            };

            // create the new chunk record
            let id  = seq_record.id();
            let seq = seq_record.normalize(false);

            let push_record = |chunk: &mut Chunk<X, O>| {
                let id_start  = O::set(id, &mut chunk.buf);
                let seq_start = chunk.buf.len();
                chunk.buf.extend_from_slice(&seq);
                let end       = chunk.buf.len();
                chunk.records.push(InnerRecord{
                    read_id, id_start,
                    seq: seq_start..end,
                    extra: X::default(),
                });
            };

            // flush the current chunk
            // if the new record would overfill its buffer
            let needed_capacity = self.chunk.buf.len() + O::len(id) + seq.len();
            if self.chunk.buf.capacity() < needed_capacity && !self.chunk.records.is_empty()
            {
                let mut new_chunk = Self::new_chunk(self.chunk_id);
                push_record(&mut new_chunk);
                break new_chunk;
            }

            push_record(&mut self.chunk);
        };
        // println!("new chunk: #{} ({} records, {} bytes)",
        //          self.chunk.chunk_id, self.chunk.records.len(), self.chunk.buf.len());
        self.chunk_id = new_chunk.chunk_id + 1;
        Some(Ok(std::mem::replace(&mut self.chunk, new_chunk)))
    }
}

/// Build an iterator for converting a fastx input into a sequence of [Chunk] values
pub fn from_fastx_reader<X, O>(reader: Box<dyn FastxReader>) -> anyhow::Result<ChunksReader<X, O>>
where X: Send,
      O: Send,
{
    Ok(ChunksReader{
        reader,
        chunk: ChunksReader::new_chunk(0),
        chunk_id: 1,
        read_id: 0,
    })
}

/// Re-assemble a sequence of [Chunk] values and call a faillible closure on each record
pub fn try_for_each<I, X, O, F, E>(iter: I, func: &mut F) -> Result<(), E>
    where I: Iterator<Item=Chunk<X, O>>,
          X: Send,
          O: Send + OptionalId,
          F: FnMut(Record<'_, X>) -> Result<(), E>

{
    for mut chunk in iter {
        chunk.try_for_each(func)?;
    }
    Ok(())
}

/// Re-assemble a sequence of [Chunk] values and call a closure on each record
pub fn for_each<I, X, O, F>(iter: I, func: &mut F)
    where I: Iterator<Item=Chunk<X, O>>,
          X: Send,
          O: Send + OptionalId,
          F: FnMut(Record<'_, X>)
{
    for mut chunk in iter {
        chunk.for_each(func);
    }
}

/// Re-assemble and re-order a sequence of [Chunk] values and call a faillible closure on each
/// record
pub fn try_for_each_reordered<I, X, O, F, E>(mut iter: I, func: &mut F) -> Result<(), E>
    where I: Iterator<Item=Chunk<X, O>>,
          X: Send,
          O: Send + OptionalId,
          F: FnMut(Record<'_, X>) -> Result<(), E>
{
    let mut buffer = HashMap::new();

    for chunk_id in 0.. {
        let mut chunk = match buffer.remove(&chunk_id) {
            Some(chunk) => chunk,
            None => loop {
                match iter.next() {
                    None => {
                        // NOTE: non-empty buffer may happen if there is an i/o error in the input thread
                        debug_assert!(buffer.is_empty());
                        return Ok(());
                    },
                    Some(chunk) if chunk.chunk_id == chunk_id => break chunk,
                    Some(chunk) => { buffer.insert(chunk.chunk_id, chunk) ;},
                }
            }
        };
        chunk.try_for_each(func)?;
    }
    Ok(())
}

/// Re-assemble and re-order a sequence of [Chunk] values and call a closure on each
/// record
pub fn for_each_reordered<I, X, O, F, E>(iter: I, func: &mut F)
    where I: Iterator<Item=Chunk<X, O>>,
          X: Send,
          O: Send + OptionalId,
          F: FnMut(Record<'_, X>),
{
    let _ = try_for_each_reordered(iter, &mut |rec| -> Result<_, Never> { func(rec); Ok(()) });
}


fn default_id() {}
fn default_op(_:(), _:()) {}

/// Type alias for the reduce function pair
pub type NoReduceType = (fn(), fn((), ()));

/// no-op `reduce` function for [Pipeline::run()]
pub const NO_REDUCE: NoReduceType = (default_id, default_op);


/// no-op `writer` function for [Pipeline::run()]
#[allow(non_snake_case)]
pub fn NO_WRITER<X>(_: Record<'_, X>) -> anyhow::Result<()> {
    Ok(())
}



/// Pipeline for processing fastx sequences
///
/// see [Self::run()]
pub struct Pipeline<X=(), O=WithId> (PhantomData<(X, O)>);

impl<X, O> Pipeline<X, O>
where X: Send + Default,
      O: Send + OptionalId
{
    /// Process a series of [needletail] fastx records over a [rayon]-parallelised map-reduce pipeline
    ///
    /// - `reader` is the [needletail::FastxReader] object providing the input sequences
    /// - `map` is a closure that will be called on each record
    ///   - it may mutate [Record::seq] (but not change its size) and [Record::extra]
    ///   - it may return a value to be processed by the `reduce` argument
    ///   - its is run in a [rayon] thread and may be used to perform cpu-intensive operations
    ///     (note: the records grouped in [Chunk]s of 64Ko in order to improve the throughput)
    /// - `reduce` is an optional pair of closures for reducing the results produced with the `map`
    ///   closure to a single value
    ///   - see [rayon::iter::ParallelIterator::reduce()]
    ///   - use [NO_REDUCE] to disable this feature
    /// - `writer` is an optional closure called on each record at the end of the pipeline
    ///   - use [NO_WRITER] ot disable this feature
    ///   - it is run in a single dedicated thread and should not perform cpu-intensive operations
    ///
    /// - if the [Record::id] is not needed then the `O` type parameter may be set to [WithoutId]
    ///   to skip it when building the chunks and speed-up the process
    ///
    pub fn run<M, W,  ID, OP, R>(
        reader: Box<dyn FastxReader>, map: M, reduce: (ID, OP), mut writer: W) -> anyhow::Result<R>

    where M: Fn(Record<'_, X>) -> R + Sync,
          W: FnMut(Record<'_, X>) -> anyhow::Result<()>,
          ID: Fn() -> R + Sync + Send,
          OP: Fn(R, R) -> R + Sync + Send,
          R: Send + Copy,
      {
          let (input_tx, input_rx) = sync_channel::<Chunk<X, O>>(INPUT_CHANNEL_SIZE);
          let (output_tx, output_rx) = sync_channel::<Chunk<X, O>>(OUTPUT_CHANNEL_SIZE);

          let (reduce_id, reduce_op) = (&reduce.0, &reduce.1);

          std::thread::scope(|s| -> anyhow::Result<R> {

              let reader_thread = s.spawn(move || -> anyhow::Result<()> {
                  for chunk in from_fastx_reader(reader)? {
                      if input_tx.send(chunk?).is_err() {
                          // abort on send error but do not report any error
                          // (at this point the receiver thread would yield a more meaningful error)
                          break;
                      }
                  }
                  Ok(())
              });

              let map_reduce_thread = s.spawn(|| -> anyhow::Result<R> {
                  let output_tx = output_tx;

                  input_rx.into_iter().par_bridge().map(|mut chunk| {
                      let mut result = reduce_id();
                      chunk.for_each(&mut |record| {
                          result = reduce_op(result, map(record));
                      });
                      output_tx.send(chunk).map_err(|_| SendError(()))?;
                      Ok(result)
                  })
                  .try_reduce(reduce_id, |a, b| Ok(reduce_op(a, b)))
              });

              try_for_each_reordered(output_rx.into_iter(), &mut writer)?;

              reader_thread.join().expect("panic in the reader thread")?;
              map_reduce_thread.join().expect("panic in the map reduce thread")
          })
      }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use biotest::Format;
    use needletail::Sequence;
    use super::*;

    #[test]
    fn pipeline() -> anyhow::Result<()>
    {
        let seq_len = 150;

        let mut rng = biotest::rand();
        let s_generator = biotest::Fasta::builder().sequence_len(seq_len).build()?;

        let mut fasta_sequences = Vec::<u8>::with_capacity(512*1024);
        s_generator.records(&mut fasta_sequences, &mut rng, 1000)?;
        let make_reader = || needletail::parse_fastx_reader(Cursor::new(fasta_sequences.clone())).unwrap();
        let make_channel = || sync_channel::<(String, String, Option<usize>)>(64);

        let to_str = |v: &[u8]| std::str::from_utf8(v).unwrap().to_owned();

        let count_a = |seq: &[u8]| seq.iter().filter(|c| **c==b'A').count();
        let compute_extra = |read_id| Some(read_id*2 + 42);

        let delay = || {
            use std::{sync::LazyLock, thread::{current, ThreadId, sleep}};
            static SLOW_THREAD: LazyLock<ThreadId> = LazyLock::new(|| current().id());
            if current().id() == *SLOW_THREAD {
                sleep(std::time::Duration::from_micros(1));
            }
        };

        std::thread::scope(|s| {
            let (hnd_pipelined, iter_pipelined) = {
                let (snd, rcv) = make_channel();

                (
                    s.spawn(|| -> anyhow::Result<usize> {
                        let snd = snd;
                        Pipeline::<Option<usize>>::run(
                            make_reader(),
                            // map
                            |rec| {
                                // make one thread much slower to verify that the chunks are
                                // correctly reordered
                                delay();

                                // modify the sequence
                                // (replace it with its reverse reverse complement)
                                rec.seq.copy_from_slice(&rec.seq.reverse_complement());

                                // set a value in the extra field
                                *rec.extra = compute_extra(rec.read_id);

                                // return the number of 'A' bases
                                count_a(rec.seq)
                            },
                            // reduce: compute the total sum of 'A' bases
                            (|| 0, |a, b| a + b),

                            // writer: send (id, seq, extra) to the channel
                            |rec| {
                                snd.send((to_str(rec.id), to_str(rec.seq), *rec.extra)).unwrap();
                                Ok(())
                            })
                    }),
                    rcv.into_iter()
                )
            };


            // single-thread implementation (used as a reference)
            let (hnd_reference, iter_reference) = {
                let (snd, rcv) = make_channel();

                (
                    s.spawn(|| -> anyhow::Result<usize> {
                        let snd = snd;
                        let mut reader = make_reader();
                        let mut read_id = 0;
                        let mut total_a = 0;

                        while let Some(rec) = reader.next() {
                            let rec = rec.unwrap();
                            let seq = rec.normalize(false).reverse_complement();
                            total_a += count_a(&seq);
                            let extra = compute_extra(read_id);
                            read_id += 1;
                            snd.send((to_str(rec.id()), to_str(&seq), extra)).unwrap();
                        }
                        Ok(total_a)
                    }),
                    rcv.into_iter()
                )
            };

            // compare each sequence computed by the writer
            let mut seq_id = 0;
            iter_pipelined.zip(iter_reference)
                .for_each(move |(pipelined, reference)| {
                    assert_eq!(pipelined, reference, "at sequence {}", seq_id);
                    seq_id += 1;
                });

            // compare the reduced value
            assert_eq!(
                hnd_pipelined.join().unwrap().unwrap(),
                hnd_reference.join().unwrap().unwrap());
        });

        Ok(())
    }
}
