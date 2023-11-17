* 0.1.0. 15/09/2023: intial version
* 0.2.0. 19/10/2023: parallelized version 
* 0.2.1. 20/10/2023: changed the option names
* 0.2.2. 20/10/2023: added the `--query_reverse` option
* 0.2.3. 24/10/2023: accepted various optimizations from Anthony Baire
* 0.2.4. 25/10/2023: added the max-threshold option
* 0.2.5. 26/10/2023: added the possibility to read input sequences from stdin
* 0.2.6. 27/10/2023: 
    * implement uppercase and base_complement operations with lookup tables
    * refactor the pipeline
        * move the output file writer into a separate thread
        * run the reader & writer threads outside the rayon threadpool (because they are not cpu-bound)
        * split the work in chunks of 32 reads (this should greatly improve the throughput on machines with many cpus)
    * Some cleaning in the error handling (remove silently ignored errors, propagate errors and avoid using panic!)
* 0.2.7. 31/10/2023: New optimizations
    * optimisation: use a faster hash function (ahash)
    * optimisation: remove redundant hashmap lookup
    * optimisation: store the canonical kmer into a fixed-size slice
* 0.2.8. 17/11/2023: 
    * Deal with lower case letters in the input sequences (all is converted to upper)
    * no not consider kmers that contain non ACGT letters (indexing and querying)
    