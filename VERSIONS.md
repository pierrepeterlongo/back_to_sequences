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
* 0.2.9. 17/11/2023: merge branch "cleaning". 
* 0.2.10. 19/11/2023: do not index low complexity kmers (e.g. AAAAAA)
* 0.2.11. 19/11/2023: 
    * add the `--no_low_complexity` option
    * optimize the way kmers containing non ACGT letters are skipped
* 0.2.12. 23/01/2024:
    * possiblity to not output filtered reads. In this case only counted kmers are output. 
* 0.2.13. 23/01/2024:
    * prints all kmers from kmer_set, whaterver their counts count
* 0.3.0 16/02/2024
    * Add the counted_kmer_threshold option
* 0.4.0 23/02/2024
    * Add the possibility to filter several input datasets
* 0.5.0 21/03/2024
    * Added the --output-kmer-positions option
    * code refactoring
* 0.5.1 25/03/2024
    * Clarify the options
* 0.6.1 03/04/2024
    * Added option output-mapping-positions
* 0.6.2 04/04/2024
    * Validate option output-mapping-positions
    * Added tests for output-mapping-positions
* 0.6.3 05/04/2024
    * Added the "thread" option
* 0.6.4 05/04/2024
    * Update benchs
* 0.6.5 15/04/2024
    * Fixed a filtering bug: As percentage of shared kmers were rounded to 2 decimals, and as we exclude the min-threshold value, 
    some reads with shared kmers were discarded while 
    they contained shared kmers. Now the tset is made before to round values + output values are rounded to 5 decimals.
* 0.6.6 26/06/2024
    * Corrected a bug related to the usage of --in-filelist when querying multiple input files
* 0.6.7 4/10/2024
    * b2s compatible with zst compressed data
* 0.6.8 6/10/2024
    * Optimized display function when using `--output-mapping-positions`
    * redirected all log messages to stderr
    * uses fxreads 0.2.14
0.7.0 25/02/2025
    * Using mapping_position also outputs the number of positions covered by at least an input kmer.
    * Outputs the total number of nucleotides in reads
    * Outputs the total number of kmers in reads
    * Outputs the number of matches of input kmers in the reads
0.7.1 05/03/2025 
    * Bug fix for kmer count
0.8.0 20/11/2025
    * Uses needlail. Enables multiline FASTA/FASTQ records