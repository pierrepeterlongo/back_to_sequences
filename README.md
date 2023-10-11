# Kmer2sequences

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Description

Given a set of kmers (fasta / fastq [.gz] format) and a set of sequences  (fasta / fastq [.gz] format), this tool will extract the sequences containing the kmers.
Each sequence is output with its original header + the ratio of shared kmers + the number of shared kmers.:

```
>original_header 0.011494253 1
GAATTGTGATACTTGCTGCCGTTAACACAGCCACTCACCTCTGTACACCACACTGGTCCGTGGAGGGTGACAAGCATAACATAGTTCGTATGTGTTGCACGCCCT
```

**key idea**: 

1. kmers (even if they are few) are indexed using kmindex. 

2. Then reads are queried against the index, again using kmindex.

3. Finally, the reads are extracted from the fasta file, given the kmindex output (that contains only headers)
   It is possible to set up a threshold on the ration of shared kmers between the read and the query.
   
   There are two commands to use: 

4. `back_to_sequences index_kmers`: indexes the kmers

5. `back_to_sequences query_sequences`: find the sequences that contain the indexed kmers
   A full example with options is given below.

## Install:

```bash
git clone https://github.com/pierrepeterlongo/kmer2sequences.git
cd kmer2sequences
cargo install --path . --locked
```

## Dependencies

* [kmindex](https://github.com/tlemane/kmindex)
  * kmindex also requires [kmtricks](https://github.com/tlemane/kmtricks). kmtricks is provided by the conda version of kmindex, but if you install kmindex from source, you will need to install kmtricks.

For compiling with mac, cf the note at the end of the file.

## Quick benchmark

Reproducible by running `bench.sh` in the `benchs` folder.
Presented results were obtained on GenOuest platform on a node with 64 cores (128 threads) Xeon 2.2 GHz (L1 = 48KB, L2 = 1.25MB, L3 = 48MB shared) with 900 GB of memory. 

Non presented results on a macbook pro Apple M2 Pro, 32Go RAM, are almost identical.

* Indexed: one million kmers of length 31 (takes 2s). The index size is 29MB. 
* Queried: from 1 to 100 million reads, each of average length 350
  * Filtering time corresponds to the `query_sequences` subcommand.
  * Exact_count Time corresponds to the `exact_count` sub command 
* `Nb filtered reads` are the number of reads containaing at least one kmer after the filtering step. 

| Number of reads | Filtering Time | Exact_count Time | Nb filtered reads | max RAM (GB) |
| --------------- | -------------- | ---------------- | ----------------- | ------------ |
| 1               | 0m1s       | 0m2s         | 0                 | 0.12         |
| 10              | 0m1s       | 0m2s         | 1                 | 0.12         |
| 100             | 0m1s       | 0m2s         | 17                | 0.12         |
| 1,000           | 0m1s       | 0m2s         | 169               | 0.12         |
| 10,000          | 0m1s       | 0m2s         | 1,596             | 0.12         |
| 100,000         | 0m2s       | 0m5s         | 16,007            | 0.53         |
| 1,000,000       | 0m20s      | 0m27s        | 160,589           | 6.8         |
| 10,000,000      | 1m45s       | 3m34s        | 1,601,178         | 21        |
| 100,000,000     | 12m53s     | 35m31s            | 16,032,079        | 21          |

## complete example:

### For testing: generate random reads and extract some of their kmers:

```bash
# Generate 1 reference sequence of random length 50000 and minimum length 100
python scripts/generate_random_fasta.py 1 50000 100 ref_seq.fasta

# Extract 1000 random "reads", each of length in [100;500] from the reference sequence
python3 scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 100 --max_size 500 --num 1000 --output reads.fasta 


# From those reads, extract 500 random sequence containing the kmers. Those kmers are stored in sequences of length in [31;70]
python3 scripts/extract_random_sequences.py --input reads.fasta --min_size 31 --max_size 70 --num 500 --output compacted_kmers.fasta

# Create the file of file, used by kmindex, containing the kmers. D is simply a prefix. 
echo ref_set:compacted_kmers.fasta > fof.txt
```

### Index kmers and extract the reads containing the kmers:

1. index the kmers: 
   
   ```bash
   back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ./bin/kmindex
   ```

2. extract the reads containing the kmers: 
   
   ```bash
   back_to_sequences query_sequences --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --kmindex_path ./bin/kmindex
   ```

**Note:** This second step can be subdivided into two commands (for debugging purpose, or for obtaining only headers of sequences containing the kmers):

    2.1 search the kmers in the reads (finds the headers of reads continaing indexed kmers): 

```bash
back_to_sequences query_sequences_get_headers --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_headers headers --kmindex_path ./bin/kmindex
```

    2.2 back to the read sequences

```bash
back_to_sequences query_sequences_to_reads --in_headers headers --in_fasta reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --threshold 0.0
```

That's all, the `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `compacted_kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the ratio of shared kmers.

## Validation

The kmindex approach is probabilistic, as it uses a bloom filter for indexing the kmers.
We propose a way to validate the results obtained by `kmer2sequences`, by using a second approach, based on a hash table. It requires more time and memory than `kmer2sequences`. It can be applied on the results of the filtered reads (as `kmer2sequences` does not suffer from false negative calls). 

Here is a small example. Suppose you estimated the number of shared kmers using the commands above: 

```bash
back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ./bin/kmindex
back_to_sequences query_sequences --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --kmindex_path ./bin/kmindex
```

You may now verify the exactness or the overestimations using: 

```bash
back_to_sequences exact_count --in_kmers compacted_kmers.fasta --in_fasta filtered_reads.fasta --out_fasta filtered_reads_exact.fasta -k 31 --out_counted_kmers counted_kmers.txt
```

The `filtered_reads_exact.fasta` file is the same as `filtered_reads.fasta`, except that each head contains an additional integer value, being the exact number of shared kmers with the original `compacted_kmers.fasta` file.
If the `--out_counted_kmers` option is used, the file `counted_kmers.txt` contains for each kmer in `compacted_kmers.fasta` the number of times it was found in `filtered_reads.fasta` (displays only kmers whose counts are higher than 0).

Example of a header in `filtered_reads_exact.fasta`:

```
>sequence453 0.003968254 1 1
CGGTTCGAGGCTGGCCTGAGCCACCGTGCCTAA...
```

The first two values (0.003968254 1) are those already computed by kmer2sequences. The third value (1) is the exact number of shared kmers. In this case the estimation is perfect. 

Example of a line in `counted_kmers.txt`:

```
CGTCATTTCCTGGGTCACAGTGAACGGACCC 1
```

Kmer `CGTCATTTCCTGGGTCACAGTGAACGGACCC` occurs once in `filtered_reads.fasta`.

## kmindex and kmtricks compilation note for mac os users.

kmindex install is a bit complex (sept 2023)

* kmtricks: does not compile on mac since 1.2.1. Solution: 
  * Installation via conda.
  * In case of conda issues between noarch and arm64. Solution: 
    
    ```bash
    conda config --env --set subdir osx-64
    conda create -p env_kmtricks_1.2.1 kmtricks=1.2.1 
    conda activate ./env_kmtricks_1.2.1
    ```
* kmindex: commit 1b018539a4a1730a51840ed5a9330c023baf3814
  * For compiling with a mac: change lines ` asm volatile("pause");` by `asm volatile("isb" : : : "memory");` in `lib/include/kmindex/spinlock.hpp`. 
  * Comment `if (!(kmv >= min_kmv_required))` line 219 of `app/kmindex/build.cpp`
    (sorry for the trouble)

# TODO

* [x] Add a validation test (04/10/2023)
* [x] Add a number of shared kmers per sequence instead of only their ratio 
  * [ ] ? add a threshold on the number of shared kmers
* [ ] Parallelize the read extraction step
* [ ] Thinks about a way to adapt this to protein sequences
* [x] Add an option to set the size of the bloom filter used by kmindex
* [ ] Estimate the FP rate (that should be null or negligible when the number of kmer to search is lower than a few thousands)