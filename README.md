# Kmer2sequences

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Description

Given a set of kmers (fasta / fastq [.gz] format) and a set of sequences  (fasta / fastq [.gz] format), this tool will extract the sequences containing the kmers.
<<<<<<< HEAD
Each sequence is output with its original header + the ratio of shared kmers + the number of shared kmers.:

```
>original_header 0.011494253 1 1
GAATTGTGATACTTGCTGCCGTTAACACAGCCACTCACCTCTGTACACCACACTGGTCCGTGGAGGGTGACAAGCATAACATAGTTCGTATGTGTTGCACGCCCT
```

**key idea**: 
`kmer2sequences` works whithin three phases:

1. kmers (even if they are few) are indexed using kmindex. 
2. Then reads are queried against the index, again using kmindex. This step enables to estimate for each read (with a non null false positive probability) the number of kmers it shares with the indexed kmers. The output of this step looks like this: 

```
>39942 0.14563107 15
CGATTTACCCCTACCTATCCCCAACTACAACCGTGAACGTTACTGAACTCCCTTTGGTCTGCATGAAACTGGGGACTTATTGGTGGTATTGAATATTAACTTCTCTCGCATTGATATGGGCACTCCATAGCGA
```
the read id 39942 shares ~14.5% of its kmers (15 kmers) with the indexed kmers.
   
3. If needed: among reads that share at least one kmer with the references, an exact count can be applied to avoid any false positive. 
The output of this step looks like this: 

```
...
>39942 0.14563107 15 15
CGATTTACCCCTACCTATCCCCAACTACAACCGTGAACGTTACTGAACTCCCTTTGGTCTGCATGAAACTGGGGACTTATTGGTGGTATTGAATATTAACTTCTCTCGCATTGATATGGGCACTCCATAGCGA
...
```
the second `15` is the exact number of shared kmers.

If requested by the user the number of occurrences of each shared kmer is provided in a text file like:
```
...
CGCGCTATAGACCGTACGCTCCACCAATTAA 4
AAGCAAGACCTACCGGCTCGTCAAAACAGAA 5
...
```
=======
Each sequence that shares at least a kmer with the indexed kmeres is output with its original header + the number of shared kmers + the ratio of shared kmers:
```
>original_header 20 6.13
TGGATAAAAAGGCTGACGAAAGGTCTAGCTAAAATTGTCAGGTGCTCTCAGATAAAGCAGTAAGCGAGTTGGTGTTCGCTGAGCGTCGACTAGGCAACGTTAAAGCTATTTTAGGC...
```
In this case 20 kmers are shared with the indexed kmers. This represents 6.13% of the kmers in the sequence.
>>>>>>> dev

## Install:

```bash
git clone https://github.com/pierrepeterlongo/kmer2sequences.git
cd kmer2sequences
cargo install --path . --locked
```
<<<<<<< HEAD

## Dependencies

* [kmindex](https://github.com/tlemane/kmindex)
  * kmindex also requires [kmtricks](https://github.com/tlemane/kmtricks). kmtricks is provided by the conda version of kmindex, but if you install kmindex from source, you will need to install kmtricks.

For compiling with mac, cf the note at the end of the file.

## Quick benchmark

Reproducible by running `bench.sh` in the `benchs` folder.
Presented results were obtained on GenOuest platform on a node with 64 cores (128 threads) Xeon 2.2 GHz (L1 = 48KB, L2 = 1.25MB, L3 = 48MB shared) with 900 GB of memory. 

Non presented results on a macbook pro Apple M2 Pro, 32Go RAM, are almost identical.

* Indexed: one million kmers of length 31 (takes 2s). The index size is 29MB. Used 32 threads
* Queried: from 1 to 100 million reads, each of average length 350
  * Filtering time corresponds to the `query_sequences` subcommand. Used 32 threads
  * Exact_count Time corresponds to the `exact_count` sub command Use 1 thread (Parellization of this step is not implemented yet)
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
| 10,000,000      | 1m45s       | 3m34s        | 1,601,178         | 12        |
| 100,000,000     | 12m53s     | 35m31s            | 16,040,388        | 12          |
| 200,000,000   | 23m30s  |  1h12m23s | 160,320,003      | 32,096,008          | 12 |

## complete example:

### For testing: generate random reads and extract some of their kmers:
=======
 
## Quick benchmark
Reproducible by running `bench.sh`` in the benchs folder. Presented results were obtained on GenOuest platform on a node with 32 cores (128 threads) Xeon 2.2 GHz.

* Indexed: one million kmers of length 31 (takes 2s). The index size is 29MB. Used 32 threads
* Queried: from 10,000 to 100 million reads, each of average length 350

| Number of reads | Time  |  max RAM |
|-----------------|----------|---|
| 10,000          | 1s    	 | 7 GB |
| 100,000         | 5s    	 | 7 GB |
| 1,000,000       | 24.0   	 | 7 GB |
| 10,000,000       | 15.0   	 | 7 GB |
>>>>>>> dev

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

<<<<<<< HEAD
1. index the kmers: 
   
```bash
back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ./bin/kmindex
```

2. extract the reads containing at least one indexed kmers: 
   
=======

>>>>>>> dev
```bash
back_to_sequences query_sequences --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --kmindex_path ./bin/kmindex
```

The `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `compacted_kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the estimated ratio of shared kmers and number of shared kmers.

<!-- **Note:** This second step can be subdivided into two commands (for debugging purpose, or for obtaining only headers of sequences containing the kmers):

2.1 search the kmers in the reads (finds the headers of reads continaing indexed kmers): 

```bash
back_to_sequences query_sequences_get_headers --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_headers headers --kmindex_path ./bin/kmindex
```

2.2 back to the read sequences

```bash
back_to_sequences query_sequences_to_reads --in_headers headers --in_fasta reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --threshold 0.0
``` -->

3. Exact count of shared kmers (optional): 

The kmindex approach is probabilistic, as it uses a bloom filter for indexing the kmers.
We may validate previous results, using an exact approach, based on a hash table. It requires more time and memory than the `query_sequences` subcommand. However, it is applied a subset of the original reads, being only those that it is estimated that they share at least one kmer with the indexed kmers.

```bash
back_to_sequences exact_count --in_kmers compacted_kmers.fasta --in_fasta filtered_reads.fasta --out_fasta filtered_reads_exact.fasta -k 31 --out_counted_kmers counted_kmers.txt
```

The `filtered_reads_exact.fasta` file is the same as `filtered_reads.fasta`, except that each head contains an additional integer value, being the exact number of shared kmers with the original `compacted_kmers.fasta` file.
If the `--out_counted_kmers` option is used, the file `counted_kmers.txt` contains for each kmer in `compacted_kmers.fasta` the number of times it was found in `filtered_reads.fasta` (displays only kmers whose counts are higher than 0).

## Summing up
Here is a summary of the commands used in the previous example:

```bash
back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ./bin/kmindex

back_to_sequences query_sequences --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --kmindex_path ./bin/kmindex

back_to_sequences exact_count --in_kmers compacted_kmers.fasta --in_fasta filtered_reads.fasta --out_fasta filtered_reads_exact.fasta -k 31 --out_counted_kmers counted_kmers.txt
```

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
<<<<<<< HEAD

* [x] Add a validation test (04/10/2023)
* [x] Add a number of shared kmers per sequence instead of only their ratio 
  * [ ] ? add a threshold on the number of shared kmers
* [ ] Parallelize the read extraction step
* [ ] Provide a way to index all kmers (not only the canonical form)
* [ ] Thinks about a way to adapt this to protein sequences
* [x] Add an option to set the size of the bloom filter used by kmindex
* [ ] Estimate the FP rate (that should be null or negligible when the number of kmer to search is lower than a few thousands)
=======
* [X] Add a validation test (04/10/2023)
* [X] Add a number of shared kmers per sequence instead of only their ratio 
* [ ] ? add a threshold on the number of shared kmers
* [X] Parallelize the read extraction step
* [ ] Thinks about a way to adapt this to protein sequences
* [X] Add an option to set the size of the bloom filter used by kmindex
>>>>>>> dev
