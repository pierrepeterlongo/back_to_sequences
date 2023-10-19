# Kmer2sequences
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Description
Given a set of kmers (fasta / fastq [.gz] format) and a set of sequences  (fasta / fastq [.gz] format), this tool will extract the sequences containing the kmers.
Each sequence that shares at least a kmer with the indexed kmeres is output with its original header + the number of shared kmers + the ratio of shared kmers:
```
>original_header 20 6.13
TGGATAAAAAGGCTGACGAAAGGTCTAGCTAAAATTGTCAGGTGCTCTCAGATAAAGCAGTAAGCGAGTTGGTGTTCGCTGAGCGTCGACTAGGCAACGTTAAAGCTATTTTAGGC...
```
In this case 20 kmers are shared with the indexed kmers. This represents 6.13% of the kmers in the sequence.

## Install:
```bash
git clone https://github.com/pierrepeterlongo/kmer2sequences.git
cargo install --path . --locked
```

## Quick benchmark
Reproducible by running `bench.sh`` in the benchs folder. Presented results were obtained on GenOuest platform on a node with 32 cores (128 threads) Xeon 2.2 GHz.

Indexed: one million kmers of length 31 (takes 2s). The index size is 29MB. Used 32 threads
Queried: from 1 to 100 million reads, each of average length 350
Reproducible by running `benchmark.sh` in the scripts folder.
Results obtained on a macbook pro Apple M2 Pro, 16Go RAM 
* Indexed: 100,000 kmers of length 31
* Queried: from 1 to 1 million reads, each of average length 500

| Number of reads | Time  |  max RAM |
|-----------------|----------|---|
| 10,000          | 1s    	 | 7 GB |
| 100,000         | 5s    	 | 7 GB |
| 1,000,000       | 15.0   	 | 7 GB |

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
* [X] Add a validation test (04/10/2023)
* [X] Add a number of shared kmers per sequence instead of only their ratio 
	* [ ] ? add a threshold on the number of shared kmers
* [ ] Parallelize the read extraction step
* [ ] Thinks about a way to adapt this to protein sequences
* [X] Add an option to set the size of the bloom filter used by kmindex
* [ ] Estimate the FP rate (that should be null or negligible when the number of kmer to search is lower than a few thousands)