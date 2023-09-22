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
 1. `back_to_sequences index_kmers`: indexes the kmers
 2. `back_to_sequences query_sequences`: find the sequences that contain the indexed kmers
 A full example with options is given below.

## Install:
```bash
git clone https://github.com/pierrepeterlongo/kmer2sequences.git
cargo install --path . --locked
```

## Dependencies
* [kmindex](https://github.com/tlemane/kmindex)
	* kmindex also requires [kmtricks](https://github.com/tlemane/kmtricks). kmtricks is provided by the conda version of kmindex, but if you install kmindex from source, you will need to install kmtricks.

For compiling with mac, cf the note at the end of the file.

## Quick benchmark
Reproducible by running `benchmark.sh` in the scripts folder.
Results obtained on a macbook pro Apple M2 Pro, 16Go RAM 
* Indexed: 100,000 kmers of length 31 (takes 2s)
* Queried: from 1 to 1 million reads, each of average length 500

| Number of reads | Time (s) |  max RAM |
|-----------------|----------|---|
| 1               | 0.2      |5.47 kb |
| 10              | 0.2      | 6.31 kb| 
| 100             | 0.2      |	14.9 kb |
| 1,000           | 0.2      | 37.3 kb |
| 10,000          | 0.2    	 | 101.6 kb |
| 100,000         | 1.2    	 | 0.74 Mb |
| 1,000,000       | 15.0   	 | 6.58 Mb |

## complete example: 
### For testing: generate random reads and extract some of their kmers: 
```bash
# Generate 100000 reads of average length 500 and minimum length 100
python scripts/generate_random_fasta.py 100000 500 100 reads.fasta

# Extract 100 random kmers of length 31 from the reads
python3 scripts/extract_random_kmers_from_a_fasta_file.py --canonical reads.fasta 31 100 kmers.fasta

# Create the file of file, used by kmindex, containing the kmers. D is simply a prefix. 
echo D:kmers.fasta > fof.txt
```

### Index kmers and extract the reads containing the kmers:

1. index the kmers: 
```bash
back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ./bin/kmindex
```

2. extract the headers of the reads containing the kmers: 
```bash
back_to_sequences query_sequences --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --kmindex_path ./bin/kmindex
```

**Note:** This second step can be subdivided into two commands (for debugging purpose, or for obtaining only headers of sequences containing the kmers):

	2.1 search the kmers in the reads: 
```bash
back_to_sequences query_sequences_get_headers --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_headers headers --kmindex_path ./bin/kmindex
```

	2.2 back to the read sequences
```bash
back_to_sequences query_sequences_to_reads --in_headers headers --in_fasta reads.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads.fasta --threshold 0.0
```

That's all, the `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the ratio of shared kmers.


### kmindex and kmtricks compilation note for mac os users.
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
* [ ] Add a validation test
* [X] Add a number of shared kmers per sequence instead of only their ratio 
	* [ ] ? add a threshold on the number of shared kmers
* [ ] Parallelize the read extraction step
* [ ] Thinks about a way to adapt this to protein sequences
* [X] Add an option to set the size of the bloom filter used by kmindex
* [ ] Estimate the FP rate (that should be null or negligible when the number of kmer to search is lower than a few thousnads)
