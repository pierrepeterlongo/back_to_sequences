# Kmer2sequences

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


<img src="k2s.jpg" alt="Old library - IA generated" width="150" height="150">


<!-- ![Old library - IA generated](k2s.jpg) -->
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
git clone https://github.com/pierrepeterlongo/back_to_sequences.git
cd back_to_sequences
cargo install --path . --locked
```
 
## Quick benchmark
Reproducible by running `bench.sh` in the benchs folder. 
Presented results were obtained on 
* the GenOuest platform on a node with 32 cores (128 threads) Xeon 2.2 GHz, denoted by "genouest" in the table below.
* and a macbook, Apple M2 pro, 16 GB RAM, denoted by "mac" in the table below.

We indexed: one million kmers (exactly 995,318) of length 31.

We queried: from 10,000 to 100 million reads, each of average length 350.

| Number of reads | Time genouest | Time mac |  max RAM |
|-----------------|----------|---|---|
| 10,000          | 2s   | 	3s | 7 GB |
| 100,000         | 4s   | 	4s | 7 GB |
| 1,000,000       | 12s  | 15s	 | 7 GB |
| 10,000,000       | 1m10  | 1m32	 | 7 GB |
| 100,000,000       | 10m48 | 23m41	 | 7 GB |

## Usage
### help
```	
Extract sequences that contain some kmers

Usage: back_to_sequences [OPTIONS] --in-sequences <IN_SEQUENCES> --in-kmers <IN_KMERS> --out-sequences <OUT_SEQUENCES>

Options:
      --in-sequences <IN_SEQUENCES>    Input fasta or fastq [.gz] file containing the original sequences (eg. reads)
      --in-kmers <IN_KMERS>            Input fasta file containing the original kmers
      --out-sequences <OUT_SEQUENCES>  Output file containing the filtered original sequences (eg. reads). It will be automatically in fasta or fastq format depending on the input file
      --out-kmers <OUT_KMERS>          If provided, output text file containing the kmers that occur in the reads with their number of occurrences [default: ]
  -k, --kmer-size <KMER_SIZE>          Size of the kmers to index and search [default: 31]
      --threshold <THRESHOLD>          Threshold of the ratio of kmers that must be found in a sequence to keep it (default 0). Thus by default, if no kmer is found in a sequence, it is not output [default: 0]
      --stranded                       Used original kmer strand (else canonical kmers are considered)
  -h, --help                           Print help
  -V, --version                        Print version
```

### Example 
```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta  --out-kmers counted_kmers.txt
```

The `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `compacted_kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the estimated ratio of shared kmers and number of shared kmers.

If the `--out-kmers` option is used, the file `counted_kmers.txt` contains for each kmer in `compacted_kmers.fasta` the number of times it was found in `filtered_reads.fasta` (displays only kmers whose counts are higher than 0).

### Generate random data for testing
```bash
# Generate 1 reference sequence of random length 50000 and minimum length 100
python scripts/generate_random_fasta.py 1 50000 100 ref_seq.fasta

# Extract 1000 random "reads", each of length in [100;500] from the reference sequence
python3 scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 100 --max_size 500 --num 1000 --output reads.fasta 

# From those reads, extract 500 random sequence containing the kmers. Those kmers are stored in sequences of length in [31;70]
python3 scripts/extract_random_sequences.py --input reads.fasta --min_size 31 --max_size 70 --num 500 --output compacted_kmers.fasta
```



# TODO
* [X] Add a validation test (04/10/2023)
* [X] Add a number of shared kmers per sequence instead of only their ratio 
* [ ] ? add a threshold on the number of shared kmers
* [X] Parallelize the read extraction step
* [ ] Thinks about a way to adapt this to protein sequences
* [X] Add an option to set the size of the bloom filter used by kmindex
