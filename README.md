# Kmer2sequences

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


![k2.png](k2s.png)
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
cd kmer2sequences
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
| 10,000          | 2s   | 	5s | 7 GB |
| 100,000         | 4s   | 	5s | 7 GB |
| 1,000,000       | 12s  | 23s	 | 7 GB |
| 10,000,000       | 1m17  | 2m59	 | 7 GB |
| 100,000,000       | 12m16 | 23m41	 | 7 GB |

## Usage
### help
```	
Extract sequences that contain some kmers

Usage: back_to_sequences [OPTIONS] --in-fasta-reads <IN_FASTA_READS> --in-fasta-kmers <IN_FASTA_KMERS> --out-fasta-reads <OUT_FASTA_READS>

Options:
      --in-fasta-reads <IN_FASTA_READS>
          Input fasta file containing the original kmers
      --in-fasta-kmers <IN_FASTA_KMERS>
          Input fasta file containing the original kmers
      --out-fasta-reads <OUT_FASTA_READS>
          Output fasta file containing the original kmers
      --out-txt-kmers <OUT_TXT_KMERS>
          If provided, output text file containing the kmers that occur in the reads with their number of occurrences [default: ]
  -k, --kmer-size <KMER_SIZE>
          Number of times to greet [default: 31]
      --stranded
          Used original kmer strand (else canonical kmers are considered)
  -h, --help
          Print help
  -V, --version
          Print version
```

### Example 
```bash
back_to_sequences --in-fasta-kmers compacted_kmers.fasta --in-fasta-reads reads.fasta --out-fasta-reads filtered_reads.fasta  --out-txt-kmers counted_kmers.txt
```

The `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `compacted_kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the estimated ratio of shared kmers and number of shared kmers.

If the `--out-txt-kmers` option is used, the file `counted_kmers.txt` contains for each kmer in `compacted_kmers.fasta` the number of times it was found in `filtered_reads.fasta` (displays only kmers whose counts are higher than 0).

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
