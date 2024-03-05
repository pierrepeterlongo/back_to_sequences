# Back to sequences

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


<img src="k2s.jpg" alt="Old library - IA generated" width="150" height="150">


<!-- ![Old library - IA generated](k2s.jpg) -->
## Description

Given a set $K$ of kmers (fasta / fastq [.gz] format) and a set of sequences  (fasta / fastq [.gz] format), this tool will extract the sequences containing some of those kmers.

A minimal ($m$) and a maximal ($M$) thresholds are proposed. A sequence whose percentage of kmers shared with $K$ are in $]m, M]$ is output with its original header + the number of shared kmers + the ratio of shared kmers:
```
>original_header 20 6.13
TGGATAAAAAGGCTGACGAAAGGTCTAGCTAAAATTGTCAGGTGCTCTCAGATAAAGCAGTAAGCGAGTTGGTGTTCGCTGAGCGTCGACTAGGCAACGTTAAAGCTATTTTAGGC...
```
In this case 20 kmers are shared with the indexed kmers. This represents 6.13% of the kmers in the sequence.

## Citation

> Anthony Baire, Pierre Peterlongo
[Back to sequences: find the origin of kmers](https://doi.org/10.1101/2023.10.26.564040). bioRxiv 
2023.10.26.564040; doi: https://doi.org/10.1101/2023.10.26.564040

```bibtex
@article{baire_back_2023,
	title = {Back to sequences: find the origin of kmers},
	url = {https://www.biorxiv.org/content/early/2023/10/29/2023.10.26.564040},
	doi = {10.1101/2023.10.26.564040},
	journaltitle = {{bioRxiv}},
	author = {Baire, Anthony and Peterlongo, Pierre},
	date = {2023},
}
```

## Install

```bash
git clone https://github.com/pierrepeterlongo/back_to_sequences.git
cd back_to_sequences
cargo install --path . 
```

A test can be performed by running `cd tiny_test; sh tiny_test.sh; cd -`.
 
## Quick benchmark (obtained with version v0.2.7)
This benchmark is reproducible by running `generate_data.sh` and then `bench.sh` in the `benchs` folder. 
Presented results were obtained on 
* the GenOuest platform on a node with 32 threads Xeon 2.2 GHz, denoted by "genouest" in the table below.
* and a MacBook, Apple M2 pro, 16 GB RAM, with 10 threads denoted by "mac" in the table below.

We indexed: one million kmers (exactly 995,318) of length 31.

We queried: from 10,000 reads to 200 million reads (+ 1 billion on the cluster), each of length 100.

| Number of reads | Time genouest | Time mac |  max RAM |
|-----------------|----------|---|---|
| 10,000          | 0.6s  | 	0.5s | 0.13 GB |
| 100,000         | 1.2s  | 	0.8s | 0.13 GB |
| 1,000,000       | 2.9s  | 3.5s	 | 0.13 GB |
| 10,000,000      | 9.0s  | 11.2s	 | 0.13 GB |
| 100,000,000     | 46.6s | 57.4	 | 0.13 GB |
| 200,000,000     | 1m24  | 1m47     | 0.13 GB |
| 1 billion       | 7m11  | -     | 0.13 GB |

## Usage
### Help
```	
Back to sequences: find the origin of kmers

Usage: back_to_sequences [OPTIONS] --in-kmers <IN_KMERS>

Options:
      --in-sequences <IN_SEQUENCES>
          Input fasta or fastq [.gz] file containing the original sequences (eg. reads). 
              The stdin is used if not provided 
              (and if `--in_filelist` is not provided neither) [default: ]
      --in-filelist <IN_FILELIST>
          Input txt file containing in each line a path to a fasta or fastq [.gz] file 
          containing the original sequences (eg. reads). 
              Note1: if this option is used, the `--out_filelist` option must be used.
                     The number of lines in out_filelist must be the same as in_filelist
              Note2: Incompatible with `--in_sequences` [default: ]
      --in-kmers <IN_KMERS>
          Input fasta file containing the original kmers
      --out-sequences <OUT_SEQUENCES>
          Output file containing the filtered original sequences (eg. reads). 
          It will be automatically in fasta or fastq format depending on the input file.
          If not provided, only the in_kmers with their count is output [default: ]
      --out-filelist <OUT_FILELIST>
          Output txt file containing in each line a path to a fasta or fastq [.gz] file 
          that will contain the related output file from the input files list  [default: ]
      --out-kmers <OUT_KMERS>
          If provided, output a text file containing the kmers that occur in the reads 
          with their number of occurrences
              Note: if `--in_filelist` is used the output counted kmers are 
              those occurring the last input file of that list [default: ]
      --counted-kmer-threshold <COUNTED_KMER_THRESHOLD>
          If out_kmers is provided, output only reference kmers whose number of occurrences 
          is at least equal to this value.
          If out_kmers is not provided, this option is ignored [default: 0]
  -k, --kmer-size <KMER_SIZE>
          Size of the kmers to index and search [default: 31]
  -m, --min-threshold <MIN_THRESHOLD>
          Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
          Minimal threshold of the ratio  (%) of kmers that must be found in a sequence to keep it (default 0%).
          Thus by default, if no kmer is found in a sequence, it is not output. [default: 0]
      --max-threshold <MAX_THRESHOLD>
          Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
          Maximal threshold of the ratio (%) of kmers that must be found in a sequence to keep it (default 100%).
          Thus by default, there is no limitation on the maximal number of kmers found in a sequence. [default: 100]
      --stranded
          Used original kmer strand (else canonical kmers are considered)
      --query-reverse
          Query the reverse complement of reads. Useless without the --stranded option
      --no-low-complexity
          Do not index low complexity kmers (ie. with a Shannon entropy < 1.0)
  -h, --help
          Print help
  -V, --version
          Print version
```

### Examples
#### Basic 
```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta  --out-kmers counted_kmers.txt
```

The `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `compacted_kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the estimated ratio of shared kmers and number of shared kmers.

If the `--out-kmers` option is used, the file `counted_kmers.txt` contains for each kmer in `compacted_kmers.fasta` the number of times it was found in `filtered_reads.fasta` (displays only kmers whose counts are higher than 0).

#### Using filters
```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta  --out-kmers counted_kmers.txt --min-threshold 50 --max-threshold 70
```

In this case only sequeces from `reads.fasta` that have more than 50% and at most 70% of their kmers in `compacted_kmers.fasta` are output.

#### Specifying strands

```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta --stranded
```
In this case, the kmers found in `compacted_kmers.fasta` are indexed in their original orientation, and kmers extracted from `reads.fasta` are queried in their original orientation. 

Note that without the `--stranded` option, all kmers (indexed and queried) are considered in their canonical form.


One may be interested in finding kmers from the reverse complement of the queried sequences. In this case we add the `--query-reverse` option together with the `--stranded` option:
```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta --stranded
```

#### Reading sequences from standard input: 

```bash
cat reads.fasta | back_to_sequences --in-kmers compacted_kmers.fasta --out-sequences filtered_reads.fasta 
```
Do not provide the `--in-sequences` if your input data are read from stdin.

#### Using several input read sets. 
Say you have three input read files `1.fa`, `2.fa`, `3.fa` to which you wish to apply `back_to_sequences`. 

1. create an input and an output files:
```bash
ls 1.fa 2.fa 3.fa > in.fof
echo 1_out.fa 2_out.fa 3_out.fa > out.fof
```
2. run `back_to_sequences`
```bash
back_to_sequences --in-filelist in.fof --in-kmers compacted_kmers.fasta --out-filelist out.fof 
```

## Generate random data for testing
You may be interested by generating a specific data set.
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
* [ ] Provide a way to index and query more than one set $K$ of kmers
* [ ] Output the strand of matched kmers
