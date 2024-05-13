# Back to sequences

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


<img src="k2s.jpg" alt="Old library - IA generated" width="150" height="150">

<!-- Add a table of content -->
* [Description](#Description)
* [Citation](#Citation)
* [Install](#Install)
* [Quick benchmark (obtained with version v0.6.4)](#Quickbenchmarkobtainedwithversionv0.6.4)
* [Usage](#Usage)
* [Generate random data for testing](#Generaterandomdatafortesting)

<!-- ![Old library - IA generated](k2s.jpg) -->
##  1. <a name='Description'></a>Description

Given a set $K$ of kmers (fasta / fastq [.gz] format) and a set of sequences  (fasta / fastq [.gz] format), this tool will extract the sequences containing some of those kmers.

A minimal ($m$) and a maximal ($M$) thresholds are proposed. A sequence whose percentage of kmers shared with $K$ are in $]m, M]$ is output with its original header + the number of shared kmers + the ratio of shared kmers:
```
>original_header 20 6.13
TGGATAAAAAGGCTGACGAAAGGTCTAGCTAAAATTGTCAGGTGCTCTCAGATAAAGCAGTAAGCGAGTTGGTGTTCGCTGAGCGTCGACTAGGCAACGTTAAAGCTATTTTAGGC...
```
In this case 20 kmers are shared with the indexed kmers. This represents 6.13% of the kmers in the sequence.

##  2. <a name='Citation'></a>Citation

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

##  3. <a name='Install'></a>Install

```bash
git clone https://github.com/pierrepeterlongo/back_to_sequences.git
cd back_to_sequences
RUSTFLAGS="-C target-cpu=native" cargo install --path .
```

A test can be performed by running `cd tiny_test; sh tiny_test.sh; cd -`.
 
##  4. <a name='Quickbenchmarkobtainedwithversionv0.6.4'></a>Quick benchmark (obtained with version v0.6.4)
This benchmark is reproducible by running `generate_data.sh` and then `bench.sh` in the `benchs` folder. 
Presented results were obtained on 
* the GenOuest platform on a node with 32 threads Xeon 2.2 GHz, denoted by "genouest" in the table below.
* and a MacBook, Apple M2 pro, 16 GB RAM, with 10 threads denoted by "mac" in the table below.

We indexed: one million kmers (exactly 995,318) of length 31.

We queried: from 10,000 reads to 200 million reads (+ 1 billion on the cluster), each of length 100.

| Number of reads | Time genouest | Time mac |  max RAM |
|-----------------|----------|---|---|
| 10,000          | 0.4s  | 0.5s | 0.13 GB |
| 100,000         | 0.7s  | 1.0s | 0.13 GB |
| 1,000,000       | 2.8s  | 5.3s | 0.13 GB |
| 10,000,000      | 11.1s  | 14.6s | 0.13 GB |
| 100,000,000     | 1m10s | 1m03s | 0.13 GB |
| 200,000,000     | 2m16s  | 2m03s | 0.13 GB |

##  5. <a name='Usage'></a>Usage
###  5.1. <a name='Help'></a>Help
```	
Back to sequences: find the origin of kmers

Usage: back_to_sequences [OPTIONS] --in-kmers <IN_KMERS>

Options:
      --in-kmers <IN_KMERS>
          Input fasta file containing the original kmers
              Note: back_to_sequences considers the content as a set of kmers
              This means that a kmer is considered only once, 
              even if it occurs multiple times in the file.
              If the stranded option is not used (default), a kmer 
              and its reverse complement are considered as the same kmer.
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
      --out-sequences <OUT_SEQUENCES>
          Output file containing the filtered original sequences (eg. reads).
          It will be automatically in fasta or fastq format depending on the input file.
          If not provided, only the in_kmers with their count is output [default: ]
      --out-filelist <OUT_FILELIST>
          Output txt file containing in each line a path to a fasta or fastq [.gz] file 
          that will contain the related output file from the input files list  [default: ]
      --out-kmers <OUT_KMERS>
          If provided, output a text file containing the kmers that occur in the reads 
          with their 
           * number of occurrences 
              or 
           * their occurrence positions if the --output_kmer_positions option is used
              Note: if `--in_filelist` is used the output counted kmers are 
              those occurring the last input file of that list [default: ]
      --counted-kmer-threshold <COUNTED_KMER_THRESHOLD>
          If out_kmers is provided, output only reference kmers whose number of occurrences 
          is at least equal to this value.
          If out_kmers is not provided, this option is ignored [default: 0]
      --output-kmer-positions
          If out_kmers is provided, either only count their number of occurrences (default)
          or output their occurrence positions (read_id, position, strand)
      --output-mapping-positions
          If provided, output matching positions on sequences in the
          out_sequence file(s) 
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
  -t, --threads <THREADS>
          Number of threads
             Note: if not provided, the number of threads is set to the number of logical cores [default: 0]
  -h, --help
          Print help
  -V, --version
          Print version
```

###  5.2. <a name='Examples'></a>Examples
####  5.2.1. <a name='Basic'></a>Basic 
```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta  --out-kmers counted_kmers.txt
```

The `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `compacted_kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the estimated ratio of shared kmers and number of shared kmers.

If the `--out-kmers` option is used, the file `counted_kmers.txt` contains for each kmer in `compacted_kmers.fasta` the number of times it was found in `filtered_reads.fasta` (displays only kmers whose counts are higher than 0).

####  5.2.2. <a name='Usingfilters'></a>Using filters
```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta  --out-kmers counted_kmers.txt --min-threshold 50 --max-threshold 70
```

In this case only sequeces from `reads.fasta` that have more than 50% and at most 70% of their kmers in `compacted_kmers.fasta` are output.

####  5.2.3. <a name='Specifyingstrands'></a>Specifying strands

```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta --stranded
```
In this case, the kmers found in `compacted_kmers.fasta` are indexed in their original orientation, and kmers extracted from `reads.fasta` are queried in their original orientation. 

Note that without the `--stranded` option, all kmers (indexed and queried) are considered in their canonical form.


One may be interested in finding kmers from the reverse complement of the queried sequences. In this case we add the `--query-reverse` option together with the `--stranded` option:
```bash
back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta --stranded
```

####  5.2.4. <a name='Readingsequencesfromstandardinput:'></a>Reading sequences from standard input: 

```bash
cat reads.fasta | back_to_sequences --in-kmers compacted_kmers.fasta --out-sequences filtered_reads.fasta 
```
Do not provide the `--in-sequences` if your input data are read from stdin.

####  5.2.5. <a name='Usingseveralinputreadsets.'></a>Using several input read sets. 
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

#### 5.2.6 Output matching kmers
##### 5.2.6.1 Output the list of matching kmers with their number of occurrences
`back_to_sequences` enables to output for each kmers in `in-kmers` set, its number of occurrences in the queried sequences. 

```bash
back_to_sequences --in-sequences sequence.fa --in-kmers kmer.fa --out-sequences /dev/null  --out-kmers out_kmers.txt
```
In this case the `out_kmers.txt` file contains, for each kmer from `kmer.fa` its number of occurrences in the `sequence.fa` file (canonical or not, depending on the usage of  the `--stranded` option). 



##### 5.2.6.2 Output the list of matching kmers with their position in sequences
`back_to_sequences` enables to output for each kmers in `in-kmers` set, its positions in the queried sequences. 

```bash
back_to_sequences --in-sequences sequence.fa --in-kmers kmer.fa --out-sequences /dev/null  --out-kmers out_kmers.txt --output-kmer-positions
```
In this case the `out_kmers.txt` file contains, for each kmer from `kmer.fa` its occurrences in the `sequence.fa` file. An occurrence is given by a triplet `(sequence_id, position, strand)`.  
- `sequence_id`: id (starting from 0) of the sequence from `sequence.fa` where the kmer occurs.
- `position`: position (starting from 0) where the kmer occurs on the sequence
- `strand`: false, except when querying the reverse complement of the sequences `--query-reverse` and using the `stranded` option.

##### Outputs for each queried sequences its location and strand of shared kmers
`back_to_sequences` enables to output for each queried sequences, the location and strand of its kmers shared with the `in-kmers` set.
```bash
 back_to_sequences --in-sequences sequence.fa --in-kmers kmer.fa --out-sequences out_sequences.fa
 ```
In this case the `out_sequences.fa` contains for each queried sequence its usual header (original header number and ratio of shared kmers with the `in-kmers` set) and additionaly, it shows the location (0-based) of shared kmers. For each location (including 0), the strand is indicated by nothing or a `-` character if the `--stranded` option is given. 





##  6. <a name='Generaterandomdatafortesting'></a>Generate random data for testing
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
* [X] Add a threshold on the number of shared kmers
* [X] Parallelize the read extraction step
* [ ] Thinks about a way to adapt this to protein sequences
* [X] Add an option to set the size of the bloom filter used by kmindex
* [ ] Provide a way to index and query more than one set $K$ of kmers
* [X] Output the strand of matched kmers
