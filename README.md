# Back to sequences

![tests](https://github.com/pierrepeterlongo/back_to_sequences/workflows/tests/badge.svg)
![lints](https://github.com/pierrepeterlongo/back_to_sequences/workflows/lints/badge.svg)
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07066/status.svg)](https://doi.org/10.21105/joss.07066)

<img src="k2s.jpg" alt="Old library - IA generated" width="150" height="150">


Given a set $K$ of kmers (fasta / fastq [.gz] format) and a set of sequences  (fasta / fastq [.gz] format), this tool will extract the sequences containing some of those kmers.

A minimal ($m$) and a maximal ($M$) thresholds are proposed. A sequence whose percentage of kmers shared with $K$ are in $]m, M]$ is output with its original header + the number of shared kmers + the ratio of shared kmers:
```
>original_header 20 6.13
TGGATAAAAAGGCTGACGAAAGGTCTAGCTAAAATTGTCAGGTGCTCTCAGATAAAGCAGTAAGCGAGTTGGTGTTCGCTGAGCGTCGACTAGGCAACGTTAAAGCTATTTTAGGC...
```
In this case 20 kmers are shared with the indexed kmers. This represents 6.13% of the kmers in the sequence.

## Install
Please see https://b2s-doc.readthedocs.io/en/latest/usage.html#installation

## Simplest usage
```bash
back_to_sequences --in-kmers kmers.fasta --in-sequences reads.fasta --out-sequences filtered_reads.fasta  --out-kmers counted_kmers.txt
```
The `filtered_reads.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers from `kmers.fasta`. The headers of each read is the same as in `reads.fasta`, plus the estimated ratio of shared kmers and number of shared kmers.

As the --out-kmers option is used, the file `counted_kmers.txt` contains for each kmer in `kmers.fasta` the number of times it was found in `filtered_reads.fasta`.

## Result example

Example results obtained on 
* the GenOuest platform on a node with 32 threads Xeon 2.2 GHz, denoted by "genouest" in the table below.
* a MacBook, Apple M2 pro, 16 GB RAM, with 10 threads, denoted by "mac" in the table below.
* AMD Ryzen 7 4.2 GHz 5800X 64 GB RAM,  with 16 threads, denoted by "AMD" in the table below.

Indexed: one million kmers eacho of length 31.
We queried: from 10,000 reads to 200 million reads each of length 100.

| Number of reads | Time genouest | Time mac | Time AMD | max RAM |
|:---------------:|:-------------:|:--------:|:--------:|:-------:|
| 10,000          | 0.7s          | 0.54s    | 0.4s     | 0.13 GB |
| 100,000         | 0.8s          | 0.8s     | 1.2s     | 0.13 GB |
| 1,000,000       | 2.0s          | 3.5s     | 7.1s     | 0.13 GB |
| 10,000,000      | 7.1s          | 11s      | 16s      | 0.13 GB |
| 100,000,000     | 47s           | 58s      | 48s      | 0.13 GB |
| 200,000,000     | 1m32s         | 1m52s    | 1m44     | 0.13 GB |

See [this page](https://b2s-doc.readthedocs.io/en/latest/results.html) for details

## Basical usages and parameters
Please reafer the specific documentation for
* [basical usages](https://b2s-doc.readthedocs.io/en/latest/use%20cases.html)
* [a complete description of parameters](https://b2s-doc.readthedocs.io/en/latest/usage.html#back-to-sequences-parameters)


## Contributions
Please check out [How to contribute](CONTRIBUTING.md)

## Citations
Baire et al., (2024). Back to sequences: Find the origin of k-mers. Journal of Open Source Software, 9(101), 7066, https://doi.org/10.21105/joss.07066

bibtex:
```bib
@article{Baire2024,
  author = {Anthony Baire and Pierre Marijon and Francesco Andreace and Pierre Peterlongo},
  title = {Back to sequences: Find the origin of k-mers}, journal = {Journal of Open Source Software},
  doi = {10.21105/joss.07066},
  url = {https://doi.org/10.21105/joss.07066},
  year = {2024},
  publisher = {The Open Journal},
  volume = {9},
  number = {101},
  pages = {7066}
}
```
## Documentation
**Full documentation** is available at [https://b2s-doc.readthedocs.io/en/latest/](https://b2s-doc.readthedocs.io/en/latest/)
