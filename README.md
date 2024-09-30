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

## Basical usages and parameters
Please reafer the specific documentation for
* [basical usages](https://b2s-doc.readthedocs.io/en/latest/use%20cases.html)
* [a complete description of parameters](https://b2s-doc.readthedocs.io/en/latest/usage.html#back-to-sequences-parameters)


## Contributions
Please check out [How to contribute](CONTRIBUTING.md)


## Documentation
**Full documentation** is available at [https://b2s-doc.readthedocs.io/en/latest/](https://b2s-doc.readthedocs.io/en/latest/)
