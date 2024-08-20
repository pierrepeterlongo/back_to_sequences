---
title: 'Back to sequences: find the origin of $k$-mers'
tags:
  - rust
  - kmer
  - indexing
  - genomic sequencing data
authors:
  - name: Anthony Baire
    affiliation: 1
  - name: Pierre Marijon 
    orcid: 0000-0002-6694-6873
    affiliation: 2
  - name: Francesco Andreace
    orcid: 0009-0008-0566-200X
    affiliation: 3
  - name: Pierre Peterlongo 
    orcid: 0000-0003-0776-6407
    affiliation: 1
affiliations:
 - name: Univ. Rennes, Inria, CNRS, IRISA - UMR 6074, Rennes, F-35000 France
   index: 1
 - name: Laboratoire de Biologie Médicale Multisites SeqOIA, Paris, France
   index: 2
 - name: Department of Computational Biology, Institut Pasteur, Université Paris Cité, Paris, F-75015, France
   index: 3
date: 10 June 2024
bibliography: paper.bib 


---

# Abstract
A vast majority of bioinformatics tools dedicated to the treatment of
  raw sequencing data heavily use the concept of $k$-mers, which are
  words of length $k$. This enables us to reduce the redundancy of data
  (and thus the memory pressure), to discard sequencing errors, and to
  dispose of objects of fixed size that can be easily manipulated and
  compared to each other. A drawback is that the link between each
  $k$-mer and the original set of sequences to which it belongs is lost.
  Given the volume of data considered in this context, finding back this
  association is costly. In this work, we present "`back_to_sequences`",
  a simple tool designed to index a set of $k$-mers of interest and to
  stream a set of sequences, extracting those containing at least one of
  the indexed $k$-mer. In addition, the occurrences of $k$-mers in the
  sequences can be provided. Our results show that `back_to_sequences`
  streams $\approx200$ short read per millisecond, allowing to search
  $k$-mers in hundreds of millions of reads in a matter of a few
  minutes.

# Statement of Need

In the 2010s, following the emergence of next-generation sequencing
technology, read assembly strategies based on the
overlap-layout-consensus paradigm (OLC) were unable to scale to tens of
millions of reads or more, prompting the usage of the *de Bruijn* graph
(dBG) data structure [@flicek2009sense; @schatz2010assembly]. The
success of dBG was due to the fact that the main difficulties associated
with the nature of the sequencing data (read redundancy, nonuniform
coverage and nonuniform overlap between reads, sequencing errors,
unknown sequencing strand) were complex to handle with OLC while being
easy to handle or simply solved with the dBG
approach [@li2012comparison].

Recall that in the dBG assembly approach, 1. all $k$-mers (words of
length $k$) from a set of reads are counted; 2. those with an abundance
lower than a threshold are considered to contain sequencing errors
and are discarded; 3. the remaining $k$-mers are organized in a dBG; 4.
the paths of the dBG form the basis of the assembly, later improved
thanks to scaffolding tools [@huson2002greedy] such as the tools provided,
for instance, by the Spades assembler [@bankevich2012spades].

The usefulness of $k$-mers did not end with their use in dBGs. A large
and redundant set of sequences, such as a sequencing read set, can be
summarized by its set of $k$-mers. Among multiple fundamental tasks, it
has been the basis for metagenome comparisons [@benoit2016multiple], for
taxonomy characterization [@wood2019improved], for indexing
purposes [@cracco2023extremely; @lemane2022kmtricks], for
genotyping [@grytten2022kage], for species
identification [@sarmashghi2019skmer], for transcript expression
estimation [@zhang2014rna], or for variant
discovery [@uricaru2015reference] to cite only a few.

One of the keys to the success of the use of $k$-mers is its low
resource needs. Whatever the sequencing coverage, once
filtered, the number of distinct $k$-mers is at most equal to the
original genome size. This offers a minimal impact on random access memory
(RAM) and/or disk needs. However, this comes at the cost of losing the link
between each $k$-mer and the sequence(s) from which it originates. Storing
these links explicitly would reintroduce the problem associated with the
abundance of original reads, as the link between each original read and
each of its $k$-mers would have to be stored. For instance, considering
$k$-mers from a sequencing experiment of a human genome ($\approx 3$
billion nucleotides) with a coverage of 50x (each $k$-mer occurs on
average in 50 distinct reads) would require more than 2TB of space
considering 64 bits for storing each link and 64 bits for storing the
associated read identifier. This is not acceptable.

There are many situations in which finding the origin of $k$-mers from a
set in a set of reads is informative. For instance, this is the case in
studies in which the output is composed of $k$-mers associated to
biological knowledge such as biological
variants [@uricaru2015reference], or $k$-mers specific to a phenotypic
trait [@lemane2022k] to cite a few. This approach can also be used for
quality control [@plaza2015quality] or contamination
removal [@gonzalez2023akmerbroom] for instance.

Finding back the link between a $k$-mer and each original read in which
it occurs can be performed by indexing the reads [@marchet2020resource]
which hardly scales hundreds of millions reads. One may also apply
*grep-like* evolved pattern matching approaches such as `pt` [@pt].
However, even if they were highly optimized in recent decades, these
approaches cannot efficiently detect thousands of $k$-mers in millions
of reads. The approaches that use $k$-mers for genotyping such as
`kage` [@grytten2022kage] may find the number of occurrences of $k$-mers
but not extract the sequences from which they originate.

In this context, we propose `back_to_sequences`, a tool specifically
dedicated to extracting from $\mathcal{S}$, a set of sequences (eg.
reads), those that contain some of the $k$-mers from a set $\mathcal{K}$
given as input. The occurrences of $k$-mers in each sequence of the
queried set $\mathcal{S}$ can also be output.

# Possible Alternatives

To the best of our knowledge, there exists no tool specifically
dedicated to this task, while indexing the set of $k$-mers$\mathcal{K}$.

Genotypers as `kage` [@grytten2022kage], using
`kmer mapper` [@kmermapper], provide a way to count the number of
occurrences of each $k$-mer from a set of reference $k$-mers in a read
file. However, they do not offer the feature to extract reads that
contain any reference $k$-mer.

Finding one unique $k$-mer of interest in a set of sequences can be done
using the classical `grep` or more recent pattern-matching tools such as
"*The Platinum Searcher*" [@pt] or "*The Silver Searcher*" [@ag].

As for testing, on the MacBook, Apple M2 pro, we queried one $k$-mer in
the $\mathcal{S}_{100M}$ dataset (see
Section [2.1](#ssec:random_banch){reference-type="ref"
reference="ssec:random_banch"}) using `grep`, *The Platinum Searcher*,
and *The Silver Searcher*.

-   `grep` required 44 seconds. Thus, by simple extrapolation, searching
    for one million $k$-mers on a single computer would require
    approximately 500 days, to be compared to 5 to 6 minutes using
    `back_to_sequences`.

-   `pt` (*The Platinum Searcher*) required 15 seconds, which can be
    extrapolated to approximately 175 days if searching for one million
    $k$-mers.

-   `ag` (*The Silver Searcher*) did not finish after 400 seconds.

Summing up, these alternative tools are not meant for querying numerous
patterns at the same time and do not scale to large problem instances.

Note also that these alternative tools are not specialized for genomic
data in which one is interested in searching for a $k$-mer and
potentially its reverse complement. Finally, these tools do not easily
provide the number of occurrences or occurrence positions of each of the
searched patterns when there are many.

# Conclusion

We believe that `back_to_sequences` is a generic and handy tool that
will be beneficial for building pipelines that require manipulating
$k$-mers and finding back the sequences from which they originate and/or
counting their number of occurrences in a set of genomic sequences. We
also believe that `back_to_sequences` will have other straightforward
applications, such as quality control, contamination removal, or
genotyping known pieces of sequences in raw sequencing datasets, all of
which being possible in real time throughout the sequencing process.

# Acknowledgements

We acknowledge the GenOuest core facility (<https://www.genouest.org>)
for providing the computing infrastructure. The work was funded by ANR
SeqDigger (ANR-19-CE45-0008), and by the Inria Challenge "OmicFinder"
(<https://project.inria.fr/omicfinder/>).

[^1]: <http://www.ebi.ac.uk/ena/data/view/ERS488262>

# References
