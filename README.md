# Kmer2sequences

## Description
Given a set of kmers (fasta format) and a set of sequences (fasta format), this tool will extract the sequences containing the kmers.

**key idea**: 
 1. kmers (even if they are few) are indexed using kmindex. 
 2. Then reads are queried against the index, again using kmindex.
 3. Finally, the reads are extracted from the fasta file, given the kmindex output (that contains only headers)
 It is possible to set up a threshold on the ration of shared kmers between the read and the query.

## Install:
```bash
git clone https://github.com/pierrepeterlongo/kmer2sequences.git
cargo install --path . --locked
```

## Dependencies
* [kmtricks](https://github.com/tlemane/kmtricks). Note that kmtrics is provided by the conda version of kmindex.
* [kmindex](https://github.com/tlemane/kmindex)

For compiling with mac, cf the note at the end of the file.

## complete example: 
1. generate random reads and extract some of their kmers: 
```bash
# Generate 100000 reads of average length 500 and minimum length 100
python scripts/generate_random_fasta.py 100000 500 100 reads.fasta

# Extract 100 random kmers of length 31 from the reads
python3 scripts/extract_random_kmers_from_a_fasta_file.py --canonical reads.fasta 31 100 kmers.fasta

# Create the file of file, used by kmindex, containing the kmers. D is simply a prefix. 
echo D:kmers.fasta > fof.txt
```

2. index the kmers: 
```bash
back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ./bin/kmindex
```
3. search the kmers in the reads: 
```bash
back_to_sequences get_headers --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_headers headers --kmindex_path ./bin/kmindex
```

4. back to the read sequences
```bash
back_to_sequences to_reads --in_tsv_dir headers --in_fasta reads.fasta --out_fasta out.fasta --threshold 0.0
```

That's all, the `out.fasta` file contains the original sequences (here reads) from `reads.fasta` that contain at least one of the kmers in `kmers.fasta`.
The headers of each read is the same as in `reads.fasta`, plus the ratio of shared kmers.

### kmindex and kmtricks compilation note for mac os users.
kmindex install is a bit complex (sept 2023)
* kmtricks: does not compile on mac since 1.2.1. Solution: 
	* Installation via conda.
	* In case of conda issues between noarch and arm64. Solution: 
```bash
 conda config --env --set subdir osx-64
 conda create -p env_kmtricks_1.2.1 kmtricks=1.2.1 
 conda activate /Users/ppeterlo/workspace/kmer2sequences/env_kmtricks_1.2.1
```
* kmindex: commit 1b018539a4a1730a51840ed5a9330c023baf3814
	* For compiling with a mac: change lines ` asm volatile("pause");` by `asm volatile("isb" : : : "memory");` in `lib/include/kmindex/spinlock.hpp`. 
	* Comment `if (!(kmv >= min_kmv_required))` line 219 of `app/kmindex/build.cpp`
(sorry for the trouble)

