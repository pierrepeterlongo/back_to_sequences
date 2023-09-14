# Kmer2sequences

## Install, struggling with mac os as of 13 set 2023
* kmtricks: compile plus depuis 1.2.1. 
	* Installation via conda.
	* conda: soucis noarch et arm64. Solution: 
```bash
 conda config --env --set subdir osx-64
 conda create -p env_kmtricks_1.2.1 kmtricks=1.2.1 
 conda activate /Users/ppeterlo/workspace/kmer2sequences/env_kmtricks_1.2.1
```

* kmindex: commit 1b018539a4a1730a51840ed5a9330c023baf3814
	* Compilation sous mac: changer les lignes ` asm volatile("pause");` par `asm volatile("isb" : : : "memory");` dans `lib/include/kmindex/spinlock.hpp`. 
	* Commenter le `if (!(kmv >= min_kmv_required))` ligne 219 de `app/kmindex/build.cpp`

ouf...


## complete example: 
1. generate random reads and extract some of their kmers: 
```bash
# Generate 100000 reads of average length 500 and minimum length 100
python scripts/generate_random_fasta.py 100000 500 100 reads.fasta

# Extract 10 random kmers of length 31 from the reads
python3 scripts/extract_random_kmers_from_a_fasta_file.py --canonical reads.fasta 31 100 kmers.fasta
# Create the fof file: 
echo D:kmers.fasta > fof.txt
```

2. index the kmers: 
```bash
cargo run -- index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ./bin/kmindex
```
3. search the kmers in the reads: 
```bash
cargo run -- get_headers --in_sequences reads.fasta --in_kmer_index indexed_kmers --out_headers headers --kmindex_path ./bin/kmindex
```

4. back to the read sequences
```bash
cargo run -- to_reads --in_tsv_dir headers --in_fasta reads.fasta --out_fasta out.fasta --threshold 0.0
```

