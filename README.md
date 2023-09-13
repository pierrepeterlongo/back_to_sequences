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
cd scripts
python generate_random_fasta.py 100 500 100 reads.fasta
python3 extract_random_kmers_from_a_fasta_file.py --canonical reads.fasta 31 10 kmers.fasta
cd -
```
2. index the kmers: 
```bash
sh index_kmers.sh -i fof.txt -k 31
```
3. search the kmers in the reads: 
```bash
sh query_reads.sh -i indexed_kmers -q scripts/reads.fasta -o output
```

