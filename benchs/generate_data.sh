echo "generate data"
python3 ../scripts/generate_random_fasta.py 1 100000000 100000000 ref_seq.fasta
# generate some kmers from the ref_seq
python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 50 --max_size 50 --num 50000 --output compacted_kmers.fasta

# generate 10000 100000 1000000 10000000 reads.
for nb_reads in 10000 100000 1000000 10000000 100000000 200000000
do
    echo "generate $nb_reads reads"
	python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 100 --max_size 100 --num ${nb_reads} --output reads_${nb_reads}.fasta 
done
