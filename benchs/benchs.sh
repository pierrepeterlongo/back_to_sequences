# indexes 100000 kmers of length 31
export TIME="\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem"
echo "generate data"
# python3 ../scripts/generate_random_fasta.py 1 100000000 100000000 ref_seq.fasta
# generate some kmers from the ref_seq
# python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 50 --max_size 50 --num 50000 --output compacted_kmers.fasta
echo ref_set:compacted_kmers.fasta > fof.txt

# Run the benchs
# generate 10000 100000 1000000 10000000 reads.
# Each time check the time to query them

for nb_reads in 10000 100000 1000000 10000000 100000000
do
    echo "generate $nb_reads reads"
#    python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 100 --max_size 500 --num ${nb_reads} --output reads_${nb_reads}.fasta 

    echo "Validation and kmer counting"
    time ./disk_mem_count.sh back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads_${nb_reads}.fasta  --out-sequences filtered_reads_${nb_reads}.fasta  -k 31 --out-kmers counted_kmers_${nb_reads}.txt
done
