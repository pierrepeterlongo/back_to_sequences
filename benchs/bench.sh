# indexes 100000 kmers of length 31
export TIME="\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem"
echo "generate data"
# python3 ../scripts/generate_random_fasta.py 1 100000000 100000000 ref_seq.fasta
# generate some kmers from the ref_seq
# python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 50 --max_size 50 --num 50000 --output compacted_kmers.fasta
echo ref_set:compacted_kmers.fasta > fof.txt

echo "Indexing kmers"
time ./disk_mem_count.sh back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ../bin/kmindex


# Run the benchs
# generate 1 10 100 1000 10000 100000 1000000 10000000 reads.
# Each time check the time to query them

for nb_reads in 1 10 100 1000 10000 100000 1000000 10000000 100000000
do
    echo "generate $nb_reads reads"
    python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 100 --max_size 500 --num ${nb_reads} --output reads_${nb_reads}.fasta 

    echo "Querying $nb_reads reads"
    time ./disk_mem_count.sh back_to_sequences query_sequences --in_sequences reads_${nb_reads}.fasta --in_kmer_index indexed_kmers --out_fasta filtered_reads_${nb_reads}.fasta -t 32 --kmindex_path ../bin/kmindex


    echo "Validation and kmer counting"
    time ./disk_mem_count.sh back_to_sequences exact_count --in_kmers compacted_kmers.fasta --in_fasta filtered_reads_${nb_reads}.fasta  --out_fasta filtered_reads_exact_${nb_reads}.fasta  -k 31 --out_counted_kmers counted_kmers_${nb_reads}.txt
done
