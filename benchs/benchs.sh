export TIME="\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem"

# Run the benchs
# Each time check the time to query them
echo
echo "Running basic benchs"
for nb_reads in 10000 100000 1000000 10000000 100000000 200000000
do
    echo
    echo "nb_reads: $nb_reads"
    time back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads_${nb_reads}.fasta  --out-sequences filtered_reads_${nb_reads}.fasta  -k 31 --out-kmers counted_kmers_${nb_reads}.txt --counted-kmer-threshold 1 > /dev/null
    rm -f filtered_reads_${nb_reads}.fasta counted_kmers_${nb_reads}.txt
done

# Run the bench using the --output-kmer-positions option
echo
echo "Running benchs with --output-kmer-positions"
for nb_reads in 10000 100000 1000000 10000000 100000000 200000000
do
    echo 
    echo "nb_reads: $nb_reads"
    time back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads_${nb_reads}.fasta  --out-sequences filtered_reads_${nb_reads}.fasta  -k 31 --out-kmers counted_kmers_${nb_reads}.txt --output-kmer-positions --counted-kmer-threshold 1 > /dev/null
    rm -f filtered_reads_${nb_reads}.fasta counted_kmers_${nb_reads}.txt
done