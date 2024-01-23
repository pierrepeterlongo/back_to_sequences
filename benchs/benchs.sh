export TIME="\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem"

# Run the benchs
# Each time check the time to query them
for nb_reads in 10000 100000 1000000 10000000 100000000 200000000
do
    echo "nb_reads: $nb_reads"
    time back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads_${nb_reads}.fasta  --out-sequences filtered_reads_${nb_reads}.fasta  -k 31 --out-kmers counted_kmers_${nb_reads}.txt
    
    time back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads_${nb_reads}.fasta  -k 31 --out-kmers counted_kmers_only_${nb_reads}.txt

    # check that counted_kmers_${nb_reads}.txt counted_kmers_only_${nb_reads}.txt are the same
    diff counted_kmers_${nb_reads}.txt counted_kmers_only_${nb_reads}.txt > /dev/null
    if [ $? -eq 0 ]
    then
        echo "OK"
    else
        echo "KO"
    fi



    
done
