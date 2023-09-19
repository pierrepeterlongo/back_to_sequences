# indexes 100000 kmers of length 31

# Gerenate random reads
python3 ../scripts/generate_random_fasta.py 100000 500 100 tmp.fasta
# Extract 100000 kmers of length 31 from the reads
python3 ../scripts/extract_random_kmers_from_a_fasta_file.py --canonical tmp.fasta 31 100000 kmers.fasta

echo D:kmers.fasta > fof.txt

echo "Indexing kmers"
time back_to_sequences index_kmers --in_kmers fof.txt --out_index indexed_kmers -k 31 --kmindex_path ../bin/kmindex


rm -rf tmp.fasta

# Run the benchs
# generate 1 10 100 1000 10000 100000 1000000 10000000 reads.
# Each time check the time to query them

for nb_reads in 1 10 100 1000 10000 100000 1000000 10000000 100000000
do
    # echo "Generating $nb_reads reads"
    python3 ../scripts/generate_random_fasta.py $nb_reads 500 100 reads_${nb_reads}.fasta
    
    echo "Querying $nb_reads reads"
    time back_to_sequences query_sequences --in_sequences reads_${nb_reads}.fasta --in_kmer_index indexed_kmers --out_fasta filtered_readsreads_${nb_reads}.fasta --kmindex_path ../bin/kmindex  
done

# clean all:
rm -f filtered_readsreads_1*.fasta
rm -f reads_1*.fasta
rm -f kmers.fasta fof.txt
rm -rf indexed_kmers
rm -rf local_kmer_index
rm -rf tmp_headers 