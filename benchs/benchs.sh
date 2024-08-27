export TIME="\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem"

## Run the benchs once data are generated

# Usage: bash benchs.sh [max_nb_reads]
# There if zero or one argument: 
# - If zero, it tests the full data set
# - If one, it tests the data set up to the specified number of reads
# check the number of arguments
if [ "$#" -gt 1 ]; then
    echo "Usage: bash benchs.sh [max_nb_reads]"
    exit 1
fi

# Get the potential max number of reads
if [ "$#" -eq 1 ]; then
    max_nb_reads=$1
else
    max_nb_reads=200000000
fi

# Run the benchs
# Each time check the time to query them
echo
echo "Running basic benchs"
for nb_reads in 10000 100000 1000000 10000000 100000000 200000000
do
    if [ $nb_reads -gt $max_nb_reads ]; then
        break
    fi
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
    if [ $nb_reads -gt $max_nb_reads ]; then
        break
    fi
    echo 
    echo "nb_reads: $nb_reads"
    time back_to_sequences --in-kmers compacted_kmers.fasta --in-sequences reads_${nb_reads}.fasta  --out-sequences filtered_reads_${nb_reads}.fasta  -k 31 --out-kmers counted_kmers_${nb_reads}.txt --output-kmer-positions --counted-kmer-threshold 1 > /dev/null
    rm -f filtered_reads_${nb_reads}.fasta counted_kmers_${nb_reads}.txt
done