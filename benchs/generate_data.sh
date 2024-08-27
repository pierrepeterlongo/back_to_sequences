## Generate data for the benchs
# This script generates the data used in the benchs.sh script.
# It generates a reference sequence, a set of kmers from the reference sequence and some set of reads from the reference sequence.
# If the user specifies a max number of reads, it will generate read files up to the specified limit of number of reads

# Usage: bash generate_data.sh [max_nb_reads]
# There if zero or one argument: 
# - If zero, it generates the full data set
# - If one, it generates the data set up to the specified number of reads
# check the number of arguments
if [ "$#" -gt 1 ]; then
    echo "Usage: bash generate_data.sh [max_nb_reads]"
    exit 1
fi

# Get the potential max number of reads
if [ "$#" -eq 1 ]; then
    max_nb_reads=$1
else
    max_nb_reads=200000000
fi


max_nb_reads=$1



echo "generate reference sequence"
python3 ../scripts/generate_random_fasta.py 1 100000000 100000000 ref_seq.fasta

echo "generate kmers"
python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 50 --max_size 50 --num 50000 --output compacted_kmers.fasta

# generate 10000 100000 1000000 10000000 100000000 200000000 reads.
for nb_reads in 10000 100000 1000000 10000000 100000000 200000000
do
    if [ $nb_reads -gt $max_nb_reads ]; then
        break
    fi
    echo "generate $nb_reads reads"
	python3 ../scripts/extract_random_sequences.py --input ref_seq.fasta --min_size 100 --max_size 100 --num ${nb_reads} --output reads_${nb_reads}.fasta 
done
