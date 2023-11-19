back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads.fasta > /dev/null
back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_0_1.fasta --min-threshold 0 --max-threshold 1 > /dev/null
back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_0_10.fasta --min-threshold 0 --max-threshold 10 > /dev/null
back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_0_100.fasta --min-threshold 5 --max-threshold 100 > /dev/null
back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_stranded.fasta --stranded > /dev/null
back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_stranded_qreverse.fasta --stranded --query-reverse > /dev/null
back_to_sequences --in-sequences sequence.fa --in-kmers lower_kmer.fa --out-sequences filtered_reads_lower.fasta > /dev/null
back_to_sequences --in-sequences sequence.fa --in-kmers lower_kmer_N.fa --out-sequences filtered_reads_lowerN.fasta > /dev/null

# Original generation of sequences
# back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_ref.fasta 
# back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_ref_0_1.fasta --min-threshold 0 --max-threshold 1
# back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_ref_0_10.fasta --min-threshold 0 --max-threshold 10
# back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_ref_0_100.fasta --min-threshold 5 --max-threshold 100
# back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_ref_stranded.fasta --stranded
# back_to_sequences --in-sequences sequence.fa --in-kmers one_kmer.fa --out-sequences filtered_reads_ref_stranded_qreverse.fasta --stranded --query-reverse
# back_to_sequences --in-sequences sequence.fa --in-kmers lower_kmer.fa --out-sequences filtered_reads_ref_lower.fasta 
# back_to_sequences --in-sequences sequence.fa --in-kmers lower_kmer_N.fa --out-sequences filtered_reads_ref_lowerN.fasta 

diff filtered_reads.fasta filtered_reads_ref.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads.fasta and filtered_reads_ref.fasta are different"
    exit 1
fi

diff filtered_reads_0_1.fasta filtered_reads_ref_0_1.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_0_1.fasta and filtered_reads_ref_0_1.fasta are different"
    exit 1
fi

diff filtered_reads_0_10.fasta filtered_reads_ref_0_10.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_0_10.fasta and filtered_reads_ref_0_10.fasta are different"
    exit 1
fi

diff filtered_reads_0_100.fasta filtered_reads_ref_0_100.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_0_100.fasta and filtered_reads_ref_0_100.fasta are different"
    exit 1
fi

diff filtered_reads_stranded.fasta filtered_reads_ref_stranded.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_stranded.fasta and filtered_reads_ref_stranded.fasta are different"
    exit 1
fi

diff filtered_reads_stranded_qreverse.fasta filtered_reads_ref_stranded_qreverse.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_stranded_qreverse.fasta and filtered_reads_ref_stranded_qreverse.fasta are different"
    exit 1
fi

diff filtered_reads_ref_lowerN.fasta filtered_reads_lowerN.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_ref_lowerN.fasta and filtered_reads_lowerN.fasta are different"
    exit 1
fi

diff filtered_reads_ref_lower.fasta filtered_reads_lower.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_ref_lower.fasta and filtered_reads_lower.fasta are different"
    exit 1
fi

rm -f filtered_reads.fasta filtered_reads_0_1.fasta filtered_reads_0_10.fasta filtered_reads_0_100.fasta filtered_reads_stranded.fasta filtered_reads_stranded_qreverse.fasta filtered_reads_lower.fasta filtered_reads_lowerN.fasta




cat sequence.fa | back_to_sequences --in-kmers one_kmer.fa --out-sequences filtered_reads.fasta > /dev/null
cat sequence.fa | back_to_sequences --in-kmers one_kmer.fa --out-sequences filtered_reads_0_1.fasta --min-threshold 0 --max-threshold 1 > /dev/null
cat sequence.fa | back_to_sequences --in-kmers one_kmer.fa --out-sequences filtered_reads_0_10.fasta --min-threshold 0 --max-threshold 10 > /dev/null
cat sequence.fa | back_to_sequences --in-kmers one_kmer.fa --out-sequences filtered_reads_0_100.fasta --min-threshold 5 --max-threshold 100 > /dev/null
cat sequence.fa | back_to_sequences --in-kmers one_kmer.fa --out-sequences filtered_reads_stranded.fasta --stranded > /dev/null
cat sequence.fa | back_to_sequences --in-kmers one_kmer.fa --out-sequences filtered_reads_stranded_qreverse.fasta --stranded --query-reverse > /dev/null


diff filtered_reads.fasta filtered_reads_ref.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads.fasta and filtered_reads_ref.fasta are different"
    exit 1
fi

diff filtered_reads_0_1.fasta filtered_reads_ref_0_1.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_0_1.fasta and filtered_reads_ref_0_1.fasta are different"
    exit 1
fi

diff filtered_reads_0_10.fasta filtered_reads_ref_0_10.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_0_10.fasta and filtered_reads_ref_0_10.fasta are different"
    exit 1
fi

diff filtered_reads_0_100.fasta filtered_reads_ref_0_100.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_0_100.fasta and filtered_reads_ref_0_100.fasta are different"
    exit 1
fi

diff filtered_reads_stranded.fasta filtered_reads_ref_stranded.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_stranded.fasta and filtered_reads_ref_stranded.fasta are different"
    exit 1
fi

diff filtered_reads_stranded_qreverse.fasta filtered_reads_ref_stranded_qreverse.fasta
if [ $? -ne 0 ]; then
    echo "filtered_reads_stranded_qreverse.fasta and filtered_reads_ref_stranded_qreverse.fasta are different"
    exit 1
fi

rm -f filtered_reads.fasta filtered_reads_0_1.fasta filtered_reads_0_10.fasta filtered_reads_0_100.fasta filtered_reads_stranded.fasta filtered_reads_stranded_qreverse.fasta


echo "All tests passed"