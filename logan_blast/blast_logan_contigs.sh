#!/bin/bash


CONTIGS_FILE=""
QUERY_FILE=""
DELETE=false
LIMIT=0
KMER_SIZE=17

print_help() {
    echo "Usage: $0 --contigs <contigs.txt> --query <query_file.fa> [--delete] [--kmer-size <k>] [--limit <n>]"
    echo "Options:"
    echo "  -c, --contigs   Path to contigs.txt file. Containing one accession per line)"
    echo "  -q, --query     Path to query fasta file"
    echo "  -k, --kmer-size K-mer size for sequence recruitment with back_to_sequences (default: ${KMER_SIZE})"
    echo "  -l, --limit     Limit number of accessions to process from contig file (default: no limit)"
    echo "  -d, --delete    Delete recruited contigs and contigs files after processing (default: keep all files)"
    echo "  -h, --help      Show this help message"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--contigs)
            CONTIGS_FILE="$2"
            shift 2
            ;;
        -q|--query)
            QUERY_FILE="$2"
            shift 2
            ;;
        -k|--kmer-size)
            KMER_SIZE="$2"
            shift 2
            ;;
        -l|--limit)
            LIMIT="$2"
            shift 2
            ;;
        -d|--delete)
            DELETE=true
            shift
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        *)
            echo  "\033[0;31mUnknown option: $1\033[0m"
            print_help
            exit 1
            ;;
    esac
done

if [[ -z "$CONTIGS_FILE" || -z "$QUERY_FILE" ]]; then
    echo  "\033[0;31mError: --contigs and --query are required.\033[0m"
    print_help
    exit 1
fi

if [ ! -f "$QUERY_FILE" ]; then
    echo  "\033[0;31mError: Query file '$QUERY_FILE' does not exist.\033[0m"
    exit 1
fi

if [ ! -f "$CONTIGS_FILE" ]; then
    echo  "\033[0;31mError: Contigs file '$CONTIGS_FILE' does not exist.\033[0m"
    exit 1
fi

if ! [[ "$KMER_SIZE" =~ ^[0-9]+$ ]] || [ "$KMER_SIZE" -le 0 ]; then
    echo  "\033[0;31mError: K-mer size must be a positive integer.\033[0m"
    exit 1
fi

if ! [[ "$LIMIT" =~ ^[0-9]+$ ]] || [ "$LIMIT" -lt 0 ]; then
    echo  "\033[0;31mError: Limit must be a non-negative integer.\033[0m"
    exit 1
fi

## check that back_to_sequences is installed
if ! command -v back_to_sequences &> /dev/null; then
    echo  "\033[0;31mError: back_to_sequences could not be found. Please install it first https://github.com/pierrepeterlongo/back_to_sequences.\033[0m"
    exit 1
fi

counter=0
while read accession; do
    if [ "$LIMIT" -ne 0 ] && [ "$counter" -ge "$LIMIT" ]; then
        echo  "\n\033[0;33mReached limit of $LIMIT accessions. Stopping further processing.\033[0m"
        break
    fi
    counter=$((counter + 1))
    echo  "\n\033[1;34m========================================\033[0m"
    echo  "\033[1;36m>>> Processing accession: ${accession} <<<\033[0m"
    echo  "\033[1;34m========================================\033[0m"
	if [ ! -f "${accession}.contigs.fa.zst" ]; then
		echo "\033[0;33mDownloading ${accession}.contigs.fa.zst...\033[0m"
		aws s3 cp s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst . --no-sign-request
        if [ $? -ne 0 ]; then
            echo "\033[0;31mError: Failed to download ${accession}.contigs.fa.zst from S3. Skipping this accession.\033[0m"
            echo "\033[0;31mThis usually occurs because some logan accessions do not have contigs files (they have only unitigs files).\033[0m"
            continue
        fi
	else
		echo "\033[0;33mUsing existing local version of ${accession}.contigs.fa.zst...\033[0m"
	fi
	echo "\033[0;33mRecruiting sequences from ${accession}.contigs.fa.zst with a match with ${QUERY_FILE}...\033[0m"
	echo "\033[0;32mback_to_sequences --kmer-size ${KMER_SIZE} --in-kmers ${QUERY_FILE} --in-sequences ${accession}.contigs.fa.zst --out-sequences T.fa\033[0m"
	back_to_sequences --kmer-size ${KMER_SIZE} --in-kmers ${QUERY_FILE} --in-sequences  ${accession}.contigs.fa.zst --out-sequences ${accession}.recruited_contigs.fa > /dev/null 2>&1

    ## check if any sequences were recruited
    if [ ! -s "${accession}.recruited_contigs.fa" ]; then
        echo "\033[0;33m\tNo sequences were recruited from ${accession}.contigs.fa.zst. Skipping BLAST step.\033[0m"
        if [ "$DELETE" = true ]; then
            echo "\033[0;33mDeleting ${accession}.recruited_contigs.fa and ${accession}.contigs.fa.zst...\033[0m"
            rm -f ${accession}.recruited_contigs.fa
            rm -f ${accession}.contigs.fa.zst
        fi
        continue
    fi

	echo "\033[0;33mAligning recruited sequences from ${accession}.contigs.fa.zst with ${QUERY_FILE}...\033[0m"
	echo "\033[0;32mpython local_blast.py --target ${QUERY_FILE} --query ${accession}.recruited_contigs.fa \033[0m"
	python local_blast.py --target ${QUERY_FILE} --query ${accession}.recruited_contigs.fa
    if [ "$DELETE" = true ]; then
        echo "\033[0;33mDeleting ${accession}.recruited_contigs.fa and ${accession}.contigs.fa.zst...\033[0m"
        rm -f ${accession}.recruited_contigs.fa
        rm -f ${accession}.contigs.fa.zst
    fi
done < ${CONTIGS_FILE}

echo

    echo  "\n\033[1;34m================\033[0m"
    echo  "\033[1;36m>>> All done <<<\033[0m"
    echo  "\033[1;34m================\033[0m\n"
# if --delete was not used, show user how to remove all intermediate files
if [ "$DELETE" = false ]; then
    echo "\033[0;33mYou  did not use --delete option. So you can manually remove all intermediate files (recruited contigs and contigs files) by running:\033[0m"
    echo "\033[0;32mrm -f *.recruited_contigs.fa *.contigs.fa.zst\033[0m"
fi

echo
echo "\033[0;33mResults can be found in directories ${QUERY_FILE}_vs_<accessions>.recruited_contigs for each accession id in the ${CONTIGS_FILE}.\033[0m"