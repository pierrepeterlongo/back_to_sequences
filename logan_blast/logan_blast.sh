#!/bin/bash


ACCESSION_FILE=""
QUERY_FILE=""
DELETE=false
UNITIGS=false
LIMIT=0
KMER_SIZE=17

print_help() {
    echo "Usage: $0 --accessions <accessions.txt> --query <query_file.fa> [--delete] [--kmer-size <k>] [--limit <n>]"
    echo "Options:"
    echo "  -a, --accessions  Path to accessions.txt file. Containing one accession per line)"
    echo "  -q, --query       Path to query fasta file"
    echo "  -u, --unitigs     Consider the unitigs verison of the accessions instead of contigs"
    echo "  -k, --kmer-size K-mer size for sequence recruitment with back_to_sequences (default: ${KMER_SIZE})"
    echo "  -l, --limit     Limit number of accessions to process from accession file (default: no limit)"
    echo "  -d, --delete    Delete recruited accessions and accessions files after processing (default: keep all files)"
    echo "  -h, --help      Show this help message"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -a|--accessions)
            ACCESSION_FILE="$2"
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
        -u|--unitigs)
            UNITIGS=true
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

if [[ -z "$ACCESSION_FILE" || -z "$QUERY_FILE" ]]; then
    echo  "\033[0;31mError: --accessions and --query are required.\033[0m"
    print_help
    exit 1
fi

if [ ! -f "$QUERY_FILE" ]; then
    echo  "\033[0;31mError: Query file '$QUERY_FILE' does not exist.\033[0m"
    exit 1
fi

if [ ! -f "$ACCESSION_FILE" ]; then
    echo  "\033[0;31mError: Accessions file '$ACCESSION_FILE' does not exist.\033[0m"
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

type=contig
if [ "$UNITIGS" = true ]; then
    type=unitig
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
	if [ ! -f "${accession}.${type}s.fa.zst" ]; then
		echo "\033[0;33mDownloading ${accession}.${type}s.fa.zst...\033[0m"
        if [ "$UNITIGS" = false ]; then
		    aws s3 cp s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst . --no-sign-request
            else 
		    aws s3 cp s3://logan-pub/u/${accession}/${accession}.unitigs.fa.zst . --no-sign-request
        fi

        if [ $? -ne 0 ]; then
            echo "\033[0;31mError: Failed to download ${accession}.${type}s.fa.zst from S3. Skipping this accession.\033[0m"
            if [ "$UNITIGS" = false ]; then
                echo "\033[0;31mThis usually occurs because some logan accessions do not have contigs files (they have only unitigs files).\033[0m"
            fi
            continue
        fi
	else
		echo "\033[0;33mUsing existing local version of ${accession}.${type}s.fa.zst...\033[0m"
	fi
	echo "\033[0;33mRecruiting sequences from ${accession}.${type}s.fa.zst with a match with ${QUERY_FILE}...\033[0m"
	echo "\033[0;32mback_to_sequences --kmer-size ${KMER_SIZE} --in-kmers ${QUERY_FILE} --in-sequences ${accession}.${type}s.fa.zst --out-sequences T.fa\033[0m"
	back_to_sequences --kmer-size ${KMER_SIZE} --in-kmers ${QUERY_FILE} --in-sequences  ${accession}.${type}s.fa.zst --out-sequences ${accession}.recruited_${type}s.fa > /dev/null 2>&1

    ## check if any sequences were recruited
    if [ ! -s "${accession}.recruited_${type}s.fa" ]; then
        echo "\033[0;33m\tNo sequences were recruited from ${accession}.${type}s.fa.zst. Skipping BLAST step.\033[0m"
        if [ "$DELETE" = true ]; then
            echo "\033[0;33mDeleting ${accession}.recruited_${type}s.fa and ${accession}.${type}s.fa.zst...\033[0m"
            rm -f ${accession}.recruited_${type}s.fa
            rm -f ${accession}.${type}s.fa.zst
        fi
        continue
    fi

	echo "\033[0;33mAligning recruited sequences from ${accession}.${type}s.fa.zst with ${QUERY_FILE}...\033[0m"
	echo "\033[0;32mpython local_blast.py --target ${QUERY_FILE} --query ${accession}.recruited_${type}s.fa \033[0m"
	python local_blast.py --target ${QUERY_FILE} --query ${accession}.recruited_${type}s.fa
    if [ "$DELETE" = true ]; then
        echo "\033[0;33mDeleting ${accession}.recruited_${type}s.fa and ${accession}.${type}s.fa.zst...\033[0m"
        rm -f ${accession}.recruited_${type}s.fa
        rm -f ${accession}.${type}s.fa.zst
    fi
done < ${ACCESSION_FILE}

echo

    echo  "\n\033[1;34m================\033[0m"
    echo  "\033[1;36m>>> All done <<<\033[0m"
    echo  "\033[1;34m================\033[0m\n"
# if --delete was not used, show user how to remove all intermediate files
if [ "$DELETE" = false ]; then
    echo "\033[0;33mYou  did not use --delete option. So you can manually remove all intermediate files (recruited ${type}s and ${type}s files) by running:\033[0m"
    echo "\033[0;32mrm -f *.recruited_${type}s.fa *.${type}s.fa.zst\033[0m"
fi

# QUERY_FILE base name without extension
# Remove all extensions from QUERY_FILE to get the base name
QUERY_BASENAME=$(basename "$QUERY_FILE")
QUERY_BASENAME="${QUERY_BASENAME%%.*}"
echo
echo "\033[0;33mResults can be found in directories <accessions>.recruited_${type}s_vs_${QUERY_BASENAME} for each accession id in the ${ACCESSION_FILE}.\033[0m"