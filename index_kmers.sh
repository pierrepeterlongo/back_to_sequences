#!/bin/bash

# function usage ()
usage()
{
  echo "Usage: "
  echo "  ./index_kmers.sh [-i str] [-o str] [-k int] [-t int] [-h]"
  echo "Options: "
  echo "  -i <path> -> input file path"
  echo "  -o <path> -> output index file path, default: indexed_kmers"
  echo "  -k <int>  -> kmer size, default: 31"
  echo "  -t <int>  -> number of threads, default: 8"
  echo "  -h        -> show help"
  exit 1
}

# set the options: 
# -i <path> -> input file path
# -o <path> -> output index file path
# -k <int>  -> kmer size
# -t <int>  -> number of threads
# -h        -> show help

# set the default values:
input_file_path=""
output_index_file_path="indexed_kmers"
kmer_size=31
nb_threads=8

# parse the options:
while getopts "i:o:k:t:h" option; do
  case "$option" in
    i)
      input_file_path=${OPTARG}
      ;;
    o)
      output_index_file_path=${OPTARG}
      ;;
    k)
      kmer_size=${OPTARG}
      ;;
    t)
      nb_threads=${OPTARG}
      ;;
    h)
      usage
      exit 1
      ;;
    *)
      usage
      exit 1
      ;;
  esac
done

# check that input_file_path is not empty:
if [ -z "${input_file_path}" ]; then
  echo "Error: input file path is empty"
  usage
  exit 1
fi

# run kmindex on the input file:
./bin/kmindex -i ${input_file_path} -o ${output_index_file_path} -k ${kmer_size} -t ${nb_threads}