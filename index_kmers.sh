#!/bin/bash

# function usage ()
usage()
{
  echo "Usage: "
  echo "  ./${0} [-i str] [-o str] [-k int] [-t int] [-h]"
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


# clean potential previous runs
rm -rf kmindex_dir ${output_index_file_path}

# get the path of the current directory:
current_dir_path=`pwd`

# run kmindex on the input file:
cmd="${current_dir_path}/bin/kmindex build -f ${input_file_path} \
        --run-dir kmindex_dir --index ${output_index_file_path} \
        --register-as kmers -k ${kmer_size} -t ${nb_threads} \
        --bloom-size 30000000 \
        --hard-min 1 "
echo "Running: ${cmd}"
${cmd}
# get the return code:
return_code=$?
# check the return code:
if [ ${return_code} -ne 0 ]; then
  echo "Error: kmindex returned with error code ${return_code}"
  exit 1
else
  echo "Indexed kmers are in file ${output_index_file_path}"
  exit 1
fi


