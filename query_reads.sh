#!/bin/bash

# usage function:
usage()
{
  echo "Usage: "
  echo "  ./${0} [-i str] [-o str] [-q str] [-h]"
  echo "Options: "
  echo "  -i <path> -> indexed kmers directory path"
  echo "  -q <path> -> queried sequences file path"
  echo "  -o <path> -> output file path, default: query_results"
  echo "  -h        -> show help"
  exit 1
}

# set the default options:
indexed_kmers_file_path="indexed_kmers"
queried_sequences_file_path=""
output_file_path="query_results"

# parse the options:
while getopts "i:o:q:h" option; do
  case "$option" in
    i)
      indexed_kmers_file_path=${OPTARG}
      ;;
    o)
      output_file_path=${OPTARG}
      ;;
    q)
      queried_sequences_file_path=${OPTARG}
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

# check that indexed_kmers_file_path is an non empty directory:
if [ ! -d "${indexed_kmers_file_path}" ]; then
  echo "Error: ${indexed_kmers_file_path} is not a non empty directory"
  usage
  exit 1
fi

# check that queried_sequences_file_path is not empty:
if [ -z "${queried_sequences_file_path}" ]; then
  echo "Error: queried sequences file path is empty"
  usage
  exit 1
fi

# clena previous runs: 
rm -rf ${output_file_path}

# get the path of the current directory:
current_dir_path=`pwd`

# create the kmindex command: 
kmindex_cmd="${current_dir_path}/bin/kmindex query --index ${indexed_kmers_file_path} -q ${queried_sequences_file_path} -o ${output_file_path} -n kmers -r 0.0001 --format matrix"

# show the command:
echo ${kmindex_cmd}

# run the command:
eval ${kmindex_cmd} 

# check that the command has succeeded:
if [ $? -ne 0 ]; then
  echo "Error: kmindex query command failed"
  exit 1
else
  echo "kmindex query command succeeded, results are in ${output_file_path}"
fi
