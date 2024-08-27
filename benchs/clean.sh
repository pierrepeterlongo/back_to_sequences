#!/bin/bash

# Clean all generated data by generate_data.sh

# List of files and directories to be removed
generated_files=(
    "ref_seq.fasta"
    "compacted_kmers.fasta"
    "reads_[0-9]*.fasta"
)

# Loop through the list and remove each item
for item in "${generated_files[@]}"; do
    cmd="rm -f $item"
    $cmd && echo "Removed $item" || echo "Failed to remove $item"
done

echo "Cleanup complete."