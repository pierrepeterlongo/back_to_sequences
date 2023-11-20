## Various scripts

* For generating simulated data: 
  
  * `generate_random_fasta.py`: generates num_sequences random sequence over the ACGT alphabet, fasta format. The user defines the average and minimal size of sequences
    `generate_random_fasta.py num_sequences average_size min_size output_file`
  
  * `extract_random_sequences.py`  Extracts sequences from a fasta of fastq file. 
    `extract_random_sequences.py [-h] --input INPUT_FILE --min_size MIN_SIZE --max_size MAX_SIZE --num NUM_SEQUENCES --output OUTPUT_FILE` 
  
  * `yield_reads.py` and `revcomp.py` are for internal use.

* For parsing back_to_sequences results: 
  
  * `filter_reads.py`on back_to_sequences results: conserve only reads whose percent of conserved kmers is at least equals to `threshold`
    `filter_reads.py file threshold`






