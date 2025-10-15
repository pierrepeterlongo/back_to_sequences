The `blast_logan.sh` script considers a query, and a file continaing list of SRA accessions (in accessions.txt). This list may be provided by a Logan-Search result.

For each accession, it will run local blast between the query and the subset of contigs or unitigs containing at least one shared k-mer (k=17 by default) with the query (uses `back_to_sequences`), and save the results in a directory named `<accessions>.recruited_[contigs/unitigs]_vs_query`.

## requires 
- biopython: `pip3 install biopython`
- blast (*eg* brew install blast or look at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- back_to_sequences (cf parent directory :) )

## Running the script

```bash
sh ./logan_blast.sh 

  -a, --accessions  Path to accessions.txt file. Containing one accession per line)
  -q, --query       Path to query fasta file
  -u, --unitigs     Consider the unitigs verison of the accessions instead of contigs
  -k, --kmer-size K-mer size for sequence recruitment with back_to_sequences (default: 17)
  -l, --limit     Limit number of accessions to process from accession file (default: no limit)
  -d, --delete    Delete recruited accessions and accessions files after processing (default: keep all files)
  -h, --help      Show this help message
```

## Creating accessions.txt
Logan-Search results can be used to create accessions.txt. Given a logan search result link or after exporting a table from Logan_search interface:

```bash
tail -n +2 kmviz-c21feaeb-4f33-4abc-b119-7db7bd47069b/query.tsv | awk '{print $1}' | tr -d "\"" > accessions.txt
```

## Example command

```bash
sh ./logan_blast.sh  -a example/accessions.txt -q example/query.fa -l 10 -k 27 -u
```