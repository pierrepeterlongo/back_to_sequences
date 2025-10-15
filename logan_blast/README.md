The `blast_logan_contigs.sh` script considers a query, and a list of contig files (in contigs.txt). This list may be provided by a logan search result.

For each contig file, it will run local blast between the query and the subset of contigs containing at least one shared k-mer (k=17 by default) with the query (uses `back_to_sequences`), and save the results in a directory named `<query_name>_vs_<accessions>.recruited_contigs`.

## requires 
- biopython: `pip3 install biopython`
- blast (*eg* brew install blast or look at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- back_to_sequences (cf parent directory :) )

## Running the script

```bash
./blast_logan_contigs.sh --contigs contigs.txt --query query.fa [--delete] [--kmer-size]
```
The `--delete` option will remove the intermediate files (recruited contigs and contigs files) after processing.

## Creating contigs.txt
Logan-Search results can be used to create contigs.txt. Given a logan search result link:

```bash
tail -n +2 kmviz-c21feaeb-4f33-4abc-b119-7db7bd47069b/query.tsv | awk '{print $1}' > contigs.txt
```