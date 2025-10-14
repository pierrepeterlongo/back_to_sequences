## requires
pip3 install biopython
need blast (eg brew install blast)

## Process to be extended
Given a contig or unitig list: contigs.txt and a query Q.fa

```bash
while read accession; do
	aws s3 cp s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst . --no-sign-request
	back_to_sequences --in-kmers Q.fa --in-sequences  ${accession}.contigs.fa.zst --out-sequences T.fa
	python local_blast.py 
done < contigs.txt
```

