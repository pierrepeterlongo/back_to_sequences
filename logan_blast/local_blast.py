from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio import SeqIO
import os

def run_blast_with_alignments(query_fasta, targets_fasta, output_dir="blast_results"):
    os.makedirs(output_dir, exist_ok=True)

    # Read query sequence
    query_record = next(SeqIO.parse(query_fasta, "fasta"))
    query_id = query_record.id

    # Create a BLAST database from the targets
    makeblastdb_cline = NcbimakeblastdbCommandline(
        dbtype="nucl",
        input_file=targets_fasta,
        out=os.path.join(output_dir, "targets_db")
    )
    makeblastdb_cline()

    # Run BLAST for each target and save alignments
    for target_record in SeqIO.parse(targets_fasta, "fasta"):
        target_id = target_record.id
        output_file = os.path.join(output_dir, f"{query_id}_vs_{target_id}.txt")

        blastn_cline = NcbiblastnCommandline(
            query=query_fasta,
            db=os.path.join(output_dir, "targets_db"),
            out=output_file,
            outfmt=0,  # Pairwise alignment format
            num_alignments=1
        )
        blastn_cline()

        # Print the alignment results
        with open(output_file, "r") as f:
            alignment_results = f.read()
            print(f"Alignment for {query_id} vs {target_id}:\n{alignment_results}\n")

if __name__ == "__main__":
    run_blast_with_alignments("Q.fa", "T.fa")
