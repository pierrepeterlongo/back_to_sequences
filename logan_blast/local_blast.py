import subprocess
from Bio import SeqIO
import os
import argparse
import shutil



def run_blast_with_alignments(query_fasta, targets_fasta, output_dir="blast_results"):
    os.makedirs(output_dir, exist_ok=True)

    # Read query sequence
    query_record = next(SeqIO.parse(query_fasta, "fasta"))
    query_id = query_record.id

    # Create a BLAST database from the targets
    makeblastdb_cmd = [
        "makeblastdb",
        "-dbtype", "nucl",
        "-in", targets_fasta,
        "-out", os.path.join(output_dir, "targets_db")
    ]
    result = subprocess.run(makeblastdb_cmd, check=True, capture_output=True, text=True)

    # Optionally, print or log the output/errors
    # if result.stdout:
    #     print("STDOUT:", result.stdout)
    # if result.stderr:
    #     print("STDERR:", result.stderr)

    # Run BLAST for each target and save alignments
    for target_record in SeqIO.parse(targets_fasta, "fasta"):
        target_id = target_record.id
        output_file = os.path.join(output_dir, f"{query_id}_vs_{target_id}.txt")

        # Define your command as a list of arguments
        blastn_cmd = [
            "blastn",
            "-query", query_fasta,
            "-db", os.path.join(output_dir, "targets_db"),
            "-out", output_file,
            "-outfmt", "0",  # Pairwise alignment format
            "-num_alignments", "1",
            "-sorthits", "0",
        ]

        # Run the command
        result = subprocess.run(blastn_cmd, check=True, capture_output=True, text=True)

        # Optionally, print or log the output/errors
        # if result.stdout:
        #     print("STDOUT:", result.stdout)
        # if result.stderr:
        #     print("STDERR:", result.stderr)

       

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BLAST alignments between query and target FASTA files.")
    parser.add_argument("-q", "--query", required=True, help="Query FASTA file")
    parser.add_argument("-t", "--target", required=True, help="Target FASTA file")
    args = parser.parse_args()
    
    if shutil.which("makeblastdb") is None or shutil.which("blastn") is None:
        raise EnvironmentError("BLAST+ executables 'makeblastdb' and/or 'blastn' not found in PATH. Please install BLAST+ locally.")
    
    query_basename = os.path.splitext(os.path.basename(args.query))[0]
    target_basename = os.path.splitext(os.path.basename(args.target))[0]
    run_blast_with_alignments(args.query, args.target, f"{query_basename}_vs_{target_basename}")
    print(f"\033[92mResults are in directory {query_basename}_vs_{target_basename}\033[0m")
