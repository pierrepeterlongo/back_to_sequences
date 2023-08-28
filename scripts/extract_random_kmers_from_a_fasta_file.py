import random
import argparse
import yield_reads
import revcomp


# check if a file is fasta or fastq by looking at the first line


def extract_random_kmers(input_file, kmer_size, num_kmers, canonical, output_file):
    
    reads = []
    for read in yield_reads.read_yielder(input_file):
        reads.append(read)
    nb_kmers_extracted = 0
    with open(output_file, "w") as f:
        while nb_kmers_extracted < num_kmers:
            # get a random read among the reads
            read = random.choice(reads)
            # check if the read is long enough
            if len(read) < kmer_size:
                continue
            # get a random kmer in the read
            start = random.randint(0, len(read) - kmer_size)
            kmer = read[start:start + kmer_size]
            if canonical:
                kmer = revcomp.canonical(kmer)
            f.write(f">{nb_kmers_extracted}\n{kmer}\n")
            nb_kmers_extracted += 1

    print(f"{num_kmers} random kmers extracted and saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract random kmers from a FASTQ file")
    parser.add_argument("input_file", help="Input FASTQ file (can be gzipped)")
    parser.add_argument("kmer_size", type=int, help="Size of kmers to extract")
    parser.add_argument("num_kmers", type=int, help="Number of kmers to extract")
    parser.add_argument("--canonical", action="store_true", help="Extract kmers in their canonical form")
    parser.add_argument("output_file", help="Output fasta file")
    args = parser.parse_args()

    extract_random_kmers(args.input_file, args.kmer_size, args.num_kmers, args.canonical, args.output_file)
