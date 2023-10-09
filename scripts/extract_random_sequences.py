import random
import argparse
import yield_reads
import revcomp


# check if a file is fasta or fastq by looking at the first line


def extract_random_sequences(input_file, min_size, max_size, num_kmers, output_file):
    
    reads = []
    for read in yield_reads.read_yielder(input_file):
        reads.append(read)
    nb_kmers_extracted = 0
    with open(output_file, "w") as f:
        while nb_kmers_extracted < num_kmers:
            # get a random read among the reads
            read = random.choice(reads)
            # check if the read is long enough
            if len(read) < min_size:
                continue
            # extract two random positions start and stop on the read, such that the 
            # length of the extracted sequence is between min_size and max_size
            start = random.randint(0, len(read) - min_size)
            length = random.randint(min_size, min(max_size, len(read) - start))
            sequence = read[start:start + length]
            f.write(f">{nb_kmers_extracted}\n{sequence}\n")
            nb_kmers_extracted += 1

    print(f"{num_kmers} random sequences extracted and saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract random kmers from a FASTQ file")
    parser.add_argument("--input", help="Input FASTQ file (can be gzipped)", dest="input_file", required=True, type=str)
    parser.add_argument("--min_size", help="Minimal size of sequences to extract", dest="min_size", required=True, type=int)
    parser.add_argument("--max_size", help="Minimal size of sequences to extract", dest="max_size", required=True, type=int)
    parser.add_argument("--num", help="Number of sequences to extract", dest="num_kmers", required=True, type=int)
    # parser.add_argument("--canonical", action="store_true", help="Extract kmers in their canonical form")
    parser.add_argument("--output", help="Output fasta file", dest="output_file", required=True, type=str)
    args = parser.parse_args()
    extract_random_sequences(args.input_file, args.min_size, args.max_size, args.num_kmers, args.output_file)
