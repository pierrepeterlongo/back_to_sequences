import random
import string
import argparse

def generate_random_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def generate_random_fasta(num_sequences, average_size, min_size, output_file):
    sequences = []

    for _ in range(num_sequences):
        seq_size = random.randint(min_size, average_size)
        sequence = generate_random_sequence(seq_size)
        sequences.append(sequence)

    with open(output_file, 'w') as f:
        for i, sequence in enumerate(sequences, start=1):
            f.write(f'>sequence{i}\n')
            f.write(sequence + '\n')

    print(f'Generated {num_sequences} random sequences of average size {average_size} and minimum size {min_size} and saved to {output_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a random FASTA file")
    parser.add_argument("num_sequences", type=int, help="Number of sequences to generate")
    parser.add_argument("average_size", type=int, help="Average size of sequences")
    parser.add_argument("min_size", type=int, help="Minimum size of sequences")
    parser.add_argument("output_file", help="Output FASTA file")
    args = parser.parse_args()

    # check that average_size >= min_size
    if args.average_size < args.min_size:
        raise ValueError("average_size must be greater than or equal to min_size")
    

    generate_random_fasta(args.num_sequences, args.average_size, args.min_size, args.output_file)
