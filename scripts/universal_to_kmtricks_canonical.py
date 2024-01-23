# Author: Pierre Peterlongo
# Last modified: 23/01/2024

# a canonical kmer is the smallest value between a kmer and its reverse complement
# the order is usually defined as A < C < G < T. 
# In kmtricks, the order is A < C < T < G.
# This script computes the canonical kmer in kmtricks order.

import sys

rev = str.maketrans("ACGTacgt", "TGCAtgca")
kmtricks_order = str.maketrans("ACGTacgt", "01320132") # A < C < T < G in kmtricks
def reverse_complement(seq: str) -> str:
    return seq.translate(rev)[::-1]

def mymin(a, b):
    kmcode_a = a.translate(kmtricks_order)
    kmcode_b = b.translate(kmtricks_order)
    if kmcode_a < kmcode_b:
        return a
    else:
        return b

def kmtricks_canonical(sequence: str):
    """Returns the smallest value between a sequence and its reverse complement

    Args:
        sequence (str): sequence

    Returns:
        (str): mallest value between a sequence and its reverse complement
    """
    return mymin(reverse_complement(sequence), sequence)

def main():
    # a line is a kmer and a count.
    # the output is the canonical kmer and the count.
    with(open(sys.argv[1], 'r')) as f:
        for line in f:
            kmer, count = line.strip().split()
            print(kmtricks_canonical(kmer), count)


if __name__ == "__main__":
    main()