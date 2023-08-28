rev = str.maketrans("ACGTacgt", "TGCAtgca")
def reverse_complement(seq: str) -> str:
    return seq.translate(rev)[::-1]



def canonical(sequence: str):
    """Returns the smallest value between a sequence and its reverse complement

    Args:
        sequence (str): sequence

    Returns:
        (str): mallest value between a sequence and its reverse complement
    """
    return min(reverse_complement(sequence), sequence)

