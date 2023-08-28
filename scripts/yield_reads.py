import gzip
import enum

"""
yield reads from fasta or fastq files, gzipped or not
"""

class FASTYPE(enum.Enum):
    FASTA = 1
    FASTQ = 2
    OTHER = 3

def type_fastaq(f):
    first_line = f.readline()
    f.seek(0)
    if first_line[0] == ">":
        return FASTYPE.FASTA
    if first_line[0] == "@":
        return FASTYPE.FASTQ
    return FASTYPE.OTHER
    


def read_yielder(input_file_name: str):
    #input_file_name = "my_file_name[.gz]"
    if input_file_name.endswith(".gz"):
        f=gzip.open(input_file_name,'rt')
    else:
        f=open(input_file_name,'rt')
    type = type_fastaq(f)
    
    if type == FASTYPE.FASTA:
        _ = f.readline()
        while True: # concatenate lines as long as they do not start with a '>' character
            seq = ""
            read_new_header = False
            while True:
                line = f.readline()
                if not line:
                    break
                if line[0] == ">":
                    read_new_header = True
                    break
                seq += line.strip()
            yield seq
            if not read_new_header:
                break
        # for line in f:
        #     if line[0]!=">":
        #         yield line.strip()
    if type == FASTYPE.FASTQ:
        i=1
        for line in f:
            if i%4 == 2:
                yield line.strip()
            i += 1
    
    if type == FASTYPE.OTHER:
        raise TypeError(f"{input_file_name} is neither a fasta or a fastq file")
    
    f.close()    

def main():
    import sys
    for read in read_yielder(sys.argv[1]):
        print(read)

if __name__ == '__main__':
    main()
    