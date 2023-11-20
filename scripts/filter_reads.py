import argparse

def print_filtered_reads(input_file, threshold):
    selected_reads = []

    with open(input_file, 'r') as file:
        current_read = None
        while True:
            header = file.readline()
            if not header:
                # End of file
                break
            sequence = file.readline().strip()
            if float(header.split()[-1]) >= threshold:
                print(f"{header}{sequence}")
        
# main function that takes in the file path and threshold value from the user cli
# uses argparse to parse the arguments
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="file path to the fasta file")
    parser.add_argument("threshold", help="threshold value to filter the reads")
    args = parser.parse_args()
    print_filtered_reads(args.file, float(args.threshold))

if __name__ == "__main__":
    main()


