import sys

from src.correct_files import correct_fasta


def main():
    original_fpath, corrected_fpath, min_contig = sys.argv[1:]
    correct_fasta(original_fpath, int(min_contig), corrected_fpath)


if __name__ == '__main__':
    main()