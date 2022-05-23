import sys

from src.correct_files import correct_fasta


def main():
    original_fpath, corrected_fpath = sys.argv[1:]
    correct_fasta(original_fpath, corrected_fpath=corrected_fpath, is_reference=True)


if __name__ == '__main__':
    main()