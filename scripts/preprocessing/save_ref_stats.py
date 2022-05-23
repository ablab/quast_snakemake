import sys

from src.common import save_ref_stats


def main():
    fasta_fpath, csv_fpath = sys.argv[1:]
    save_ref_stats(fasta_fpath, csv_fpath)


if __name__ == '__main__':
    main()

