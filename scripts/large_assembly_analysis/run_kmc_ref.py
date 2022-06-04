import os
import sys
from os.path import join, isdir

from src import qconfig
from src.kmers_analyzer import *
from src.logger import *


def main():
    output_dirpath, tmp_dirpath, reference, threads = sys.argv[1:]
    threads = int(threads)

    print_timestamp()
    kmer_len = qconfig.unique_kmer_len
    print_info('Running analysis based on unique ' + str(kmer_len) + '-mers...')

    print_info('  Running KMC on reference...')

    ref_kmc_out_fpath = count_kmers(tmp_dirpath, reference, kmer_len, threads)
    unique_kmers = get_kmers_cnt(tmp_dirpath, ref_kmc_out_fpath, threads)
    if not unique_kmers:
        print_warning('KMC failed. Skipping...')
        return

    print_info('    Downsampling k-mers...')
    downsample_kmers(tmp_dirpath, reference, ref_kmc_out_fpath, kmer_len, threads)


if __name__ == '__main__':
    main()

