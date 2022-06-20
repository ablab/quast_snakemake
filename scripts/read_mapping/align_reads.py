import argparse
import os
import shlex
import sys
from os.path import join, isfile, dirname, abspath, exists

from src import qutils, qconfig
from src.common import get_labels_from_paths
from src.logger import print_info
from src.qutils import is_non_empty_file, add_suffix, assert_file_exists, get_chr_len_fpath
from src.ra_utils.misc import *


def align_single_file(ref_fpath, label, reads_fpaths, reads_types, output_dirpath, max_threads, reads_type='all'):
    # use absolute paths because we will change workdir
    fpath = abspath(ref_fpath)
    output_dirpath = abspath(output_dirpath)
    sam_fpath = join(output_dirpath, label + '.sam')
    bam_fpath = join(output_dirpath, label + '.bam')

    for i in range(len(reads_fpaths)):
        reads = reads_fpaths[i].replace(',', ' ').split()
        reads_fpaths[i] = ' '.join(abspath(r) for r in reads)

    print_info('  Running BWA...')

    prev_dir = os.getcwd()
    os.chdir(output_dirpath)
    bwa_index(fpath)
    sam_fpaths = align_reads(fpath, sam_fpath, reads_fpaths, reads_types, output_dirpath, max_threads)

    if len(sam_fpaths) > 1:
        merge_sam_files(sam_fpaths, sam_fpath, bam_fpath, max_threads)
    elif len(sam_fpaths) == 1:
        shutil.move(sam_fpaths[0], sam_fpath)
        tmp_bam_fpath = sam_fpaths[0].replace('.sam', '.bam')
        if is_non_empty_file(tmp_bam_fpath):
            shutil.move(tmp_bam_fpath, bam_fpath)

    print_info('  Done.')
    os.chdir(prev_dir)
    if not is_non_empty_file(sam_fpath):
        print_error('  Failed running BWA for ' + fpath)
        return None, None, None

    correct_sam_fpath = join(output_dirpath, label + '.' + reads_type + '.correct.sam')  # write in output dir
    clean_read_names(sam_fpath, correct_sam_fpath)
    convert_sam(correct_sam_fpath, bam_fpath, max_threads)

    assert_file_exists(bam_fpath, 'bam file')


def align_reads(ref_fpath, sam_fpath, reads_fpaths, reads_types, output_dir, max_threads):
    out_sam_fpaths = []
    for idx, (reads_fpath, reads_type) in enumerate(zip(reads_fpaths, reads_types)):
        reads = reads_fpath.split()
        if all(exists(r) for r in reads):
            run_aligner(reads_fpath, ref_fpath, sam_fpath, out_sam_fpaths, output_dir, max_threads, idx, reads_type=reads_type)
    return out_sam_fpaths


def run_aligner(reads, ref_fpath, sam_fpath, out_sam_fpaths, output_dir, max_threads, idx, reads_type):
    bwa_cmd = 'bwa mem -t ' + str(max_threads)
    insert_sizes = []
    temp_sam_fpaths = []
    if reads_type == 'pacbio' or reads_type == 'nanopore':
        if reads_type == 'pacbio':
            preset = ' -ax map-pb '
        else:
            preset = ' -ax map-ont '
        cmdline = 'minimap2' + ' -t ' + str(max_threads) + preset + ref_fpath + ' ' + reads
    else:
        cmdline = bwa_cmd + (' -p ' if reads_type == 'interleaved' else ' ') + ref_fpath + ' ' + reads
    output_fpath = add_suffix(sam_fpath, reads_type + str(idx + 1))
    bam_fpath = output_fpath.replace('.sam', '.bam')
    qutils.call_subprocess(shlex.split(cmdline), stdout=open(output_fpath, 'w'))
    convert_sam(output_fpath, bam_fpath, max_threads)
    if reads_type == 'pe':
        insert_size, _, _ = calculate_insert_size(output_fpath, output_dir, qutils.name_from_fpath(sam_fpath))
        if insert_size is not None and insert_size < qconfig.optimal_assembly_max_IS:
            insert_sizes.append(insert_size)
    temp_sam_fpaths.append(output_fpath)

    if len(temp_sam_fpaths) == 1:
        final_sam_fpath = add_suffix(sam_fpath, reads_type)
        final_bam_fpath = final_sam_fpath.replace('.sam', '.bam')
        shutil.move(temp_sam_fpaths[0], final_sam_fpath)
        shutil.move(temp_sam_fpaths[0].replace('.sam', '.bam'), final_bam_fpath)
        out_sam_fpaths.append(final_sam_fpath)
    else:
        out_sam_fpaths.extend(temp_sam_fpaths)

    if insert_sizes:
        ref_name = qutils.name_from_fpath(ref_fpath)
        insert_size_fpath = join(output_dir, ref_name + '.is.txt')
        with open(insert_size_fpath, 'w') as out:
            out.write(str(max(insert_sizes)))


def merge_sam_files(tmp_sam_fpaths, sam_fpath, bam_fpath, max_threads):
    tmp_bam_fpaths = []
    for tmp_sam_fpath in tmp_sam_fpaths:
        if is_non_empty_file(tmp_sam_fpath):
            tmp_bam_fpath = tmp_sam_fpath.replace('.sam', '.bam')
            tmp_bam_sorted_fpath = add_suffix(tmp_bam_fpath, 'sorted')
            if not is_non_empty_file(tmp_bam_sorted_fpath):
                sort_bam(tmp_bam_fpath, tmp_bam_sorted_fpath)
            tmp_bam_fpaths.append(tmp_bam_sorted_fpath)
    qutils.call_subprocess(['samtools', 'merge', '-@', str(max_threads), bam_fpath] + tmp_bam_fpaths)
    convert_sam(bam_fpath, sam_fpath, max_threads)
    return sam_fpath


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads_fpaths', nargs='+')
    parser.add_argument('--reads_types', nargs='+')
    parser.add_argument('--threads', type=int)
    parser.add_argument('--reference')
    parser.add_argument('--output_dirpath')

    args = parser.parse_args()

    label = get_labels_from_paths([args.reference])[0]
    align_single_file(args.reference, label, args.reads_fpaths, args.reads_types, args.output_dirpath, args.threads, reads_type='all')


if __name__ == '__main__':
    main()