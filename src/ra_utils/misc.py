############################################################################
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import gzip
import os
import re
import shutil

from src import qutils, qconfig
from src.common import get_chr_lengths_from_fastafile
from src.logger import print_info, print_warning, print_error, print_notice

from collections import OrderedDict
from urllib.request import urlopen

from os.path import join, isfile, basename, dirname, getsize

from src.qutils import is_non_empty_file, get_free_memory, correct_name, calc_median, percentile

gridss_dirpath = None
gridss_version = '1.4.1'
gridss_fname = 'gridss-' + gridss_version + '.jar'


def get_safe_fpath(output_dirpath, fpath):  # reuse file if it exists; else write in output_dir
    if not isfile(fpath):
        return join(output_dirpath, basename(fpath))
    return fpath


def correct_paired_reads_names(fpath, name_ending, output_dir):
    name, ext = os.path.splitext(fpath)
    try:
        if ext in ['.gz', '.gzip']:
            handler = gzip.open(fpath, mode='rt')
            corrected_fpath = join(output_dir, basename(name))
        else:
            handler = open(fpath)
            corrected_fpath = join(output_dir, basename(fpath))
    except IOError:
        return False
    if is_non_empty_file(corrected_fpath):
        print_info('Using existing FASTQ file ' + corrected_fpath)
        return corrected_fpath
    with handler as f:
        with open(corrected_fpath, 'w') as out_f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    full_read_name = line.split()[0] + name_ending
                    out_f.write(full_read_name + '\n')
                elif i % 2 == 0:
                    out_f.write('+\n')
                else:
                    out_f.write(line)
    return corrected_fpath


def paired_reads_names_are_equal(reads_fpaths, temp_output_dir):
    first_read_names = []

    for idx, fpath in enumerate(reads_fpaths):  # Note: will work properly only with exactly two files
        name_ending = '/%d' % (idx + 1)
        reads_type = 'forward'
        if idx:
            reads_type = 'reverse'

        _, ext = os.path.splitext(fpath)
        try:
            if ext in ['.gz', '.gzip']:
                handler = gzip.open(fpath, mode='rt')
            else:
                handler = open(fpath)
        except IOError:
            print_notice('Cannot check equivalence of paired reads names, BWA may fail if reads are discordant')
            return True
        first_line = handler.readline()
        handler.close()
        full_read_name = first_line.strip().split()[0]
        if len(full_read_name) < 3 or not full_read_name.endswith(name_ending):
            print_info()
            print_notice('Improper read names in ' + fpath + ' (' + reads_type + ' reads)! '
                           'Names should end with /1 (for forward reads) or /2 (for reverse reads) '
                           'but ' + full_read_name + ' was found!\nQUAST will attempt to fix read names.')
            corrected_fpath = correct_paired_reads_names(fpath, name_ending, temp_output_dir)
            if not corrected_fpath:
                print_warning('Failed correcting read names. ')
                return False
            if reads_type == 'forward':
                qconfig.forward_reads[qconfig.forward_reads.index(fpath)] = corrected_fpath
            elif reads_type == 'reverse':
                qconfig.reverse_reads[qconfig.reverse_reads.index(fpath)] = corrected_fpath
        first_read_names.append(full_read_name[1:-2])  # truncate trailing /1 or /2 and @/> prefix (Fastq/Fasta)
    if len(first_read_names) != 2:  # should not happen actually
        print_warning('Something bad happened and we failed to check paired reads names!')
        return False
    if first_read_names[0] != first_read_names[1]:
        print_warning('Paired read names do not match! Check %s and %s!' % (first_read_names[0], first_read_names[1]))
        return False
    return True


def get_correct_names_for_chroms(output_dirpath, fasta_fpath, alignment_fpath, reads_csv):
    correct_chr_names = dict()
    fasta_chr_lengths = get_chr_lengths_from_fastafile(fasta_fpath)
    alignment_chr_lengths = OrderedDict()
    header_fpath = join(dirname(output_dirpath), basename(alignment_fpath) + '.header')
    if not isfile(alignment_fpath) and not isfile(header_fpath):
        return None
    if isfile(alignment_fpath):
        if alignment_fpath.endswith(".sam"):
            qutils.call_subprocess(['samtools', 'view', '-H', '-S', alignment_fpath],
                                   stdout=open(header_fpath, 'w'))
        else:
            qutils.call_subprocess(['samtools', 'view', '-H', alignment_fpath],
                                   stdout=open(header_fpath, 'w'))
    chr_name_pattern = 'SN:(\S+)'
    chr_len_pattern = 'LN:(\d+)'

    with open(header_fpath) as f:
        for l in f:
            if l.startswith('@SQ'):
                chr_name = re.findall(chr_name_pattern, l)[0]
                chr_len = re.findall(chr_len_pattern, l)[0]
                alignment_chr_lengths[chr_name] = int(chr_len)

    inconsistency = ''
    if len(fasta_chr_lengths) != len(alignment_chr_lengths):
        inconsistency = 'Number of chromosomes'
    else:
        for fasta_chr, sam_chr in zip(fasta_chr_lengths.keys(), alignment_chr_lengths.keys()):
            if correct_name(sam_chr) == fasta_chr[:len(sam_chr)] and alignment_chr_lengths[sam_chr] == fasta_chr_lengths[fasta_chr]:
                correct_chr_names[sam_chr] = fasta_chr
            elif alignment_chr_lengths[sam_chr] != fasta_chr_lengths[fasta_chr]:
                inconsistency = 'Chromosome lengths'
                break
            else:
                inconsistency = 'Chromosome names'
                break
    if inconsistency:
        if reads_csv:
            print_warning(inconsistency + ' in ' + fasta_fpath + ' and corresponding file ' + alignment_fpath + ' do not match. ' +
                           'QUAST will try to realign reads')
        else:
            print_error(inconsistency + ' in ' + fasta_fpath + ' and corresponding file ' + alignment_fpath + ' do not match. ' +
                         'Use SAM file obtained by aligning reads')
        return None
    return correct_chr_names


def all_read_names_correct(sam_fpath):
    with open(sam_fpath) as sam_in:
        for i, l in enumerate(sam_in):
            if i > 1000000:
                return True
            if not l:
                continue
            fs = l.split('\t')
            read_name = fs[0]
            if read_name[-2:] == '/1' or read_name[-2:] == '/2':
                return False
    return True


def clean_read_names(sam_fpath, correct_sam_fpath):
    with open(sam_fpath) as sam_in:
        with open(correct_sam_fpath, 'w') as sam_out:
            for l in sam_in:
                if not l:
                    continue
                fs = l.split('\t')
                read_name = fs[0]
                if read_name[-2:] == '/1' or read_name[-2:] == '/2':
                    fs[0] = read_name[:-2]
                    l = '\t'.join(fs)
                sam_out.write(l)
    return correct_sam_fpath


def sort_bam(bam_fpath, sorted_bam_fpath, threads=None, sort_rule=None):
    if not threads:
        threads = qconfig.max_threads
    cmd = ['samtools', 'sort', '-@', str(threads), '-o', sorted_bam_fpath, bam_fpath]
    if sort_rule:
        cmd += [sort_rule]
    qutils.call_subprocess(cmd)


def bwa_index(ref_fpath):
    cmd = ['bwa', 'index', '-p', ref_fpath, ref_fpath]
    if getsize(ref_fpath) > 2 * 1024 ** 3:  # if reference size bigger than 2GB
        cmd += ['-a', 'bwtsw']
    if not is_non_empty_file(ref_fpath + '.bwt'):
        qutils.call_subprocess(cmd)


def get_gridss_memory():
    free_mem = get_free_memory()
    if free_mem >= 64:
        return 31
    elif free_mem >= 32:
        return 16
    elif free_mem >= 16:
        return 8
    return 2


def reformat_bedpe(raw_bed_fpath, bed_fpath):
    header = None
    with open(raw_bed_fpath) as f:
        with open(bed_fpath, 'w') as out_f:
            out_f.write('\t'.join(['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B', 'TYPE\n']))
            for i, line in enumerate(f):
                if i == 0:
                    header = line[1:].split('\t')
                    continue
                if not line.startswith('#'):
                    fs = line.split('\t')
                    try:
                        # chrom1 start1  end1    chrom2  start2  end2    name    score   strand1 strand2
                        sv = dict(zip(header, fs))
                        sv_type = 'BND'
                        if sv['strand1'] == sv['strand2']:
                            sv_type = 'INV'
                        out_f.write('\t'.join([sv['chrom1'], sv['start1'], sv['end1'],
                                               sv['chrom2'], sv['start2'], sv['end2'], sv_type]) + '\n')
                    except ValueError:
                        pass


def check_cov_file(cov_fpath):
    raw_cov_fpath = cov_fpath + '_raw'
    with open(cov_fpath, 'r') as coverage:
        for line in coverage:
            if len(line.split()) != 2:
                shutil.copy(cov_fpath, raw_cov_fpath)
                os.remove(cov_fpath)
                return False
            else:
                return True


def bam_to_bed(output_dirpath, name, bam_fpath, bedpe=False):
    raw_bed_fpath = join(output_dirpath, name + '.bed')
    if bedpe:
        bedpe_fpath = join(output_dirpath, name + '.bedpe')
        if not is_non_empty_file(bedpe_fpath) and not is_non_empty_file(bedpe_fpath):
            qutils.call_subprocess(['bedtools', 'bamtobed', '-i', bam_fpath, '-bedpe'],
                                   stdout=open(bedpe_fpath, 'w'))
            with open(bedpe_fpath, 'r') as bedpe:
                with open(raw_bed_fpath, 'w') as bed_file:
                    for line in bedpe:
                        fs = line.split()
                        start, end = fs[1], fs[5]
                        bed_file.write('\t'.join([fs[0], start, end + '\n']))
    else:
        if not is_non_empty_file(raw_bed_fpath):
            qutils.call_subprocess(['bedtools', 'bamtobed', '-i', bam_fpath],
                               stdout=open(raw_bed_fpath, 'w'))

    sorted_bed_fpath = join(output_dirpath, name + '.sorted.bed')
    if not is_non_empty_file(sorted_bed_fpath):
        qutils.call_subprocess(['sort', '-k1,1', '-k2,2n', raw_bed_fpath],
                            stdout=open(sorted_bed_fpath, 'w'))
    return sorted_bed_fpath


def calculate_read_len(sam_fpath):
    read_lengths = []
    mapped_flags = ['99', '147', '83', '163']  # reads mapped in correct order
    with open(sam_fpath) as sam_in:
        for i, l in enumerate(sam_in):
            if i > 1000000:
                break
            if l.startswith('@'):
                continue
            fs = l.split('\t')
            flag = fs[1]
            if flag not in mapped_flags:
                continue
            read_lengths.append(len(fs[9]))
    return sum(read_lengths) * 1.0 / len(read_lengths)


def calculate_genome_cov(in_fpath, out_fpath, chr_len_fpath, print_all_positions=True):
    cmd = ['bedtools', 'genomecov', '-ibam' if in_fpath.endswith('.bam') else '-i', in_fpath, '-g', chr_len_fpath]
    if print_all_positions:
        cmd += ['-bga']
    qutils.call_subprocess(cmd, stdout=open(out_fpath, 'w'))


def get_max_min_is(insert_sizes):
    decile_1 = percentile(insert_sizes, 10)
    decile_9 = percentile(insert_sizes, 90)
    return decile_1, decile_9


def convert_sam(in_fpath, out_fpath, max_threads, filter_rule=None):
    cmd = ['samtools', 'view', '-@', str(max_threads), '-h']
    if in_fpath.endswith('.sam'):
        cmd += ['-S']
    if out_fpath.endswith('.bam'):
        cmd += ['-b']
    if filter_rule:
        cmd += filter_rule.split()
    cmd.append(in_fpath)
    qutils.call_subprocess(cmd, stdout=open(out_fpath, 'w'))


def is_valid_bed(bed_fpath):
    # check last 10 lines
    with open(bed_fpath, 'r') as f:
        lines_found = []
        block_counter = -1
        _buffer = 1024
        while len(lines_found) < 10:
            try:
                f.seek(block_counter * _buffer, os.SEEK_END)
            except IOError:
                f.seek(0)
                lines_found = f.readlines()
                break
            lines_found = f.readlines()
        block_counter -= 1
    for line in lines_found:
        if not line.startswith('#'):
            fs = line.split('\t')
            try:
                align1 = (int(fs[1]), int(fs[2]), correct_name(fs[0]), fs[6])
                align2 = (int(fs[4]), int(fs[5]), correct_name(fs[3]), fs[6])
            except IndexError:
                return False
    return True


def calculate_insert_size(sam_fpath, output_dir, ref_name, reads_suffix=''):
    insert_size_fpath = join(output_dir, ref_name + ('.' + reads_suffix if reads_suffix else '') + '.is.txt')
    if is_non_empty_file(insert_size_fpath):
        try:
            with open(insert_size_fpath) as f:
                insert_size = int(f.readline())
                min_insert_size = int(f.readline())
                max_insert_size = int(f.readline())
            if insert_size:
                return insert_size, min_insert_size, max_insert_size
        except:
            pass
    insert_sizes = []
    mapped_flags = ['99', '147', '83', '163']  # reads mapped in correct orientation and within insert size
    with open(sam_fpath) as sam_in:
        for i, l in enumerate(sam_in):
            if i > 1000000:
                break
            if l.startswith('@'):
                continue
            fs = l.split('\t')
            flag = fs[1]
            if flag not in mapped_flags:
                continue
            insert_size = abs(int(fs[8]))
            insert_sizes.append(insert_size)

    if insert_sizes:
        insert_sizes.sort()
        median_is = calc_median(insert_sizes)
        if median_is <= 0:
            return None, None, None
        min_insert_size, max_insert_size = get_max_min_is(insert_sizes)
        insert_size = max(qconfig.optimal_assembly_min_IS, median_is)
        with open(insert_size_fpath, 'w') as out_f:
            out_f.write(str(insert_size) + '\n')
            out_f.write(str(min_insert_size) + '\n')
            out_f.write(str(max_insert_size) + '\n')
        return insert_size, min_insert_size, max_insert_size
    return None, None, None
