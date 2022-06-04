############################################################################
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
import argparse
import os
import shutil
from collections import defaultdict
from os.path import join, abspath, exists, basename, isdir

from src import qconfig, qutils, reporting
from src.common import read_fasta
from src.logger import *
from src.qutils import md5, get_free_memory
from src.reporting import save_kmers

KMER_FRACTION = 0.001
KMERS_INTERVAL = 1000
MAX_CONTIGS_NUM = 10000
MAX_REF_CONTIGS_NUM = 200
MIN_CONTIGS_LEN = 10000
EXT_RELOCATION_SIZE = 100000


def create_kmc_stats_file(output_dirpath, contigs_fpath, label, completeness,
                          corr_len, mis_len, undef_len, total_len, translocations, relocations):
    kmc_check_fpath = join(output_dirpath, label + '.sf')
    kmc_stats_fpath = join(output_dirpath, label + '.stat')
    with open(kmc_check_fpath, 'w') as check_f:
        check_f.write("Assembly md5 checksum: %s\n" % md5(contigs_fpath))
    with open(kmc_stats_fpath, 'w') as stats_f:
        stats_f.write("Completeness: %s\n" % completeness)
        if corr_len or mis_len:
            stats_f.write("K-mer-based correct length: %d\n" % corr_len)
            stats_f.write("K-mer-based misjoined length: %d\n" % mis_len)
            stats_f.write("K-mer-based undefined length: %d\n" % undef_len)
            stats_f.write("Total length: %d\n" % total_len)
            stats_f.write("# translocations: %d\n" % translocations)
            stats_f.write("# 100 kbp relocations: %d\n" % relocations)


def check_kmc_successful_check(output_dirpath, contigs_fpath, label, ref_fpath):
    kmc_check_fpath = join(output_dirpath, label + '.sf')
    if not exists(kmc_check_fpath):
        return False
    successful_check_content = open(kmc_check_fpath).read().split('\n')
    if len(successful_check_content) < 2:
        return False
    if successful_check_content[0].strip().split()[-1] != str(md5(contigs_fpath)):
        return False
    if successful_check_content[1].strip().split()[-1] != str(md5(ref_fpath)):
        return False
    return True


def get_kmers_cnt(tmp_dirpath, kmc_db_fpath, threads, run_histo=True):
    histo_fpath = join(tmp_dirpath, basename(kmc_db_fpath) + '.histo.txt')
    if run_histo:
        run_kmc(['transform', kmc_db_fpath, 'histogram', histo_fpath], threads)
    kmers_cnt = 0
    if exists(histo_fpath):
        kmers_cnt = int(open(histo_fpath).read().split()[-1])
    return kmers_cnt


def count_kmers(tmp_dirpath, fpath, kmer_len, threads):
    kmc_out_fpath = join(tmp_dirpath, 'reference.kmc')
    max_mem = max(2, get_free_memory())
    run_kmc(['-m' + str(max_mem), '-n128', '-k' + str(kmer_len), '-fm', '-cx1', '-ci1', fpath, kmc_out_fpath, tmp_dirpath],
             threads, use_kmc_tools=False)
    return kmc_out_fpath


def align_kmers(label, output_dirpath, ref_fpath, kmers_fpath, max_threads):
    out_fpath = join(output_dirpath, label + '.kmers.coords')
    cmdline = ['minimap2', '-cx', 'sr', '-s' + str(qconfig.unique_kmer_len * 2), '--frag=no',
               '-t', str(max_threads), ref_fpath, kmers_fpath]
    qutils.call_subprocess(cmdline, indent='  ', stdout=open(out_fpath, "w"))
    kmers_pos_by_chrom = defaultdict(list)
    kmers_by_chrom = defaultdict(list)
    with open(out_fpath) as f:
        for line in f:
            fs = line.split('\t')
            if len(fs) < 10:
                continue
            contig, chrom, pos = fs[0], fs[5], fs[7]
            kmers_pos_by_chrom[chrom].append(int(pos))
            kmers_by_chrom[chrom].append(contig)
    return kmers_by_chrom, kmers_pos_by_chrom


def downsample_kmers(tmp_dirpath, ref_fpath, kmc_db_fpath, kmer_len, threads):
    downsampled_txt_fpath = join(tmp_dirpath, 'kmc.downsampled.txt')
    open(downsampled_txt_fpath, 'w').close()
    prev_kmer_idx = 0
    for record in read_fasta(ref_fpath):
        chrom, seq = record.id, str(record.seq)
        kmc_fasta_fpath = join(tmp_dirpath, 'kmers_' + chrom + '.fasta')
        num_kmers_in_seq = len(seq) - kmer_len + 1
        with open(kmc_fasta_fpath, 'w') as out_f:
            for i in range(num_kmers_in_seq):
                out_f.write('>' + str(i) + '\n')
                out_f.write(seq[i: i + kmer_len] + '\n')
        filtered_fpath = join(tmp_dirpath, 'kmers_' + chrom + '.filtered.fasta')
        filter_contigs(kmc_fasta_fpath, filtered_fpath, kmc_db_fpath, threads, min_kmers=1)
        filtered_kmers = set()
        for record in read_fasta(filtered_fpath):
            filtered_kmers.add(record.id)
        with open(downsampled_txt_fpath, 'a') as out_f:
            kmer_i = 0
            for record in read_fasta(kmc_fasta_fpath):
                if record.id in filtered_kmers:
                    if not kmer_i or int(record.id) - kmer_i >= KMERS_INTERVAL:
                        kmer_i = int(record.id)
                        out_f.write('>%s_%d\n' % (chrom, kmer_i))
                        out_f.write(str(record.seq) + '\n')
        prev_kmer_idx += num_kmers_in_seq
        if qconfig.space_efficient:
            os.remove(kmc_fasta_fpath)
    return downsampled_txt_fpath


def get_clear_name(fpath):
    return basename(fpath).replace('.kmc', '')


def intersect_kmers(tmp_dirpath, kmc_out_fpaths, threads):
    intersect_out_fpath = join(tmp_dirpath, '_'.join([get_clear_name(kmc_out_fpath)[:30] for kmc_out_fpath in kmc_out_fpaths]) + '.kmc')
    if len(kmc_out_fpaths) == 2:
        run_kmc(['simple'] + kmc_out_fpaths + ['intersect', intersect_out_fpath], threads)
    else:
        prev_kmc_out_fpath = kmc_out_fpaths[0]
        for i in range(1, len(kmc_out_fpaths)):
            tmp_out_fpath = join(tmp_dirpath, get_clear_name(prev_kmc_out_fpath) + '_' + str(i) + '.kmc')
            run_kmc(['simple', prev_kmc_out_fpath, kmc_out_fpaths[i], 'intersect', tmp_out_fpath], threads)
            prev_kmc_out_fpath = tmp_out_fpath
        intersect_out_fpath = prev_kmc_out_fpath
    return intersect_out_fpath


def filter_contigs(input_fpath, output_fpath, db_fpath, threads, min_kmers=1):
    if input_fpath.endswith('.txt'):
        input_fpath = '@' + input_fpath
    run_kmc(['filter', db_fpath, input_fpath, '-ci' + str(min_kmers), '-fa', output_fpath], threads)


def parse_kmer_results(reports, output_dirpath, labels):
    for index, label in enumerate(labels):
        report = reports[label]
        kmc_stats_fpath = join(output_dirpath, label + '.stat')
        stats_content = open(kmc_stats_fpath).read().split('\n')
        if len(stats_content) < 1:
            return
        report.add_field(reporting.Fields.KMER_COMPLETENESS, '%.2f' % float(stats_content[0].strip().split(': ')[-1]))
        if len(stats_content) >= 7:
            corr_len = int(stats_content[1].strip().split(': ')[-1])
            mis_len = int(stats_content[2].strip().split(': ')[-1])
            undef_len = int(stats_content[3].strip().split(': ')[-1])
            total_len = int(stats_content[4].strip().split(': ')[-1])
            translocations = int(stats_content[5].strip().split(': ')[-1])
            relocations = int(stats_content[6].strip().split(': ')[-1])
            report.add_field(reporting.Fields.KMER_CORR_LENGTH, '%.2f' % (corr_len * 100.0 / total_len))
            report.add_field(reporting.Fields.KMER_MIS_LENGTH, '%.2f' % (mis_len * 100.0 / total_len))
            report.add_field(reporting.Fields.KMER_UNDEF_LENGTH, '%.2f' % (undef_len * 100.0 / total_len))
            report.add_field(reporting.Fields.KMER_TRANSLOCATIONS, translocations)
            report.add_field(reporting.Fields.KMER_RELOCATIONS, relocations)
            report.add_field(reporting.Fields.KMER_MISASSEMBLIES, translocations + relocations)
    save_kmers(reports, output_dirpath)


def run_kmc(params, threads, use_kmc_tools=True):
    tool_fpath = 'kmc_tools' if use_kmc_tools else 'kmc'
    qutils.call_subprocess([tool_fpath, '-t' + str(threads), '-hp'] + params)


def _get_dist_inconstistency(pos, prev_pos, ref_pos, prev_ref_pos, cyclic_ref_lens):
    dist = abs(abs(pos - prev_pos) - abs(ref_pos - prev_ref_pos))
    if cyclic_ref_lens and ref_pos < prev_ref_pos and cyclic_ref_lens - dist < dist:
        dist = abs(abs(pos - prev_pos) - abs(ref_pos - prev_ref_pos + cyclic_ref_lens))
    return dist

