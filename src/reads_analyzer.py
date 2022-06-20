############################################################################
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
from __future__ import division
import os
import re
import shutil
import shlex
from collections import defaultdict
from math import sqrt
from os.path import isfile, join, basename, abspath, isdir, dirname, exists

import pysam as pysam

from src import qconfig, reporting
from src.align_reads import align_single_file
from src.qutils import *
from src.ra_utils.misc import *
from src.reporting import save_reads

ref_sam_fpaths = {}
COVERAGE_FACTOR = 10


class Mapping(object):
    MIN_MAP_QUALITY = 20  # for distiguishing "good" reads and "bad" ones

    def __init__(self, ref, start, mapq, ref_next, query_len):
        #self.ref, self.start, self.mapq, self.ref_next, self.len = \
        #    fields[2], int(fields[3]), int(fields[4]), fields[6], len(fields[9])
        self.ref, self.start, self.mapq, self.ref_next, self.len = ref, start, mapq, ref_next, query_len
        self.end = self.start + self.len - 1  # actually not always true because of indels
        self.mapq = int(mapq)

    @staticmethod
    def parse(line):
        if line.startswith('@'):  # comment
            return None
        if len(line.split('\t')) < 11:  # not valid line
            return None
        mapping = Mapping(line.split('\t'))
        return mapping


class QuastDeletion(object):
    ''' describes situtations: GGGGBBBBBNNNNNNNNNNNNBBBBBBGGGGGG, where
    G -- "good" read (high mapping quality)
    B -- "bad" read (low mapping quality)
    N -- no mapped reads
    size of Ns fragment -- "deletion" (not less than MIN_GAP)
    size of Bs fragment -- confidence interval (not more than MAX_CONFIDENCE_INTERVAL,
        fixing last/first G position otherwise)
    '''

    MAX_CONFIDENCE_INTERVAL = 150
    MIN_GAP = qconfig.extensive_misassembly_threshold - 2 * MAX_CONFIDENCE_INTERVAL

    def __init__(self, ref, prev_good=None, prev_bad=None, next_bad=None, next_good=None, next_bad_end=None):
        self.ref, self.prev_good, self.prev_bad, self.next_bad, self.next_good, self.next_bad_end = \
            ref, prev_good, prev_bad, next_bad, next_good, next_bad_end
        self.id = 'QuastDEL'

    def is_valid(self):
        return self.prev_good is not None and self.prev_bad is not None and \
               self.next_bad is not None and self.next_good is not None and \
               (self.next_bad - self.prev_bad > QuastDeletion.MIN_GAP)

    def set_prev_good(self, mapping):
        self.prev_good = mapping.end
        self.prev_bad = self.prev_good  # prev_bad cannot be earlier than prev_good!
        return self  # to use this function like "deletion = QuastDeletion(ref).set_prev_good(coord)"

    def set_prev_bad(self, mapping=None, position=None):
        self.prev_bad = position if position else mapping.end
        if self.prev_good is None or self.prev_good + QuastDeletion.MAX_CONFIDENCE_INTERVAL < self.prev_bad:
            self.prev_good = max(1, self.prev_bad - QuastDeletion.MAX_CONFIDENCE_INTERVAL)
        return self  # to use this function like "deletion = QuastDeletion(ref).set_prev_bad(coord)"

    def set_next_good(self, mapping=None, position=None):
        self.next_good = position if position else mapping.start
        if self.next_bad is None:
            self.next_bad = self.next_good
        elif self.next_good - QuastDeletion.MAX_CONFIDENCE_INTERVAL > self.next_bad:
            self.next_good = self.next_bad + QuastDeletion.MAX_CONFIDENCE_INTERVAL

    def set_next_bad(self, mapping):
        self.next_bad = mapping.start
        self.next_bad_end = mapping.end
        self.next_good = self.next_bad  # next_good is always None at this moment (deletion is complete otherwise)

    def set_next_bad_end(self, mapping):
        if self.next_bad is None:
            self.next_bad = mapping.start
        self.next_bad_end = mapping.end
        self.next_good = min(mapping.start, self.next_bad + QuastDeletion.MAX_CONFIDENCE_INTERVAL)

    def __str__(self):
        return '\t'.join(map(str, [self.ref, self.prev_good, self.prev_bad,
                          self.ref, self.next_bad, self.next_good,
                          self.id]))


def create_fai_file(cur_ref_fpath):
    qutils.call_subprocess(['samtools', 'faidx', cur_ref_fpath])


def process_one_ref(cur_ref_fpath, output_dirpath, max_threads, bam_fpath=None, bed_fpath=None):
    ref_name = qutils.name_from_fpath(cur_ref_fpath)
    if not bam_fpath:
        sam_fpath = join(output_dirpath, ref_name + '.sam')
        bam_fpath = join(output_dirpath, ref_name + '.bam')
        bam_sorted_fpath = join(output_dirpath, ref_name + '.sorted.bam')
    else:
        sam_fpath = bam_fpath.replace('.bam', '.sam')
        bam_sorted_fpath = add_suffix(bam_fpath, 'sorted')
    bed_fpath = bed_fpath or join(output_dirpath, ref_name + '.bed')
    if is_non_empty_file(bed_fpath):
        if not is_valid_bed(bed_fpath):
            print_warning('  Existing BED-file: ' + bed_fpath + ' may be corrupted. Bed file will be re-created. ')
            bed_fpath = join(output_dirpath, ref_name + '.bed')
        else:
            print_info('  Using existing BED-file: ' + bed_fpath)
            return bed_fpath

    if not isfile(bam_sorted_fpath):
        convert_sam(sam_fpath, bam_fpath, max_threads, filter_rule='-f 2') # not unmapped and proper_pair
        sort_bam(bam_fpath, bam_sorted_fpath, threads=max_threads)
    if not is_non_empty_file(bam_sorted_fpath + '.bai'):
        call_subprocess(['samtools', 'index', bam_sorted_fpath])
    create_fai_file(cur_ref_fpath)
    vcf_output_dirpath = join(output_dirpath, ref_name + '_gridss')
    vcf_fpath = join(vcf_output_dirpath, ref_name + '.vcf')
    if not is_non_empty_file(vcf_fpath):
        if isdir(vcf_output_dirpath):
            shutil.rmtree(vcf_output_dirpath, ignore_errors=True)
        os.makedirs(vcf_output_dirpath)
        max_mem = get_gridss_memory()
        bwa_index(cur_ref_fpath)
        qutils.call_subprocess(['gridss', '--reference', cur_ref_fpath, '--output', vcf_fpath, '--assembly',
                                join(vcf_output_dirpath, ref_name + '.gridss.bam'),
                                '--threads', str(max_threads), '--workingdir', vcf_output_dirpath, bam_sorted_fpath])
    if is_non_empty_file(vcf_fpath):
        raw_bed_fpath = add_suffix(bed_fpath, 'raw')
        filtered_bed_fpath = add_suffix(bed_fpath, 'filtered')
        reformat_bedpe(raw_bed_fpath, bed_fpath)
    return bed_fpath


def search_sv_with_gridss(main_ref_fpath, bam_fpath, output_dirpath, max_threads):
    print_info('  Searching structural variations with GRIDSS...')
    final_bed_fpath = join(output_dirpath, qutils.name_from_fpath(main_ref_fpath) + '_' + qconfig.sv_bed_fname)
    if isfile(final_bed_fpath):
        print_info('    Using existing file: ' + final_bed_fpath)
        return final_bed_fpath

    if not get_path_to_program('java') or not check_java_version(1.8):
        print_warning('Java 1.8 (Java version 8) or later is required to run GRIDSS. Please install it and rerun QUAST.')
        return None
    if not get_path_to_program('Rscript'):
        print_warning('R is required to run GRIDSS. Please install it and rerun QUAST.')
        return None

    process_one_ref(main_ref_fpath, output_dirpath, max_threads, bam_fpath=bam_fpath, bed_fpath=final_bed_fpath)
    print_info('    Saving to: ' + final_bed_fpath)
    return final_bed_fpath


def search_trivial_deletions(temp_output_dir, bam_fpath, seq_lengths, need_ref_splitting=False):
    deletions = []
    trivial_deletions_fpath = join(temp_output_dir, qconfig.trivial_deletions_fname)
    print_info('  Looking for trivial deletions (long zero-covered fragments)...')
    need_trivial_deletions = True
    if need_trivial_deletions or need_ref_splitting:
        call_subprocess(['samtools', 'index', bam_fpath])
        bamfile = pysam.AlignmentFile(bam_fpath, "rb")
        cur_deletion = None
        for record in bamfile.fetch():
            mapping = Mapping(record.reference_name, record.reference_start, record.mapping_quality, record.next_reference_name, record.query_alignment_length)
            if mapping:
                if mapping.ref == '*':
                    continue
                # common case: continue current deletion (potential) on the same reference
                if cur_deletion and cur_deletion.ref == mapping.ref:
                    if cur_deletion.next_bad is None:  # previous mapping was in region BEFORE 0-covered fragment
                        # just passed 0-covered fragment
                        if mapping.start - cur_deletion.prev_bad > QuastDeletion.MIN_GAP:
                            cur_deletion.set_next_bad(mapping)
                            if mapping.mapq >= Mapping.MIN_MAP_QUALITY:
                                cur_deletion.set_next_good(mapping)
                                if cur_deletion.is_valid():
                                    deletions.append(cur_deletion)
                                cur_deletion = QuastDeletion(mapping.ref).set_prev_good(mapping)
                        # continue region BEFORE 0-covered fragment
                        elif mapping.mapq >= Mapping.MIN_MAP_QUALITY:
                            cur_deletion.set_prev_good(mapping)
                        else:
                            cur_deletion.set_prev_bad(mapping)
                    else:  # previous mapping was in region AFTER 0-covered fragment
                        # just passed another 0-cov fragment between end of cur_deletion BAD region and this mapping
                        if mapping.start - cur_deletion.next_bad_end > QuastDeletion.MIN_GAP:
                            if cur_deletion.is_valid():  # add previous fragment's deletion if needed
                                deletions.append(cur_deletion)
                            cur_deletion = QuastDeletion(mapping.ref).set_prev_bad(position=cur_deletion.next_bad_end)
                        # continue region AFTER 0-covered fragment (old one or new/another one -- see "if" above)
                        elif mapping.mapq >= Mapping.MIN_MAP_QUALITY:
                            cur_deletion.set_next_good(mapping)
                            if cur_deletion.is_valid():
                                deletions.append(cur_deletion)
                            cur_deletion = QuastDeletion(mapping.ref).set_prev_good(mapping)
                        else:
                            cur_deletion.set_next_bad_end(mapping)
                # special case: just started or just switched to the next reference
                else:
                    if cur_deletion and cur_deletion.ref in seq_lengths:  # switched to the next ref
                        cur_deletion.set_next_good(position=seq_lengths[cur_deletion.ref])
                        if cur_deletion.is_valid():
                            deletions.append(cur_deletion)
                    cur_deletion = QuastDeletion(mapping.ref).set_prev_good(mapping)

                '''if need_ref_splitting:
                    cur_ref = ref_labels[mapping.ref]
                    if mapping.ref_next.strip() == '=' or cur_ref == ref_labels[mapping.ref_next]:
                        if ref_files[cur_ref] is not None:
                            ref_files[cur_ref].write(line)'''
        if cur_deletion and cur_deletion.ref in seq_lengths:  # switched to the next ref
            cur_deletion.set_next_good(position=seq_lengths[cur_deletion.ref])
            if cur_deletion.is_valid():
                deletions.append(cur_deletion)
    # if need_ref_splitting:
    #    for ref_handler in ref_files.values():
    #        if ref_handler is not None:
    #            ref_handler.close()
    if need_trivial_deletions:
        print_info('  Trivial deletions: %d found' % len(deletions))
        print_info('    Saving to: ' + trivial_deletions_fpath)
        with open(trivial_deletions_fpath, 'w') as f:
            for deletion in deletions:
                f.write(str(deletion) + '\n')
    return trivial_deletions_fpath


def align_reference(ref_fpath, output_dir, using_reads='all', max_threads=4, calculate_coverage=False):
    required_files = []
    ref_name = qutils.name_from_fpath(ref_fpath)
    cov_fpath = qconfig.cov_fpath or join(output_dir, ref_name + '.cov')
    uncovered_fpath = add_suffix(cov_fpath, 'uncovered')
    if using_reads != 'all':
        cov_fpath = add_suffix(cov_fpath, using_reads)
        uncovered_fpath = add_suffix(uncovered_fpath, using_reads)
    insert_size_fpath = join(output_dir, ref_name + '.is.txt')
    if not is_non_empty_file(uncovered_fpath):
        required_files.append(uncovered_fpath)
    if not is_non_empty_file(insert_size_fpath) and (using_reads == 'all' or using_reads == 'pe'):
        required_files.append(insert_size_fpath)

    temp_output_dir = join(output_dir, 'tmp')
    if not isdir(temp_output_dir):
        os.makedirs(temp_output_dir)

    correct_chr_names, sam_fpath, bam_fpath = align_single_file(ref_fpath, output_dir, temp_output_dir,
                                                                max_threads, sam_fpath=qconfig.reference_sam,
                                                                bam_fpath=qconfig.reference_bam, required_files=required_files,
                                                                alignment_only=True, reads_type=using_reads)
    if not qconfig.optimal_assembly_insert_size or qconfig.optimal_assembly_insert_size == 'auto':
        if using_reads == 'pe' and sam_fpath:
            insert_size, _, _ = calculate_insert_size(sam_fpath, output_dir, ref_name)
            if not insert_size:
                print_info('  Failed calculating insert size.')
            else:
                qconfig.optimal_assembly_insert_size = insert_size
        elif using_reads == 'all' and is_non_empty_file(insert_size_fpath):
            try:
                insert_size = int(open(insert_size_fpath).readline())
                if insert_size:
                    qconfig.optimal_assembly_insert_size = insert_size
            except:
                pass

    if not required_files:
        return sam_fpath, bam_fpath, uncovered_fpath
    if not sam_fpath:
        print_info('  Failed detecting uncovered regions.')
        return None, None, None

    if calculate_coverage:
        bam_mapped_fpath = get_safe_fpath(temp_output_dir, add_suffix(bam_fpath, 'mapped'))
        bam_sorted_fpath = get_safe_fpath(temp_output_dir, add_suffix(bam_mapped_fpath, 'sorted'))

        if is_non_empty_file(bam_sorted_fpath):
            print_info('  Using existing sorted BAM-file: ' + bam_sorted_fpath)
        else:
            convert_sam(bam_fpath, bam_mapped_fpath, max_threads, filter_rule='-F 4') #not unmapped
            sort_bam(bam_mapped_fpath, bam_sorted_fpath)
        if not is_non_empty_file(uncovered_fpath) and calculate_coverage:
            get_coverage(temp_output_dir, ref_fpath, ref_name, bam_fpath, bam_sorted_fpath,
                         correct_chr_names, max_threads, cov_fpath, uncovered_fpath=uncovered_fpath, create_cov_files=False)
    return sam_fpath, bam_fpath, uncovered_fpath


def parse_reads_stats(stats_fpath):
    reads_stats = defaultdict(int)
    reads_stats['coverage_thresholds'] = []
    with open(stats_fpath) as f:
        for line in f:
            value = line.split()[0]
            if 'total' in line:
                reads_stats['total'] = int(value)
            elif 'secondary' in line:
                reads_stats['total'] -= int(value)
            elif 'supplementary' in line:
                reads_stats['total'] -= int(value)
            elif 'duplicates' in line:
                reads_stats['total'] -= int(value)
            elif 'read1' in line:
                reads_stats['right'] = value
            elif 'read2' in line:
                reads_stats['left'] = value
            elif 'mapped' in line and '%' in line:
                reads_stats['mapped'] = value
                reads_stats['mapped_pcnt'] = get_pcnt_reads(value, reads_stats['total'])
            elif 'properly paired' in line:
                reads_stats['paired'] = value
                reads_stats['paired_pcnt'] = get_pcnt_reads(value, reads_stats['total'])
            elif 'singletons' in line:
                reads_stats['singletons'] = value
                reads_stats['singletons_pcnt'] = get_pcnt_reads(value, reads_stats['total'])
            elif 'different chr' in line and 'mapQ' not in line:
                reads_stats['misjoint'] = value
                reads_stats['misjoint_pcnt'] = get_pcnt_reads(value, reads_stats['total'])
            elif 'depth' in line:
                reads_stats['depth'] = value
            elif 'coverage' in line:
                reads_stats['coverage_thresholds'].append(float(value))
    return reads_stats


def get_pcnt_reads(reads, total_reads):
    return float('%.2f' % (int(reads) * 100.0 / total_reads)) if total_reads != 0 else None


def parse_read_stats(labels, reports, output_dir, ref_fpath):
    ref_reads_stats = None
    if ref_fpath:
        ref_name = qutils.name_from_fpath(ref_fpath)
        stats_fpath = join(output_dir, ref_name + '.stat')
        if isfile(stats_fpath):
            ref_reads_stats = parse_reads_stats(stats_fpath)
            if int(ref_reads_stats['mapped']) == 0:
                print_info('  BWA: nothing aligned for reference.')

    # process all contigs files
    for label in labels:
        report = reports[label]
        stats_fpath = join(output_dir, label + '.stat')
        if ref_reads_stats:
            report.add_field(reporting.Fields.REF_MAPPED_READS, ref_reads_stats['mapped'])
            report.add_field(reporting.Fields.REF_MAPPED_READS_PCNT, ref_reads_stats['mapped_pcnt'])
            report.add_field(reporting.Fields.REF_PROPERLY_PAIRED_READS, ref_reads_stats['paired'])
            report.add_field(reporting.Fields.REF_PROPERLY_PAIRED_READS_PCNT, ref_reads_stats['paired_pcnt'])
            report.add_field(reporting.Fields.REF_SINGLETONS, ref_reads_stats['singletons'])
            report.add_field(reporting.Fields.REF_SINGLETONS_PCNT, ref_reads_stats['singletons_pcnt'])
            report.add_field(reporting.Fields.REF_MISJOINT_READS, ref_reads_stats['misjoint'])
            report.add_field(reporting.Fields.REF_MISJOINT_READS_PCNT, ref_reads_stats['misjoint_pcnt'])
            report.add_field(reporting.Fields.REF_DEPTH, ref_reads_stats['depth'])
            if ref_reads_stats['coverage_thresholds'] and len(ref_reads_stats['coverage_thresholds']) == len(qconfig.coverage_thresholds):
                report.add_field(reporting.Fields.REF_COVERAGE__FOR_THRESHOLDS,
                                [ref_reads_stats['coverage_thresholds'][i] for i, threshold in enumerate(qconfig.coverage_thresholds)])
                report.add_field(reporting.Fields.REF_COVERAGE_1X_THRESHOLD, ref_reads_stats['coverage_thresholds'][0])
        if not isfile(stats_fpath):
            continue
        reads_stats = parse_reads_stats(stats_fpath)
        report.add_field(reporting.Fields.TOTAL_READS, reads_stats['total'])
        report.add_field(reporting.Fields.LEFT_READS, reads_stats['left'])
        report.add_field(reporting.Fields.RIGHT_READS, reads_stats['right'])
        report.add_field(reporting.Fields.MAPPED_READS, reads_stats['mapped'])
        report.add_field(reporting.Fields.MAPPED_READS_PCNT, reads_stats['mapped_pcnt'])
        report.add_field(reporting.Fields.PROPERLY_PAIRED_READS, reads_stats['paired'])
        report.add_field(reporting.Fields.PROPERLY_PAIRED_READS_PCNT, reads_stats['paired_pcnt'])
        if int(reads_stats['mapped']) == 0:
            print_info('  BWA: nothing aligned for ' + '\'' + label + '\'.')
        report.add_field(reporting.Fields.SINGLETONS, reads_stats['singletons'])
        report.add_field(reporting.Fields.SINGLETONS_PCNT, reads_stats['singletons_pcnt'])
        report.add_field(reporting.Fields.MISJOINT_READS, reads_stats['misjoint'])
        report.add_field(reporting.Fields.MISJOINT_READS_PCNT, reads_stats['misjoint_pcnt'])
        report.add_field(reporting.Fields.DEPTH, reads_stats['depth'])
        if reads_stats['coverage_thresholds'] and len(reads_stats['coverage_thresholds']) == len(qconfig.coverage_thresholds):
            report.add_field(reporting.Fields.COVERAGE__FOR_THRESHOLDS,
                            [reads_stats['coverage_thresholds'][i] for i, threshold in enumerate(qconfig.coverage_thresholds)])
            report.add_field(reporting.Fields.COVERAGE_1X_THRESHOLD, reads_stats['coverage_thresholds'][0])
    save_reads(reports, output_dir)


def get_physical_coverage(output_dirpath, ref_name, bam_fpath, cov_fpath, chr_len_fpath, threads):
    raw_cov_fpath = add_suffix(cov_fpath, 'raw')
    if not is_non_empty_file(raw_cov_fpath):
        print_info('  Calculating physical coverage...')
        ## keep properly mapped, unique, non-duplicate paired-end reads only
        bam_filtered_fpath = join(output_dirpath, ref_name + '.physical.bam')
        convert_sam(bam_fpath, bam_filtered_fpath, threads,
                      filter_rule='-f 2 -F 256') #not supplementary and proper_pair
        ## sort by read names
        bam_filtered_sorted_fpath = join(output_dirpath, ref_name + '.physical.sorted.bam')
        sort_bam(bam_filtered_fpath, bam_filtered_sorted_fpath, sort_rule='-n')
        bed_fpath = bam_to_bed(output_dirpath, ref_name + '.physical', bam_filtered_sorted_fpath, bedpe=True)
        calculate_genome_cov(bed_fpath, raw_cov_fpath, chr_len_fpath)
    return raw_cov_fpath


def get_coverage(output_dirpath, ref_fpath, ref_name, bam_fpath, bam_sorted_fpath, correct_chr_names, max_threads,
                 cov_fpath, physical_cov_fpath=None, uncovered_fpath=None, create_cov_files=True):
    raw_cov_fpath = cov_fpath + '_raw'
    chr_len_fpath = get_chr_len_fpath(ref_fpath, correct_chr_names)
    if not is_non_empty_file(cov_fpath):
        print_info('  Calculating reads coverage...')
        if not is_non_empty_file(raw_cov_fpath):
            if not is_non_empty_file(bam_sorted_fpath):
                sort_bam(bam_fpath, bam_sorted_fpath)
            calculate_genome_cov(bam_sorted_fpath, raw_cov_fpath, chr_len_fpath)
            assert_file_exists(raw_cov_fpath, 'coverage file')
        if uncovered_fpath:
            print_uncovered_regions(raw_cov_fpath, uncovered_fpath, correct_chr_names)
        if create_cov_files:
            proceed_cov_file(raw_cov_fpath, cov_fpath, correct_chr_names)
    if not is_non_empty_file(physical_cov_fpath) and create_cov_files:
        raw_cov_fpath = get_physical_coverage(output_dirpath, ref_name, bam_fpath,
                                              physical_cov_fpath, chr_len_fpath, max_threads)
        proceed_cov_file(raw_cov_fpath, physical_cov_fpath, correct_chr_names)
    return cov_fpath, physical_cov_fpath


def proceed_cov_file(raw_cov_fpath, cov_fpath, correct_chr_names):
    chr_depth = defaultdict(list)
    used_chromosomes = dict()
    chr_index = 0
    with open(raw_cov_fpath, 'r') as in_coverage:
        with open(cov_fpath, 'w') as out_coverage:
            for line in in_coverage:
                fs = list(line.split())
                name = fs[0]
                depth = int(float(fs[-1]))
                if name not in used_chromosomes:
                    chr_index += 1
                    used_chromosomes[name] = str(chr_index)
                    correct_name = correct_chr_names[name] if correct_chr_names else name
                    out_coverage.write('#' + correct_name + ' ' + used_chromosomes[name] + '\n')
                if len(fs) > 3:
                    start, end = int(fs[1]), int(fs[2])
                    chr_depth[name].extend([depth] * (end - start))
                else:
                    chr_depth[name].append(depth)
                if len(chr_depth[name]) >= COVERAGE_FACTOR:
                    max_index = len(chr_depth[name]) - (len(chr_depth[name]) % COVERAGE_FACTOR)
                    for index in range(0, max_index, COVERAGE_FACTOR):
                        cur_depth = sum(chr_depth[name][index: index + COVERAGE_FACTOR]) // COVERAGE_FACTOR
                        out_coverage.write(' '.join([used_chromosomes[name], str(cur_depth) + '\n']))
                    chr_depth[name] = chr_depth[name][index + COVERAGE_FACTOR:]
            if not qconfig.debug:
                os.remove(raw_cov_fpath)


def print_uncovered_regions(raw_cov_fpath, uncovered_fpath, correct_chr_names):
    uncovered_regions = defaultdict(list)
    with open(raw_cov_fpath) as in_coverage:
        for line in in_coverage:
            fs = list(line.split())
            name = fs[0]
            depth = int(float(fs[-1]))
            correct_name = correct_chr_names[name] if correct_chr_names else name
            if len(fs) > 3 and depth == 0:
                uncovered_regions[correct_name].append((fs[1], fs[2]))
    with open(uncovered_fpath, 'w') as out_f:
        for chrom, regions in uncovered_regions.items():
            for start, end in regions:
                out_f.write('\t'.join([chrom, start, end]) + '\n')

