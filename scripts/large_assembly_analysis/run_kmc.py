import os
import sys
from os.path import join, isdir

from src import qconfig
from src.common import parse_ref_stats
from src.kmers_analyzer import *
from src.kmers_analyzer import _get_dist_inconstistency
from src.logger import *


def main():
    output_dirpath, tmp_dirpath, ref_kmc_out_fpath, reference_csv, is_cyclic, \
    downsampled_kmers_fpath, contigs_fpath, label, threads = sys.argv[1:]
    threads = int(threads)

    is_cyclic = True if is_cyclic == 'True' else False
    genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(reference_csv, skip_ns=True)

    print_info('  Analyzing assemblies completeness...')
    kmer_len = qconfig.unique_kmer_len

    print_info('    ' + label)

    kmc_out_fpath = count_kmers(tmp_dirpath, contigs_fpath, kmer_len, threads)
    intersect_out_fpath = intersect_kmers(tmp_dirpath, [ref_kmc_out_fpath, kmc_out_fpath], threads)
    matched_kmers = get_kmers_cnt(tmp_dirpath, intersect_out_fpath, threads)
    unique_kmers = get_kmers_cnt(tmp_dirpath, ref_kmc_out_fpath, threads, run_histo=False)
    completeness = matched_kmers * 100.0 / unique_kmers

    print_info('  Analyzing assemblies correctness...')

    corr_len = None
    mis_len = None
    undef_len = None
    translocations, relocations = None, None
    total_len = 0
    contig_lens = dict()
    for record in read_fasta(contigs_fpath):
        total_len += len(record.seq)
        contig_lens[record.id] = len(record.seq)

    if not downsampled_kmers_fpath or not exists(downsampled_kmers_fpath):
        print_warning('Scaffolding accuracy will not be assessed.')
    else:
        corr_len = 0
        mis_len = 0
        kmers_by_contig, kmers_pos_by_contig = align_kmers(label, tmp_dirpath, contigs_fpath, downsampled_kmers_fpath,
                                                           threads)
        cyclic_ref_lens = genome_size if is_cyclic else None
        translocations = 0
        relocations = 0
        with open(join(tmp_dirpath, label + '.misjoins.txt'), 'w') as out:
            for contig in kmers_by_contig.keys():
                contig_markers = []
                prev_pos, prev_ref_pos, prev_chrom, marker = None, None, None, None
                for pos, kmer_name in sorted(zip(kmers_pos_by_contig[contig], kmers_by_contig[contig]), key=lambda x: x[0]):
                    ref_chrom, ref_pos = kmer_name.rsplit('_', 1)
                    ref_pos = int(ref_pos)
                    if prev_pos and prev_chrom:
                        if prev_chrom == ref_chrom and abs(
                                abs(pos - prev_pos) / abs(ref_pos - prev_ref_pos) - 1) <= 0.05:
                            marker = (pos, ref_pos, ref_chrom)
                        elif marker:
                            contig_markers.append(marker)
                            pos, ref_pos, ref_chrom, marker = None, None, None, None
                    prev_pos, prev_ref_pos, prev_chrom = pos, ref_pos, ref_chrom
                if marker:
                    contig_markers.append(marker)
                prev_pos, prev_ref_pos, prev_chrom = None, None, None
                is_misassembled = False
                for marker in contig_markers:
                    pos, ref_pos, ref_chrom = marker
                    if prev_pos and prev_chrom:
                        if ref_chrom != prev_chrom:
                            translocations += 1
                            out.write('Translocation in %s: %s %d | %s %d\n' %
                                      (contig, prev_chrom, prev_pos, ref_chrom, pos))
                            is_misassembled = True
                        elif _get_dist_inconstistency(pos, prev_pos, ref_pos, prev_ref_pos,
                                                      cyclic_ref_lens) > EXT_RELOCATION_SIZE:
                            relocations += 1
                            out.write('Relocation in %s: %d (%d) | %d (%d)\n' %
                                      (contig, prev_pos, prev_ref_pos, pos, ref_pos))
                            is_misassembled = True
                    prev_pos, prev_ref_pos, prev_chrom = pos, ref_pos, ref_chrom
                if is_misassembled:
                    mis_len += contig_lens[contig]
                elif len(contig_markers) > 0:
                    corr_len += contig_lens[contig]
        undef_len = total_len - corr_len - mis_len

    create_kmc_stats_file(output_dirpath, contigs_fpath, label, completeness,
                              corr_len, mis_len, undef_len, total_len, translocations, relocations)


if __name__ == '__main__':
    main()