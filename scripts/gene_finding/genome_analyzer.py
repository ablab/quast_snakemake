############################################################################
# Copyright (c) 2015-2020 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
from collections import defaultdict
from os.path import join, exists

from scripts.gene_finding.prepare_genome_analyzer import parse_results
from src import qconfig
from src.common import read_fasta, parse_ref_stats, save_csv_from_dict
from src.logger import *


class AlignedBlock():
    def __init__(self, seqname=None, start=None, end=None, contig=None, start_in_contig=None, end_in_contig=None):
        self.seqname = seqname
        self.start = start
        self.end = end
        self.contig = contig
        self.start_in_contig = start_in_contig
        self.end_in_contig = end_in_contig

    def format_gene_info(self, region):
        start, end = self.start_in_contig, self.end_in_contig
        if self.start < region.start:
            region_shift = region.start - self.start
            if start < end:
                start += region_shift
            else:
                start -= region_shift
        if region.end < self.end:
            region_size = region.end - max(region.start, self.start)
            if start < end:
                end = start + region_size
            else:
                end = start - region_size
        return self.contig + ':' + str(start) + '-' + str(end)


def main():
    output_dirpath, reference_csv, contigs_fpath, assembly_label, coords_dirpath = sys.argv[1:6]
    container_fnames = sys.argv[6:]

    genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(reference_csv)

    result_fpath = join(output_dirpath, assembly_label + '_info.txt')
    contig_result_fpath = join(output_dirpath, assembly_label + '_contig_info.txt')
    ref_result_fpath = join(output_dirpath, assembly_label + '_ref_info.txt')
    print_info('  ' + assembly_label)

    containers = []
    for container_fname in container_fnames:
        containers.append(parse_results(container_fname))

    results = dict()
    ref_lengths = defaultdict(int)

    coords_base_fpath = os.path.join(coords_dirpath, assembly_label + '.coords')
    if qconfig.use_all_alignments:
        coords_fpath = coords_base_fpath
    else:
        coords_fpath = coords_base_fpath + '.filtered'

    if not os.path.isfile(coords_fpath):
        print_error('File with alignment coords (' + coords_fpath + ') not found! Try to restart QUAST.',
            indent='  ')
        return None, None

    # EXAMPLE:
    #    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
    #=====================================================================================
    #  338980   339138  |     2298     2134  |      159      165  |    79.76  | gi|48994873|gb|U00096.2|	NODE_0_length_6088
    #  374145   374355  |     2306     2097  |      211      210  |    85.45  | gi|48994873|gb|U00096.2|	NODE_0_length_6088

    genome_mapping = {}
    for chr_name, chr_len in reference_chromosomes.items():
        genome_mapping[chr_name] = [0] * (chr_len + 1)

    contig_tuples = list(read_fasta(contigs_fpath))  # list of FASTA entries (in tuples: name, seq)
    sorted_contig_tuples = sorted(enumerate(contig_tuples), key=lambda x: len(x[1].seq), reverse=True)
    sorted_contigs_names = []
    contigs_order = []
    for idx, record in sorted_contig_tuples:
        sorted_contigs_names.append(record.id)
        contigs_order.append(idx)

    features_in_contigs = [0] * len(sorted_contigs_names)  # for cumulative plots: i-th element is the number of genes in i-th contig
    operons_in_contigs = [0] * len(sorted_contigs_names)
    aligned_blocks_by_contig_name = {} # for gene finding: contig_name --> list of AlignedBlock

    gene_searching_enabled = len(containers)
    if qconfig.memory_efficient and gene_searching_enabled:
        print_warning('Analysis of genes and/or operons files (provided with -g and -O) requires extensive RAM usage, consider running QUAST without them if memory consumption is critical.')
    if gene_searching_enabled:
        for name in sorted_contigs_names:
            aligned_blocks_by_contig_name[name] = []
    with open(coords_fpath) as coordfile:
        for line in coordfile:
            s1 = int(line.split('|')[0].split()[0])
            e1 = int(line.split('|')[0].split()[1])
            s2 = int(line.split('|')[1].split()[0])
            e2 = int(line.split('|')[1].split()[1])
            contig_name = line.split()[12].strip()
            chr_name = line.split()[11].strip()

            if chr_name not in genome_mapping:
                print_error("Something went wrong and chromosome names in your coords file (" + coords_base_fpath + ") " \
                             "differ from the names in the reference. Try to remove the file and restart QUAST.")
                return None

            if gene_searching_enabled:
                aligned_blocks_by_contig_name[contig_name].append(AlignedBlock(seqname=chr_name, start=s1, end=e1,
                                                                               contig=contig_name, start_in_contig=s2, end_in_contig=e2))
            for i in range(s1, e1 + 1):
                genome_mapping[chr_name][i] = 1

    for chr_name in genome_mapping.keys():
        for i in ns_by_chromosomes[chr_name]:
            genome_mapping[chr_name][i] = 0
        ref_lengths[chr_name] = sum(genome_mapping[chr_name])

    if qconfig.space_efficient and coords_fpath.endswith('.filtered'):
        os.remove(coords_fpath)

    # counting genome coverage and gaps number
    gaps_count = 0
    if qconfig.analyze_gaps:
        gaps_fpath = os.path.join(output_dirpath, assembly_label + '_gaps.txt') if not qconfig.space_efficient else '/dev/null'
        with open(gaps_fpath, 'w') as gaps_file:
            for chr_name, chr_len in reference_chromosomes.items():
                gaps_file.write(chr_name + '\n')
                cur_gap_size = 0
                for i in range(1, chr_len + 1):
                    if genome_mapping[chr_name][i] == 1 or i in ns_by_chromosomes[chr_name]:
                        if cur_gap_size >= qconfig.min_gap_size:
                            gaps_count += 1
                            gaps_file.write(str(i - cur_gap_size) + ' ' + str(i - 1) + '\n')
                        cur_gap_size = 0
                    else:
                        cur_gap_size += 1
                if cur_gap_size >= qconfig.min_gap_size:
                    gaps_count += 1
                    gaps_file.write(str(chr_len - cur_gap_size + 1) + ' ' + str(chr_len) + '\n')

    # finding genes and operons
    stats = dict()
    contig_stats = dict()
    contig_stats['features'] = [0] * len(sorted_contigs_names)
    contig_stats['operons'] = [0] * len(sorted_contigs_names)
    ref_genes_num = 0
    ref_operons_num = 0
    for container in containers:
        if not container.region_list:
            continue

        if container.kind == 'operon':
            ref_operons_num = len(container.region_list)
        else:
            ref_genes_num += len(container.region_list)

        total_full = 0
        total_partial = 0
        found_fpath = os.path.join(output_dirpath, assembly_label + '_genomic_features_' + container.kind.lower() + '.txt')
        found_file = open(found_fpath, 'w')
        found_file.write('%s\t\t%s\t%s\t%s\t%s\n' % ('ID or #', 'Start', 'End', 'Type', 'Contig'))
        found_file.write('=' * 50 + '\n')

        # 0 - gene is not found,
        # 1 - gene is found,
        # 2 - part of gene is found
        found_list = [0] * len(container.region_list)
        for i, region in enumerate(container.region_list):
            found_list[i] = 0
            gene_blocks = []
            if region.id is None:
                region.id = '# ' + str(region.number + 1)
            for contig_id, name in enumerate(sorted_contigs_names):
                cur_feature_is_found = False
                for cur_block in aligned_blocks_by_contig_name[name]:
                    if cur_block.seqname != region.seqname:
                        continue
                    if region.end <= cur_block.start or cur_block.end <= region.start:
                        continue
                    elif cur_block.start <= region.start and region.end <= cur_block.end:
                        if found_list[i] == 2:  # already found as partial gene
                            total_partial -= 1
                        found_list[i] = 1
                        total_full += 1
                        contig_info = cur_block.format_gene_info(region)
                        found_file.write('%s\t\t%d\t%d\tcomplete\t%s\n' % (region.id, region.start, region.end, contig_info))
                        if container.kind == 'operon':
                            contig_stats['operons'][contig_id] += 1  # inc number of found genes/operons in id-th contig
                        else:
                            contig_stats['features'][contig_id] += 1

                        cur_feature_is_found = True
                        break
                    elif min(region.end, cur_block.end) - max(region.start, cur_block.start) >= qconfig.min_gene_overlap:
                        if found_list[i] == 0:
                            found_list[i] = 2
                            total_partial += 1
                        gene_blocks.append(cur_block)
                    if cur_feature_is_found:
                        break
                if cur_feature_is_found:
                    break
            # adding info about partially found genes/operons
            if found_list[i] == 2:  # partial gene/operon
                contig_info = ','.join([block.format_gene_info(region) for block in sorted(gene_blocks, key=lambda block: block.start)])
                found_file.write('%s\t\t%d\t%d\tpartial\t%s\n' % (region.id, region.start, region.end, contig_info))
        stats[container.kind] = [total_full, total_partial]

    if ref_genes_num:
        for kind, d in stats.items():
            if kind == 'operon':
                stats['operon'].append(ref_operons_num)
            elif kind == 'gene':
                stats['gene'].append(ref_genes_num)
            else: stats[kind].append(0)

    contig_stats['gaps_count'] = [gaps_count]
    contig_stats['features_unsorted'] = [features_in_contigs[idx] for idx in contigs_order]
    contig_stats['operons_unsorted'] = [operons_in_contigs[idx] for idx in contigs_order]
    save_csv_from_dict(stats, result_fpath)
    save_csv_from_dict(contig_stats, contig_result_fpath)
    save_csv_from_dict(ref_lengths, ref_result_fpath, add_header=False)


if __name__ == '__main__':
    main()
