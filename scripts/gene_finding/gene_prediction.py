############################################################################
# Copyright (c) 2015-2020 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import tempfile
import itertools
import csv
import shutil
from os.path import exists

from src import qutils, qconfig
from src.common import open_gzipsafe, read_fasta, write_fasta, save_csv_from_dict, rev_comp
from src.genes_parser import Gene
from src.logger import *
from src.qutils import get_path_to_program


def merge_gffs(gffs, out_path):
    '''Merges all GFF files into a single one, dropping GFF header.'''
    out_file = open(out_path, 'w')
    out_file.write('##gff-version 3\n')
    for gff_path in gffs:
        with open(gff_path) as gff_file:
            out_file.writelines(itertools.islice(gff_file, 2, None))
    out_file.close()
    return out_path


def parse_gff(gff_path):
    gff_file = open_gzipsafe(gff_path)
    r = csv.reader(list(filter(lambda l: not l.startswith("#"), gff_file)),
        delimiter='\t')
    for index, _source, type, start, end, score, strand, phase, extra in r:
        if type != 'mRNA':
            continue  # We're only interested in genes here.

        attrs = dict(kv.split("=") for kv in extra.split(";"))
        yield index, attrs.get('Name'), int(start), int(end), strand
    gff_file.close()


def analyze_gff(out_gff_path, contigs):
    unique, total = set(), 0
    genes = []
    for contig, gene_id, start, end, strand in parse_gff(out_gff_path):
        total += 1
        if strand == '+':
            gene_seq = contigs[contig][start - 1:end]
        else:
            gene_seq = rev_comp(contigs[contig][start - 1:end])
        if gene_seq not in unique:
            unique.add(gene_seq)
        gene = Gene(contig=contig, start=start, end=end, strand=strand, seq=gene_seq)
        gene.is_full = gene.start > 1 and gene.end < len(contigs[contig])
        genes.append(gene)

    gene_lengths = [int(x) for x in qconfig.gene_lengths.split(",")]
    full_cnt = [sum([gene.end - gene.start >= threshold for gene in genes if gene.is_full]) for threshold in gene_lengths]
    partial_cnt = [sum([gene.end - gene.start >= threshold for gene in genes if not gene.is_full]) for threshold in gene_lengths]

    return genes, len(unique), total, full_cnt, partial_cnt


def run_glimmerHMM(assembly_label, contigs_fpath, gff_fpath, tmp_dir):
    tool_exec_fpath = get_path_to_program('glimmerhmm')
    if not tool_exec_fpath:
        return None, None

    print_info('  ' + assembly_label)

    # Note: why arabidopsis? for no particular reason, really.
    trained_dir = tool_exec_fpath.replace('bin/glimmerhmm', os.path.join('share', 'glimmerhmm', 'trained_dir', 'arabidopsis'))
    if not exists(trained_dir):
        trained_dir = tool_exec_fpath.replace('bin/glimmerhmm',
                                              os.path.join('GlimmerHMM', 'trained_dir', 'arabidopsis'))
        if not exists(trained_dir):
            print_error('Trained dir ' + trained_dir + 'was not found, please check the installation')
    contigs = {}
    gffs = []
    base_dir = tempfile.mkdtemp(dir=tmp_dir)
    for seq_num, record in enumerate(read_fasta(contigs_fpath)):
        seq_num = str(seq_num)
        record.id = record.id[:qutils.MAX_CONTIG_NAME_GLIMMER]
        contig_path = os.path.join(base_dir, seq_num + '.fasta')
        gff_path = os.path.join(base_dir, seq_num + '.gff')

        write_fasta(contig_path, [record])
        return_code = qutils.call_subprocess(
            [tool_exec_fpath, contig_path, '-d', trained_dir, '-g', '-o', gff_path],
            stdout=sys.stderr,
            indent='    ')
        if return_code == 0:
            gffs.append(gff_path)
            contigs[record.id] = str(record.seq)

    if not gffs:
        return None, None

    merge_gffs(gffs, gff_fpath)
    return gff_fpath, contigs


def run_prodigal(contigs_fpath, gff_fpath):
    tool_exec_fpath = get_path_to_program('prodigal')
    if not tool_exec_fpath:
        print_error('Cannot find Prodigal tool. Please install it or turn off gene prediction.')
        return None, None

    contigs = {}
    for seq_num, record in enumerate(read_fasta(contigs_fpath)):
        contigs[record.id] = str(record.seq)
    return_code = qutils.call_subprocess(
        [tool_exec_fpath, '-i', contigs_fpath, '-f', 'gff', '-o', gff_fpath],
            stdout=sys.stderr,
            indent='    ')
    if return_code == 0:
        return gff_fpath, contigs
    else:
        # try to run prodigal in anonymous mode (not enough data)
        return_code = qutils.call_subprocess(
            [tool_exec_fpath, '-i', contigs_fpath, '-f', 'gff', '-o', gff_fpath, '-p', 'anon'],
                stdout=sys.stderr,
                indent='    ')
        if return_code == 0:
            return gff_fpath, contigs
        else:
            return None, None


def main():
    output_dirpath, contigs_fpath, assembly_label, tmp_dirpath, tool = sys.argv[1:]
    tmp_dirpath = os.path.join(output_dirpath, 'tmp')

    gff_fpath = os.path.join(output_dirpath, assembly_label + '.gff')
    if tool == 'prodigal':
        gff_fpath, contigs = run_prodigal(contigs_fpath, gff_fpath)
    else:
        gff_fpath, contigs = run_glimmerHMM(assembly_label, contigs_fpath, gff_fpath, tmp_dirpath)

    if gff_fpath:
        genes, unique, total, full_genes, partial_genes = analyze_gff(gff_fpath, contigs)
        print_info('    Genes = ' + str(unique) + ' unique, ' + str(total) + ' total')
        print_info('    Predicted genes (GFF): ' + gff_fpath)

        stats = {'genes': [g.asdict() for g in genes], 'unique': [unique], 'full': full_genes, 'partial': partial_genes}
        result_fpath = os.path.join(tmp_dirpath, assembly_label + '.csv')
        save_csv_from_dict(stats, result_fpath)


if __name__ == '__main__':
    main()
