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


def glimmerHMM(tool_exec_fpath, fasta_fpath, out_fpath, err_path, tmp_dir):
    def run(contig_path, tmp_path):
        with open(err_path, 'a') as err_file:
            return_code = qutils.call_subprocess(
                [tool_exec_fpath, contig_path, '-d', trained_dir, '-g', '-o', tmp_path],
                stdout=err_file,
                stderr=err_file,
                indent='    ')
            return return_code

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
    for seq_num, record in enumerate(read_fasta(fasta_fpath)):
        seq_num = str(seq_num)
        record.id = record.id[:qutils.MAX_CONTIG_NAME_GLIMMER]
        contig_path = os.path.join(base_dir, seq_num + '.fasta')
        gff_path = os.path.join(base_dir, seq_num + '.gff')

        write_fasta(contig_path, [record])
        if run(contig_path, gff_path) == 0:
            gffs.append(gff_path)
            contigs[record.id] = str(record.seq)

    if not gffs:
        return None, None, None, None, None, None

    out_gff_fpath = out_fpath + '.gff'
    out_gff_path = merge_gffs(gffs, out_gff_fpath)
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

    if not qconfig.debug:
        shutil.rmtree(base_dir)

    return out_gff_path, genes, len(unique), total, full_cnt, partial_cnt


def main():
    output_dirpath, contigs_fpath, assembly_label, tmp_dirpath = sys.argv[1:]
    tmp_dirpath = os.path.join(output_dirpath, 'tmp')
    tool_exec_fpath = get_path_to_program('glimmerhmm')
    if not tool_exec_fpath:
        return

    print_info('  ' + assembly_label)

    out_fpath = os.path.join(output_dirpath, assembly_label + '_glimmer')
    err_fpath = os.path.join(output_dirpath, assembly_label + '_glimmer.stderr')

    out_gff_path, genes, unique, total, full_genes, partial_genes = glimmerHMM(tool_exec_fpath,
        contigs_fpath, out_fpath, err_fpath, tmp_dirpath)

    if out_gff_path:
        print_info('    Genes = ' + str(unique) + ' unique, ' + str(total) + ' total')
        print_info('    Predicted genes (GFF): ' + out_gff_path)

    stats = {'genes': [g.asdict() for g in genes], 'unique': [unique], 'full': full_genes, 'partial': partial_genes}
    result_fpath = os.path.join(tmp_dirpath, assembly_label + '_glimmer.csv')
    save_csv_from_dict(stats, result_fpath)


if __name__ == '__main__':
    main()
