############################################################################
# Copyright (c) 2015-2020 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
#
# This is auxiliary file for contigs_analyzer.py
#
############################################################################

from __future__ import with_statement
import gzip
import os
import re
from itertools import repeat
from os.path import isdir, join, basename

from collections import OrderedDict

from src import qconfig
from src.qutils import val_to_str

contig_aligner_dirpath = join(qconfig.LIBS_LOCATION, 'minimap2')
ref_labels_by_chromosomes = OrderedDict()
intergenomic_misassemblies_by_asm = {}
contigs_aligned_lengths = {}


def is_same_reference(chr1, chr2):
    return ref_labels_by_chromosomes[chr1] == ref_labels_by_chromosomes[chr2]


def get_ref_by_chromosome(chrom):
    return ref_labels_by_chromosomes[chrom] if chrom in ref_labels_by_chromosomes else 'ref'


def parse_cs_tag(cigar):
    cs_pattern = re.compile(r':\d+|\*[acgtn]+|\-[acgtn]+|\+[acgtn]+')
    return cs_pattern.findall(cigar)


def print_file(all_rows, fpath, append_to_existing_file=False):
    colwidths = repeat(0)
    for row in all_rows:
        colwidths = [max(len(v), w) for v, w in zip([row['metricName']] + [val_to_str(v) for v in row['values']], colwidths)]
    txt_file = open(fpath, 'a' if append_to_existing_file else 'w')
    for row in all_rows:
        txt_file.write('  '.join('%-*s' % (colwidth, cell) for colwidth, cell
                                     in zip(colwidths, [row['metricName']] + [val_to_str(v) for v in row['values']])) + '\n')


def create_minimap_output_dir(output_dir):
    minimap_output_dir = join(output_dir, qconfig.aligner_output_dirname)
    '''if qconfig.is_combined_ref:
        from quast_libs import search_references_meta
        if search_references_meta.is_quast_first_run:
            minimap_output_dir = join(minimap_output_dir, 'raw')
            if not os.path.isdir(minimap_output_dir):
                os.mkdir(minimap_output_dir)'''
    return minimap_output_dir


def close_handlers(ca_output):
    for handler in vars(ca_output).values():  # assume that ca_output does not have methods (fields only)
        if handler:
            handler.close()

