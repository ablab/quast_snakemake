import ast
from collections import OrderedDict
from os.path import splitext

import pandas as pd
import gzip
from functools import partial

from Bio import SeqIO


def open_gzipsafe(filename):
    f = partial(gzip.open, mode='rt') if filename.endswith('.gz') else open
    return f(filename)


def read_fasta(input_file):
    with open_gzipsafe(input_file) as f:
        for record in SeqIO.parse(f, 'fasta'):
            yield record


def write_fasta(output_file, sequences):
    SeqIO.write(sequences, output_file, "fasta")


def get_chr_lengths_from_fastafile(input_file):
    chr_lengths = OrderedDict()
    for record in read_fasta(input_file):
        chr_lengths[record.id] = len(record.seq)
    return chr_lengths


def parse_ref_stats(reference_csv, skip_ns=True):
    df = pd.read_csv(reference_csv, index_col=0)
    genome_size = sum(df['length'])
    if skip_ns:
        genome_size -= sum(df['Ns'].str.len())
    reference_chromosomes = df['length'].to_dict()
    ns_by_chromosomes = df['Ns'].apply(lambda x: ast.literal_eval(x)).to_dict()
    return genome_size, reference_chromosomes, ns_by_chromosomes


def save_ref_stats(fasta_fpath, csv_fpath):
    ref_stats = dict()
    for record in read_fasta(fasta_fpath):
        chr_name = record.id.split()[0]
        ns_coords = [x + 1 for x, s in enumerate(record.seq) if s == 'N']
        ref_stats[chr_name] = [len(record.seq), ns_coords]
    pd.DataFrame.from_dict(ref_stats, orient='index', columns=['length','Ns']).to_csv(csv_fpath)


def get_labels_from_paths(contigs_fpaths):
    return [splitext(f)[0].split('/')[-1] for f in contigs_fpaths]


def save_csv(df, filename):
    df.to_csv(filename)


def save_csv_from_dict(d, filename):
    pd.DataFrame.from_dict(dict([(k,pd.Series(v)) for k,v in d.items()]),orient='index').to_csv(filename)
