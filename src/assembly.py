import os
from ast import literal_eval
from os.path import join, isdir

import pandas as pd

from src.common import save_csv_from_dict


class Assembly:
    def __init__(self, output_dir, contigs_fpath, label):
        self.contigs_fpath = contigs_fpath
        self.label = label
        if output_dir:
            self.aux_dir = join(output_dir, 'aux')
            self.status_file = join(self.aux_dir, label + ".status.csv")
            self.result_file = join(self.aux_dir, label + ".result.csv")
            self.aligned_stats_file = join(self.aux_dir, label + ".aligned_stats.csv")

    def add_results(self, status, results, aligned_lengths, misassemblies_in_contigs, aligned_lengths_by_contigs):
        self.status, self.results, self.aligned_lengths, self.misassemblies_in_contigs, self.aligned_lengths_by_contigs = \
            status, results, aligned_lengths, misassemblies_in_contigs, aligned_lengths_by_contigs

    def save_results(self):
        if not isdir(self.aux_dir):
            os.makedirs(self.aux_dir)
        save_csv_from_dict(self.results, self.result_file)
        pd.DataFrame.from_dict({'status': [self.status], 'aligned_lengths': self.aligned_lengths,
                                'misassemblies_in_contigs':self.misassemblies_in_contigs,
                                'aligned_lengths_by_contigs': self.aligned_lengths_by_contigs},
                               orient='index').to_csv(self.aligned_stats_file)

    def literal_converter(self, val):
        try:
            return literal_eval(val)
        except ValueError:
            return val

    def parse_results(self):
        self.results = pd.read_csv(self.result_file, index_col=0, header=0, names=['val'], converters={'val': self.literal_converter}).squeeze("columns")
        df = pd.read_csv(self.aligned_stats_file, index_col=0)
        self.status = df.loc['status'][0]
        self.aligned_lengths = df.loc['aligned_lengths'].dropna().to_list()
        self.misassemblies_in_contigs = df.loc['misassemblies_in_contigs'].dropna().to_list()
        self.aligned_lengths_by_contigs = df.loc['aligned_lengths_by_contigs'].dropna().to_list()

