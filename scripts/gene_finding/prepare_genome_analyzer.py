import argparse
import os
from ast import literal_eval
from os.path import join, isdir

import pandas as pd

from src import genes_parser
from src.common import parse_ref_stats

# reading genes and operons
from src.genes_parser import Gene
from src.logger import *


class FeatureContainer:
    def __init__(self, fpaths, kind=''):
        self.kind = kind
        self.fpaths = fpaths
        self.region_list = []
        self.chr_names_dict = {}


def chromosomes_names_dict(feature, regions, chr_names):
    """
    returns dictionary to translate chromosome name in list of features (genes or operons) to
    chromosome name in reference file.
    """
    region_2_chr_name = {}

    for region in regions:
        if region.seqname in chr_names:
            region_2_chr_name[region.seqname] = region.seqname
        else:
            region_2_chr_name[region.seqname] = None

    if len(chr_names) == 1 and len(region_2_chr_name) == 1 and region_2_chr_name[regions[0].seqname] is None:
        chr_name = chr_names.pop()
        print_notice(
            'Reference name in file with genomic features of type "%s" (%s) does not match the name in the reference file (%s). '
            'QUAST will ignore this issue and count as if they match.' %
            (feature, regions[0].seqname, chr_name),
            indent='  ')
        for region in regions:
            region.seqname = chr_name
            region_2_chr_name[region.seqname] = chr_name
    elif all(chr_name is None for chr_name in region_2_chr_name.values()):
        print_warning(
            'Reference names in file with genomic features of type "%s" do not match any chromosome. Check your genomic feature file(s).' % (
                feature),
            indent='  ')
    elif None in region_2_chr_name.values():
        print_warning(
            'Some of the reference names in file with genomic features of type "%s" does not match any chromosome. '
            'Check your genomic feature file(s).' % (feature), indent='  ')

    return region_2_chr_name


def parse_results(containers_file):
    df = pd.read_csv(containers_file, index_col=0)
    container = FeatureContainer(fpaths=None, kind=df.loc['kind'][0])
    genes = df.loc['region_list'].apply(literal_eval).dropna().to_list()
    container.region_list = [Gene(**g) for g in genes]

    df = pd.read_csv(containers_file.replace('.csv', '.chroms.csv'), header=None, index_col=0, squeeze=True)
    container.chr_names_dict = df.to_dict()
    return container


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--features', nargs='+')
    parser.add_argument('--features_fpaths', nargs='+')
    parser.add_argument('--labels', nargs='+')
    parser.add_argument('--reference_csv')
    parser.add_argument('--output_dirpath')

    args = parser.parse_args()
    # from quast_libs import search_references_meta
    # if search_references_meta.is_quast_first_run:
    #    coords_dirpath = join(coords_dirpath, 'raw')

    print_timestamp()
    print_info('Running Genome analyzer...')

    output_dirpath = args.output_dirpath
    if not isdir(output_dirpath):
        os.mkdir(output_dirpath)

    genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(args.reference_csv)

    result_fpath = join(output_dirpath, 'genome_info.txt')
    res_file = open(result_fpath, 'w')

    containers = []
    for feature, feature_fpath in zip(args.features, args.features_fpaths):
        containers.append(FeatureContainer([feature_fpath], feature))
    if not containers:
        print_notice('No file with genomic features were provided. '
                     'Use the --features option if you want to specify it.\n', indent='  ')
    for container in containers:
        if not container.fpaths:
            continue

        for fpath in container.fpaths:
            container.region_list += genes_parser.get_genes_from_file(fpath, container.kind)

        if len(container.region_list) == 0:
            print_warning('No genomic features of type "' + container.kind + '" were loaded.', indent='  ')
            res_file.write('Genomic features of type "' + container.kind + '" loaded: ' + 'None' + '\n')
        else:
            print_info('  Loaded ' + str(
                len(container.region_list)) + ' genomic features of type "' + container.kind + '"')
            res_file.write('Genomic features of type "' + container.kind + '" loaded: ' + str(
                len(container.region_list)) + '\n')
            container.chr_names_dict = chromosomes_names_dict(container.kind, container.region_list,
                                                              list(reference_chromosomes.keys()))
        pd.DataFrame.from_dict({'kind': [container.kind], 'region_list': [g.asdict() for g in container.region_list]},
                               orient='index').to_csv(join(output_dirpath, container.kind + '.csv'))
        pd.DataFrame.from_dict(container.chr_names_dict,
                               orient='index').to_csv(join(output_dirpath, container.kind + '.chroms.csv'), header=False)


if __name__ == '__main__':
    main()