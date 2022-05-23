############################################################################
# Copyright (c) 2015-2020 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os


######## MAIN ############
from os.path import join, dirname, isdir

from src import plotter, qconfig, reporting, N50
from src.align_contigs import AlignerStatus
from src.common import get_chr_lengths_from_fastafile
from src.logger import print_info, print_timestamp
from src.save_results import save_result_for_unaligned, save_result


def parse_aligner_stats(reports, output_dirpath, assemblies, labels, ref_fpath, reference_chromosomes, genome_size, aligned_stats_dirpath):
    if not isdir(aligned_stats_dirpath):
        os.makedirs(aligned_stats_dirpath)

    for asm in assemblies:
        asm.parse_results()
    aligner_statuses = [asm.status for asm in assemblies]
    '''if AlignerStatus.OK in aligner_statuses:
        if qconfig.is_combined_ref:
            save_combined_ref_stats(results, contigs_fpaths, ref_labels_by_chromosomes, output_dir, logger)'''

    saved_reports = []
    for index, label in enumerate(labels):
        if assemblies[index].status == AlignerStatus.OK:
            saved_reports.append(save_result(assemblies[index].results, reports[label], label, ref_fpath, genome_size))
        elif assemblies[index].status == AlignerStatus.NOT_ALIGNED:
            save_result_for_unaligned(assemblies[index].results, reports[label])

    if AlignerStatus.OK in aligner_statuses:
        reporting.save_misassemblies(reports, output_dirpath)
        reporting.save_unaligned(reports, output_dirpath)
        if qconfig.draw_plots:
            plotter.draw_misassemblies_plot(saved_reports, join(output_dirpath, 'misassemblies_plot'), 'Misassemblies')
        if qconfig.draw_plots or qconfig.html_report:
            misassemblies_in_contigs = dict(
                (labels[i], assemblies[i].misassemblies_in_contigs) for i in range(len(labels)))
            plotter.frc_plot(output_dirpath, reference_chromosomes, labels,
                             [asm.aligned_lengths_by_contigs for asm in assemblies],
                             misassemblies_in_contigs,
                             join(aligned_stats_dirpath, 'misassemblies_frcurve_plot'), 'misassemblies')

    oks = aligner_statuses.count(AlignerStatus.OK)
    not_aligned = aligner_statuses.count(AlignerStatus.NOT_ALIGNED)
    failed = aligner_statuses.count(AlignerStatus.FAILED)
    errors = aligner_statuses.count(AlignerStatus.ERROR)
    problems = not_aligned + failed + errors
    all = len(aligner_statuses)

    if oks == all:
        print_info('Done.')
    if oks < all and problems < all:
        print_info('Done for ' + str(all - problems) + ' out of ' + str(
            all) + '. For the rest, only basic stats are going to be evaluated.')
    if problems == all:
        print_info(
            'Failed aligning the contigs for all the assemblies. Only basic stats are going to be evaluated.')
        return assemblies, aligner_statuses.count(AlignerStatus.OK)

    print_timestamp()
    print_info('Running NA-NGA calculation...')

    reference_length = sum(reference_chromosomes.values())
    assembly_lengths = []
    for asm in assemblies:
        assembly_lengths.append(sum(get_chr_lengths_from_fastafile(asm.contigs_fpath).values()))

    for i, (asm, assembly_len) in enumerate(
            zip(assemblies, assembly_lengths)):
        lens = asm.aligned_lengths_by_contigs
        sorted_lengths = sorted(lens, reverse=True)
        na50, la50 = N50.NG50_and_LG50(sorted_lengths, assembly_len)
        auNA = N50.au_metric(sorted_lengths, assembly_len)
        nax, lax = N50.NG50_and_LG50(sorted_lengths, assembly_len, qconfig.x_for_additional_Nx)
        if not qconfig.is_combined_ref:
            nga50, lga50 = N50.NG50_and_LG50(sorted_lengths, reference_length)
            ngax, lgax = N50.NG50_and_LG50(sorted_lengths, reference_length, qconfig.x_for_additional_Nx)
            auNGA = N50.au_metric(sorted_lengths, reference_length)

        label = asm.label
        print_info('  ' + label +
                 ', Largest alignment = ' + str(max(lens)) +
                 ', NA50 = ' + str(na50) +
                 (', NGA50 = ' + str(nga50) if not qconfig.is_combined_ref and nga50 else '') +
                 ', LA50 = ' + str(la50) +
                 (', LGA50 = ' + str(lga50) if not qconfig.is_combined_ref and lga50 else ''))
        reports[label].add_field(reporting.Fields.LARGALIGN, max(lens))
        reports[label].add_field(reporting.Fields.TOTAL_ALIGNED_LEN, sum(lens))
        reports[label].add_field(reporting.Fields.NA50, na50)
        reports[label].add_field(reporting.Fields.NAx, nax)
        reports[label].add_field(reporting.Fields.auNA, auNA)
        reports[label].add_field(reporting.Fields.LA50, la50)
        reports[label].add_field(reporting.Fields.LAx, lax)
        if not qconfig.is_combined_ref:
            reports[label].add_field(reporting.Fields.NGA50, nga50)
            reports[label].add_field(reporting.Fields.NGAx, ngax)
            reports[label].add_field(reporting.Fields.auNGA, auNGA)
            reports[label].add_field(reporting.Fields.LGA50, lga50)
            reports[label].add_field(reporting.Fields.LGAx, lgax)

    ########################################################################
    num_contigs = max([len(asm.aligned_lengths_by_contigs) for asm in assemblies])

    # saving to html
    if qconfig.html_report:
        from src.html_saver import html_saver
        html_saver.save_assembly_lengths(output_dirpath, labels, assembly_lengths)

    aligned_lengths_lists = [asm.aligned_lengths_by_contigs for asm in assemblies]
    if qconfig.draw_plots:
        # Drawing cumulative plot (aligned contigs)...
        plotter.cumulative_plot(reference_chromosomes, labels, aligned_lengths_lists,
                                os.path.join(aligned_stats_dirpath, 'cumulative_plot'),
                                'Cumulative length (aligned contigs)')

        # Drawing NAx and NGAx plots...
    plotter.Nx_plot(output_dirpath, num_contigs > qconfig.max_points, labels, aligned_lengths_lists, aligned_stats_dirpath + '/NAx_plot', 'NAx',
                    assembly_lengths)
    if not qconfig.is_combined_ref:
        plotter.Nx_plot(output_dirpath, num_contigs > qconfig.max_points, labels, aligned_lengths_lists,
                        aligned_stats_dirpath + '/NGAx_plot', 'NGAx', [reference_length for i in range(len(labels))])

    print_info('Done.')
    return assemblies, aligner_statuses.count(AlignerStatus.OK)
