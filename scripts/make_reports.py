import sys
from os.path import dirname, exists

from src import basic_stats, reporting
from src.align_contigs import AlignerStatus
from src.assembly import Assembly
from src.common import *
from src.html_saver import html_saver
from src.icarus import icarus
from src.logger import print_info, print_timestamp
from src.plotter_aux import dict_color_and_ls, save_colors_and_ls
from src.save_results import *


def parse_genome_stats(reports, reference_csv, labels, output_dirpath, genome_analyzer_dirpath):
    # from quast_libs import search_references_meta
    # if search_references_meta.is_quast_first_run:
    #    coords_dirpath = join(coords_dirpath, 'raw')

    genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(reference_csv, skip_ns=True)

    result_fpath = join(genome_analyzer_dirpath, 'genome_info.txt')
    res_file = open(result_fpath, 'a')

    results = defaultdict(int)

    # for cumulative plots:
    files_features_in_contigs = {}  # "filename" : [ genes in sorted contigs (see below) ]
    files_unsorted_features_in_contigs = {}  # "filename" : [ genes in sorted contigs (see below) ]
    files_operons_in_contigs = {}
    files_unsorted_operons_in_contigs = {}

    # for histograms
    genome_mapped = []
    full_found_genes = []
    full_found_operons = []

    res_file.write('reference chromosomes:\n')
    for chr_name, chr_len in reference_chromosomes.items():
        res_file.write('\t' + chr_name + ' (total length: ' + str(chr_len) + ' bp, ' +
                       'total length without N\'s: ' + str(chr_len - len(ns_by_chromosomes[chr_name])) +
                       ' bp)\n')
    res_file.write('\n')
    res_file.write('total genome size: ' + str(genome_size) + '\n\n')
    res_file.write('gap min size: ' + str(qconfig.min_gap_size) + '\n')
    res_file.write('partial gene/operon min size: ' + str(qconfig.min_gene_overlap) + '\n\n')
    # header
    # header
    res_file.write('\n\n')
    res_file.write('%-25s| %-10s| %-12s| %-10s| %-10s| %-10s| %-10s| %-10s|\n'
                   % ('assembly', 'genome', 'duplication', 'gaps', 'genes', 'partial', 'operons', 'partial'))
    res_file.write('%-25s| %-10s| %-12s| %-10s| %-10s| %-10s| %-10s| %-10s|\n'
                   % ('', 'fraction', 'ratio', 'number', '', 'genes', '', 'operons'))
    res_file.write('=' * 120 + '\n')

    for label in labels:
        result_fpath = join(genome_analyzer_dirpath, label + '_info.txt')
        if exists(result_fpath):
            df = pd.read_csv(result_fpath, index_col=0)
            # TODO: add non-gene features
            gene_stats = df.loc['gene'].to_list()
            operon_stats = df.loc['operon'].to_list()
            results[reporting.Fields.OPERONS + "_full"] = operon_stats[0]
            results[reporting.Fields.OPERONS + "_partial"] = operon_stats[1]
            results[reporting.Fields.GENES + "_full"] += gene_stats[0]
            results[reporting.Fields.GENES + "_partial"] += gene_stats[1]
            reports[label].add_field(reporting.Fields.REF_OPERONS, operon_stats[2])
            reports[label].add_field(reporting.Fields.REF_GENES, gene_stats[2])

    if not results:
        print_info('Genome analyzer failed for all the assemblies.')
        res_file.close()
        return

    for i, label in enumerate(labels):
        contig_result_fpath = join(genome_analyzer_dirpath, label + '_contig_info.txt')
        if exists(result_fpath):
            df = pd.read_csv(contig_result_fpath, index_col=0)

            files_features_in_contigs[label] = df.loc['features'].dropna().to_list()
            files_unsorted_features_in_contigs[label] = df.loc['features_unsorted'].dropna().to_list()
            files_operons_in_contigs[label] = df.loc['operons'].dropna().to_list()
            files_unsorted_operons_in_contigs[label] = df.loc['operons_unsorted'].dropna().to_list()
            full_found_genes.append(sum(files_features_in_contigs[label]))
            full_found_operons.append(sum(files_operons_in_contigs[label]))

            gaps_count = results["gaps_count"]
            res_file.write('%-25s| %-10s| %-12s| %-10s|'
                           % (label[:24], reports[label].get_field(reporting.Fields.MAPPEDGENOME),
                              reports[label].get_field(reporting.Fields.DUPLICATION_RATIO), gaps_count))

            genome_mapped.append(float(reports[label].get_field(reporting.Fields.MAPPEDGENOME)))

            for (field, full, part) in [(reporting.Fields.GENES, results[reporting.Fields.GENES + "_full"], results[reporting.Fields.GENES + "_partial"]),
                                        (reporting.Fields.OPERONS, results[reporting.Fields.OPERONS + "_full"], results[reporting.Fields.OPERONS + "_partial"])]:
                if full is None and part is None:
                    res_file.write(' %-10s| %-10s|' % ('-', '-'))
                else:
                    res_file.write(' %-10s| %-10s|' % (full, part))
                    reports[label].add_field(field, '%s + %s part' % (full, part))
            res_file.write('\n')
    res_file.close()

    ref_genes_num = results[reporting.Fields.GENES + "_full"] + results[reporting.Fields.GENES + "_partial"]
    ref_operons_num = results[reporting.Fields.OPERONS + "_full"] + results[reporting.Fields.OPERONS + "_partial"]
    if qconfig.html_report:
        if ref_genes_num:
            html_saver.save_features_in_contigs(output_dirpath, labels, 'features',
                                                files_features_in_contigs, ref_genes_num)
        if ref_operons_num:
            html_saver.save_features_in_contigs(output_dirpath, labels, 'operons',
                                                files_operons_in_contigs, ref_operons_num)

    if qconfig.draw_plots:
        # cumulative plots:
        if ref_genes_num:
            plotter.genes_operons_plot(ref_genes_num, labels, files_features_in_contigs,
                                       genome_analyzer_dirpath + '/features_cumulative_plot', 'genomic features')
            # TODO: fix
            # plotter.frc_plot(output_dirpath, ref_fpath, aligned_contigs_fpaths, contigs_aligned_lengths,
            #                 files_unsorted_features_in_contigs,
            #                 genome_stats_dirpath + '/features_frcurve_plot', 'genomic features')
            plotter.histogram(labels, full_found_genes,
                              genome_analyzer_dirpath + '/complete_features_histogram',
                              '# complete genomic features')
        if ref_operons_num:
            plotter.genes_operons_plot(ref_operons_num, labels, files_operons_in_contigs,
                                       genome_analyzer_dirpath + '/operons_cumulative_plot', 'operons')
            # TODO: fix
            # plotter.frc_plot(output_dirpath, ref_fpath, aligned_contigs_fpaths, contigs_aligned_lengths,
            #                 files_unsorted_operons_in_contigs,
            #                 genome_stats_dirpath + '/operons_frcurve_plot', 'operons')
            plotter.histogram(labels, full_found_operons,
                              genome_analyzer_dirpath + '/complete_operons_histogram',
                              '# complete operons')
        plotter.histogram(labels, genome_mapped,
                          genome_analyzer_dirpath + '/genome_fraction_histogram',
                          'Genome fraction, %', top_value=100)

    print_info('Done.')
    return reports


def parse_aligner_stats(reports, output_dirpath, assemblies, labels, ref_fpath, genome_size):
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
            plotter.frc_plot(dirname(output_dirpath), ref_fpath, labels,
                             [asm.aligned_lengths_by_contigs for asm in assemblies],
                             misassemblies_in_contigs,
                             join(output_dirpath, 'misassemblies_frcurve_plot'), 'misassemblies')

    oks = aligner_statuses.count(AlignerStatus.OK)
    not_aligned = aligner_statuses.count(AlignerStatus.NOT_ALIGNED)
    failed = aligner_statuses.count(AlignerStatus.FAILED)
    errors = aligner_statuses.count(AlignerStatus.ERROR)
    problems = not_aligned + failed + errors
    all = len(aligner_statuses)

    # logger._num_nf_errors = num_nf_errors + errors

    if oks == all:
        print_info('Done.')
    if oks < all and problems < all:
        print_info('Done for ' + str(all - problems) + ' out of ' + str(
            all) + '. For the rest, only basic stats are going to be evaluated.')
    if problems == all:
        print_info(
            'Failed aligning the contigs for all the assemblies. Only basic stats are going to be evaluated.')
    return assemblies, aligner_statuses.count(AlignerStatus.OK)


def main():
    min_contig, reference_csv, ref_fpath, output_dirpath, \
    contig_analyzer_dirpath, genome_analyzer_dirpath = sys.argv[1:7]
    contigs_fpaths = sys.argv[7:]

    if isdir(output_dirpath):
        qutils.remove_reports(output_dirpath)

    qconfig.min_contig = int(min_contig)
    genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(reference_csv, skip_ns=True)

    labels = get_labels_from_paths(contigs_fpaths)
    save_colors_and_ls(labels)
    assemblies = [Assembly(contig_analyzer_dirpath, contigs_fpath, label) for contigs_fpath, label in zip(contigs_fpaths, labels)]

    reports = dict((label, reporting.get(label)) for label in labels)

    icarus_gc_fpath, circos_gc_fpath = basic_stats.do(reports, ref_fpath, assemblies, output_dirpath, join(output_dirpath, 'basic_stats'))

    assemblies, successful_runs = parse_aligner_stats(reports, output_dirpath, assemblies, labels, ref_fpath, genome_size)
    reports = parse_genome_stats(reports, reference_csv, labels, output_dirpath, genome_analyzer_dirpath)

    html_saver.save_colors(output_dirpath, labels, dict_color_and_ls)
    html_saver.save_total_report(reports, output_dirpath, labels, qconfig.min_contig, ref_fpath)

    reports_fpaths, transposed_reports_fpaths = reporting.save_total(reports, output_dirpath)
    all_pdf_fpath = join(output_dirpath, "report.pdf")

    ########################################################################
    ### LARGE DRAWING TASKS
    ########################################################################
    if qconfig.draw_plots or qconfig.create_icarus_html or qconfig.draw_circos:
        print_timestamp()
        print_info('Creating large visual summaries...')
        print_info('This may take a while: press Ctrl-C to skip this step..')
        detailed_contigs_reports_dirpath = os.path.join(contig_analyzer_dirpath, qconfig.detailed_contigs_reports_dirname)
        try:
            if detailed_contigs_reports_dirpath:
                report_for_icarus_fpath_pattern = os.path.join(detailed_contigs_reports_dirpath,
                                                               qconfig.icarus_report_fname_pattern)
                stdout_pattern = os.path.join(detailed_contigs_reports_dirpath, qconfig.contig_report_fname_pattern)
            else:
                report_for_icarus_fpath_pattern = None
                stdout_pattern = None
            draw_alignment_plots = qconfig.create_icarus_html
            draw_circos_plot = qconfig.draw_circos and ref_fpath and successful_runs
            number_of_steps = sum(
                [int(bool(value)) for value in [draw_alignment_plots, draw_circos_plot, all_pdf_fpath]])
            if qconfig.draw_plots:
                # full report in PDF format: all tables and plots
                print_info(
                    '  1 of %d: Creating PDF with all tables and plots...' % number_of_steps)
                plotter.fill_all_pdf_file(all_pdf_fpath)

            if draw_alignment_plots:
                ########################################################################
                ### VISUALIZE CONTIG ALIGNMENT
                ########################################################################
                print_info(
                    '  %d of %d: Creating Icarus viewers...' % (2 if all_pdf_fpath else 1, number_of_steps))

                icarus_html_fpath = icarus.do(reports,
                    contigs_fpaths, report_for_icarus_fpath_pattern, output_dirpath, ref_fpath,
                    stdout_pattern=stdout_pattern, gc_fpath=icarus_gc_fpath, json_output_dir=qconfig.json_output_dirpath)
                    #features=features_containers, genes_by_labels=genes_by_labels)
                    #cov_fpath=cov_fpath, physical_cov_fpath=physical_cov_fpath,

            '''if draw_circos_plot:
                print_info('  %d of %d: Creating Circos plot...' % (number_of_steps, number_of_steps))
                from quast_libs import circos
                circos_png_fpath, circos_legend_fpath = circos.do(ref_fpath, contigs_fpaths,
                                                                  report_for_icarus_fpath_pattern, circos_gc_fpath,
                                                                  features_containers, cov_fpath,
                                                                  os.path.join(output_dirpath, 'circos'), logger)'''

            print_info('Done')
        except KeyboardInterrupt:
            print_info('..step skipped!')


if __name__ == '__main__':
    main()