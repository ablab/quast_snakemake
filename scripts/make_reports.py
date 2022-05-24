import argparse
from ast import literal_eval
from os.path import exists

from src.aligned_stats import parse_aligner_stats
from scripts.gene_finding.prepare_genome_analyzer import parse_results
from src import basic_stats
from src.assembly import Assembly
from src.common import *
from src.genes_parser import Gene
from src.html_saver import html_saver
from src.icarus import icarus
from src.logger import print_info, print_timestamp, print_error
from src.plotter_aux import dict_color_and_ls, save_colors_and_ls
from src.save_results import *


def parse_genome_stats(reports, reference_csv, assemblies, labels, output_dirpath, genome_analyzer_dirpath):
    # from quast_libs import search_references_meta
    # if search_references_meta.is_quast_first_run:
    #    coords_dirpath = join(coords_dirpath, 'raw')

    genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(reference_csv, skip_ns=True)

    result_fpath = join(genome_analyzer_dirpath, 'genome_info.txt')
    res_file = open(result_fpath, 'a')

    results = defaultdict()

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
            for container_kind in df.index:
                stats = df.loc[container_kind].to_list()
                if container_kind == 'operon':
                    results[reporting.Fields.OPERONS + "_full"] = stats[0]
                    results[reporting.Fields.OPERONS + "_partial"] = stats[1]
                    reports[label].add_field(reporting.Fields.REF_OPERONS, stats[2])
                else:
                    if reporting.Fields.GENES + "_full" not in results:
                        results[reporting.Fields.GENES + "_full"] = 0
                        results[reporting.Fields.GENES + "_partial"] = 0
                    results[reporting.Fields.GENES + "_full"] += stats[0]
                    results[reporting.Fields.GENES + "_partial"] += stats[1]
                    reports[label].add_field(reporting.Fields.REF_GENES, stats[2])

    if not results:
        print_info('Genome analyzer failed for all the assemblies.')
        res_file.close()
        return

    ref_genes_num, ref_operons_num = 0, 0
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

            gaps_count = df.loc["gaps_count"].dropna().to_list()[0]
            res_file.write('%-25s| %-10s| %-12s| %-10s|'
                           % (label[:24], reports[label].get_field(reporting.Fields.MAPPEDGENOME),
                              reports[label].get_field(reporting.Fields.DUPLICATION_RATIO), gaps_count))

            genome_mapped.append(float(reports[label].get_field(reporting.Fields.MAPPEDGENOME)))

            for (field, full, part) in [(reporting.Fields.GENES, results.get(reporting.Fields.GENES + "_full",None), results.get(reporting.Fields.GENES + "_partial",None)),
                                        (reporting.Fields.OPERONS, results.get(reporting.Fields.OPERONS + "_full",None), results.get(reporting.Fields.OPERONS + "_partial",None))]:
                if full is None and part is None:
                    res_file.write(' %-10s| %-10s|' % ('-', '-'))
                else:
                    res_file.write(' %-10s| %-10s|' % (full, part))
                    reports[label].add_field(field, '%s + %s part' % (full, part))
                    if field == reporting.Fields.OPERONS:
                        ref_operons_num = results[reporting.Fields.OPERONS + "_full"] + results[reporting.Fields.OPERONS + "_partial"]

                    if field == reporting.Fields.GENES:
                        ref_genes_num = results[reporting.Fields.GENES + "_full"] + results[reporting.Fields.GENES + "_partial"]

            res_file.write('\n')
    res_file.close()

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

            plotter.frc_plot(output_dirpath, reference_chromosomes, labels, [asm.aligned_lengths_by_contigs for asm in assemblies],
                             files_unsorted_features_in_contigs,
                             genome_analyzer_dirpath + '/features_frcurve_plot', 'genomic features')
            plotter.histogram(labels, full_found_genes,
                              genome_analyzer_dirpath + '/complete_features_histogram',
                              '# complete genomic features')
        if ref_operons_num:
            plotter.genes_operons_plot(ref_operons_num, labels, files_operons_in_contigs,
                                       genome_analyzer_dirpath + '/operons_cumulative_plot', 'operons')

            plotter.frc_plot(output_dirpath, reference_chromosomes, labels, [asm.aligned_lengths_by_contigs for asm in assemblies],
                             files_unsorted_operons_in_contigs,
                             genome_analyzer_dirpath + '/operons_frcurve_plot', 'operons')
            plotter.histogram(labels, full_found_operons,
                              genome_analyzer_dirpath + '/complete_operons_histogram',
                              '# complete operons')
        plotter.histogram(labels, genome_mapped,
                          genome_analyzer_dirpath + '/genome_fraction_histogram',
                          'Genome fraction, %', top_value=100)

    print_info('Done.')


def parse_glimmer(labels, reports, tmp_glimmer_dirpath):
    genes_by_labels = dict()
    # saving label
    for label in labels:
        glimmer_csv = join(tmp_glimmer_dirpath, label + '_glimmer.csv')
        if not exists(glimmer_csv): continue
        df = pd.read_csv(glimmer_csv, index_col=0)
        genes_list = df.loc['genes'].apply(literal_eval).dropna().to_list()
        genes_list = [Gene(**g) for g in genes_list]
        genes_by_labels[label] = genes_list
        full_genes = df.loc['full'].dropna().to_list()
        partial_genes = df.loc['partial'].dropna().to_list()
        unique = df.loc['unique'].dropna().to_list()
        if unique[0] is not None:
            reports[label].add_field(reporting.Fields.PREDICTED_GENES_UNIQUE, unique[0])
        if full_genes is not None:
            genes = ['%s + %s part' % (full_cnt, partial_cnt) for full_cnt, partial_cnt in zip(full_genes, partial_genes)]
            reports[label].add_field(reporting.Fields.PREDICTED_GENES, genes)
        if unique is None and full_genes is None:
            print_error(
                'Failed running Glimmer for %s. ' % label + ('Run with the --debug option'
                ' to see the command line.' if not qconfig.debug else ''))

    #if not qconfig.debug:
    #    shutil.rmtree(tmp_glimmer_dirpath)
    return genes_by_labels


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', dest='reference_csv')
    parser.add_argument('-o', dest='output_dirpath')
    parser.add_argument('-m', dest='min_contig')
    parser.add_argument('-r', dest='reference')
    parser.add_argument('--features', nargs='+')
    parser.add_argument('--contigs_fpaths', nargs='+')
    parser.add_argument('--contig_analyzer_dirpath')
    parser.add_argument('--genome_analyzer_dirpath')
    parser.add_argument('--glimmer_dirpath')
    args = parser.parse_args()

    output_dirpath = args.output_dirpath
    if isdir(output_dirpath):
        qutils.remove_reports(output_dirpath)

    qconfig.min_contig = int(args.min_contig)
    ref_fpath = args.reference
    contigs_fpaths = args.contigs_fpaths

    labels = get_labels_from_paths(contigs_fpaths)
    save_colors_and_ls(labels)
    assemblies = [Assembly(args.contig_analyzer_dirpath, contigs_fpath, label) for contigs_fpath, label in zip(contigs_fpaths, labels)]

    reports = dict((label, reporting.get(label)) for label in labels)

    reference_chromosomes = []
    if ref_fpath:
        genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(args.reference_csv, skip_ns=True)
        assemblies, successful_runs = parse_aligner_stats(reports, output_dirpath, assemblies, labels, ref_fpath, reference_chromosomes,
                                                      genome_size, join(output_dirpath, 'aligned_stats'))
        parse_genome_stats(reports, args.reference_csv, assemblies, labels, output_dirpath, args.genome_analyzer_dirpath)

    icarus_gc_fpath, circos_gc_fpath = basic_stats.do(reports, ref_fpath, reference_chromosomes, assemblies, output_dirpath, join(output_dirpath, 'basic_stats'))

    features_containers = [parse_results(c) for c in args.features] if args.features else []

    genes_by_labels = parse_glimmer(labels, reports, args.glimmer_dirpath)

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
        try:
            if args.contig_analyzer_dirpath:
                detailed_contigs_reports_dirpath = os.path.join(args.contig_analyzer_dirpath,
                                                                qconfig.detailed_contigs_reports_dirname)
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
                    stdout_pattern=stdout_pattern, gc_fpath=icarus_gc_fpath, json_output_dir=qconfig.json_output_dirpath,
                    features=features_containers, genes_by_labels=genes_by_labels)
                    #cov_fpath=cov_fpath, physical_cov_fpath=physical_cov_fpath,

            if draw_circos_plot:
                print_info('  %d of %d: Creating Circos plot...' % (number_of_steps, number_of_steps))
                from src import circos
                circos_png_fpath, circos_legend_fpath = circos.do(ref_fpath, contigs_fpaths, labels,
                                                                  report_for_icarus_fpath_pattern, circos_gc_fpath,
                                                                  features_containers,
                                                                  os.path.join(output_dirpath, 'circos'))

            print_info('Done')
        except KeyboardInterrupt:
            print_info('..step skipped!')


if __name__ == '__main__':
    main()