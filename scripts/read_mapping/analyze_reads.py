import sys
from os.path import join

from src import qutils, qconfig
from src.common import get_labels_from_paths
from src.qutils import get_chr_len_fpath
from src.ra_utils.misc import bam_to_bed, calculate_genome_cov


def main():
    output_dirpath, ref_fpath, bam_fpath, max_threads = sys.argv[1:]
    label = get_labels_from_paths([ref_fpath])[0]
    bed_fpath = bam_to_bed(output_dirpath, label, bam_fpath)
    chr_len_fpath = get_chr_len_fpath(ref_fpath)
    stats_fpath = join(output_dirpath, label + '.stat')
    cov_fpath = join(output_dirpath, label + '.genomecov')
    calculate_genome_cov(bed_fpath, cov_fpath, chr_len_fpath, print_all_positions=False)

    qutils.call_subprocess(['samtools', 'flagstat', '-@', str(max_threads), bam_fpath],
                           stdout=open(stats_fpath, 'w'))
    avg_depth = 0
    coverage_for_thresholds = [0 for threshold in qconfig.coverage_thresholds]
    with open(cov_fpath) as f:
        for line in f:
            l = line.split()  # genome; depth; number of bases; size of genome; fraction of bases with depth
            depth, genome_fraction = int(l[1]), float(l[4])
            if l[0] == 'genome':
                avg_depth += depth * genome_fraction
                for i, threshold in enumerate(qconfig.coverage_thresholds):
                    if depth >= threshold:
                        coverage_for_thresholds[i] += genome_fraction

    with open(stats_fpath, 'a') as out_f:
        out_f.write('%s depth\n' % int(avg_depth))
        for i, threshold in enumerate(qconfig.coverage_thresholds):
            out_f.write('%.2f coverage >= %sx\n' % (coverage_for_thresholds[i] * 100, threshold))


if __name__ == '__main__':
    main()