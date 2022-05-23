import sys

from src.contigs_analyzer import align_and_analyze
from src.assembly import Assembly
from src.common import parse_ref_stats


def main():
    ref_fpath, asm_label, contigs_fpath, output_dirpath, is_cyclic, reference_csv, threads = sys.argv[1:]
    is_cyclic = True if is_cyclic == 'True' else False
    genome_size, reference_chromosomes, ns_by_chromosomes = parse_ref_stats(reference_csv, skip_ns=True)
    status, results, aligned_lengths, misassemblies_in_contigs, aligned_lengths_by_contigs = \
        align_and_analyze(ref_fpath, contigs_fpath, asm_label, output_dirpath, is_cyclic,
                          reference_chromosomes, ns_by_chromosomes, threads=int(threads))
    asm = Assembly(output_dirpath, contigs_fpath, asm_label)
    asm.add_results(status, results, aligned_lengths, misassemblies_in_contigs, aligned_lengths_by_contigs)
    asm.save_results()


if __name__ == '__main__':
    main()
