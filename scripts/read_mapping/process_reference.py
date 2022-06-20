import os
from os.path import isdir, join

import pysam

from src import qutils, qconfig
from src.logger import *
from src.qutils import is_non_empty_file, add_suffix
from src.ra_utils.misc import *
from src.reads_analyzer import search_trivial_deletions, get_coverage, search_sv_with_gridss


def main():
    ref_fpath, bam_fpath, output_dir, max_threads, search_sv = sys.argv[1:]
    no_sv = not bool(search_sv)
    ref_name = "reference"
    temp_output_dir = join(output_dir, 'tmp')
    if not isdir(temp_output_dir):
        os.makedirs(temp_output_dir)

    bed_fpath = join(output_dir, ref_name + '.bed')
    cov_fpath = join(output_dir, ref_name + '.cov')
    physical_cov_fpath = join(output_dir, ref_name + '.physical.cov')
    required_files = []

    if no_sv:
        print_info('  Will not search Structural Variations ')
        bed_fpath = None
    else:
        required_files.append(bed_fpath)
    '''elif is_non_empty_file(bed_fpath):
        if not is_valid_bed(bed_fpath):
            print_warning('  Existing BED-file: ' + bed_fpath + ' may be corrupted. Bed file will be re-created. ')
            required_files.append(join(output_dir, ref_name + '.bed'))
        else: print_info('  Using existing BED-file: ' + bed_fpath)'''
    if qconfig.create_icarus_html:
        '''if is_non_empty_file(cov_fpath):
            is_correct_file = check_cov_file(cov_fpath)
            if is_correct_file:
                print_info('  Using existing reads coverage file: ' + cov_fpath)
            else:
                required_files.append(cov_fpath)
        if is_non_empty_file(physical_cov_fpath):
            print_info('  Using existing physical coverage file: ' + physical_cov_fpath)
        else:'''
        required_files.append(physical_cov_fpath)
    else:
        print_info('  Will not calculate coverage (--fast or --no-html, or --no-icarus, or --space-efficient is specified)')
        cov_fpath = None
        physical_cov_fpath = None

    bam_mapped_fpath = get_safe_fpath(temp_output_dir, add_suffix(bam_fpath, 'mapped'))
    bam_sorted_fpath = get_safe_fpath(temp_output_dir, add_suffix(bam_mapped_fpath, 'sorted'))

    convert_sam(bam_fpath, bam_mapped_fpath, max_threads, filter_rule='-F 4') #not unmapped
    sort_bam(bam_mapped_fpath, bam_sorted_fpath)
    if qconfig.create_icarus_html and (not is_non_empty_file(cov_fpath) or not is_non_empty_file(physical_cov_fpath)):
        correct_chr_names = get_correct_names_for_chroms(output_dir, ref_fpath, bam_fpath)
        cov_fpath, physical_cov_fpath = get_coverage(temp_output_dir, ref_fpath, ref_name, bam_fpath, bam_sorted_fpath,
                                                     correct_chr_names, max_threads, cov_fpath, physical_cov_fpath)
    if not no_sv:
        bamfile = pysam.AlignmentFile(bam_sorted_fpath, "rb")
        seq_names = bamfile.references
        seq_lengths = bamfile.lengths
        seq_lengths_dict = dict(zip(seq_names, seq_lengths))
        need_ref_splitting = False

        trivial_deletions_fpath = \
            search_trivial_deletions(temp_output_dir, bam_sorted_fpath, seq_lengths_dict, need_ref_splitting)
        try:
            gridss_sv_fpath = search_sv_with_gridss(ref_fpath, bam_mapped_fpath, temp_output_dir, max_threads)
            qutils.cat_files([gridss_sv_fpath, trivial_deletions_fpath], bed_fpath)
        except:
            pass
        if isfile(trivial_deletions_fpath) and not is_non_empty_file(bed_fpath):
            shutil.copy(trivial_deletions_fpath, bed_fpath)

    if not no_sv:
        if is_non_empty_file(bed_fpath):
            print_info('  Structural variations are in ' + bed_fpath)
        else:
            if isfile(bed_fpath):
                print_info('  No structural variations were found.')
            else:
                print_info('  Failed searching structural variations.')
    if is_non_empty_file(cov_fpath):
        print_info('  Coverage distribution along the reference genome is in ' + cov_fpath)
    else:
        if not qconfig.create_icarus_html:
            print_info('  Failed to calculate coverage distribution')


if __name__ == '__main__':
    main()