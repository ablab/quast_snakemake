############################################################################
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from os.path import join

import shutil

from src import qutils
from src.qutils import get_path_to_program


def get_lineage(is_prokaryote, is_fungus=False):
    if is_prokaryote:
        return 'bacteria_odb10'
    elif is_fungus:
        return 'fungi_odb10'
    else:
        return 'eukaryota_odb10'


def cleanup(busco_output_dir):
    do_not_remove_exts = ['.log', '.txt']
    for dirpath, dirnames, files in os.walk(busco_output_dir):
        for dirname in dirnames:
            shutil.rmtree(join(dirpath, dirname))
        for filename in files:
            if os.path.splitext(filename)[1] not in do_not_remove_exts:
                os.remove(join(dirpath, filename))
        break


def main():
    label, contigs_fpath, out_dirpath, lineage, threads = sys.argv[1:]
    cmdline = [get_path_to_program('busco'), '-i', contigs_fpath, '-m', 'genome', '--offline',
               '-o', label, '--out_path', out_dirpath, '-l', lineage, '-c', threads, '-f']
    qutils.call_subprocess(cmdline)


if __name__ == '__main__':
    main()
