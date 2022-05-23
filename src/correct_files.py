import re
from collections import defaultdict

from Bio.Seq import Seq

from src.common import *
from src.logger import *
from src.qconfig import MAX_CONTIG_NAME


def correct_seq(seq, original_fpath):
    corr_seq = seq.upper()

    # correcting alternatives
    # dic = {'M': 'A', 'K': 'G', 'R': 'A', 'Y': 'C', 'W': 'A', 'S': 'C', 'V': 'A', 'B': 'C', 'H': 'A', 'D': 'A'}
    dic = {'M': 'N', 'K': 'N', 'R': 'N', 'Y': 'N', 'W': 'N', 'S': 'N', 'V': 'N', 'B': 'N', 'H': 'N', 'D': 'N'}
    pat = "(%s)" % "|".join(map(re.escape, dic.keys()))
    corr_seq = re.sub(pat, lambda m: dic[m.group()], corr_seq)

    # make sure that only A, C, G, T or N are in the sequence
    if re.compile(r'[^ACGTN]').search(corr_seq):
        print_error('Skipping ' + original_fpath + ' because it contains non-ACGTN characters.', indent='    ')
        return None
    return Seq(corr_seq)


def correct_name(name, max_name_len=MAX_CONTIG_NAME):
    name = re.sub(r'[^\w\._\-+|]', '_', name.strip())[:max_name_len]
    name = re.sub(r'[\.+]$', '', name)
    return re.sub(r"[\|\+\-=\/]", '_', name)


def correct_asm_label(name, max_name_len=MAX_CONTIG_NAME):
    return name.strip()[:max_name_len]


def get_uniq_name(name, used_names):
    if name in used_names:
        name += '_' + str(used_names[name])
    return name


def correct_fasta(original_fpath, min_contig=0, corrected_fpath=None, is_reference=False, no_check=False):
    modified_fasta_entries = []
    used_seq_names = defaultdict(int)
    for record in read_fasta(original_fpath):
        seq = str(record.seq)
        if not record.id:
            print_error('Skipping ' + original_fpath + ' because >sequence_name field is empty '
                                                        'for the entry starting with "%s".' % seq[:20], indent='    ')
            return False

        if (len(seq) >= min_contig) or is_reference:
            corr_name = correct_name(record.id)
            record.id = get_uniq_name(corr_name, used_seq_names)
            used_seq_names[corr_name] += 1

            if not no_check:
                # seq to uppercase, because we later looking only uppercase letters
                record.seq = correct_seq(seq, original_fpath)
                if not record.seq:
                    return False
            else:
                if re.compile(r'[^ACGTN]').search(seq):
                    print_error('File ' + original_fpath + ' contains non-ACGTN characters. '
                                    'Please re-run QUAST without --no-check.', indent='    ', exit_code=1)
                    return False
                record.seq = seq
            modified_fasta_entries.append(record)

    if not modified_fasta_entries:
        print_warning('Skipping ' + original_fpath + ' because file is empty.', indent='    ')
        return False
    if corrected_fpath:
        write_fasta(corrected_fpath, modified_fasta_entries)
    return True
