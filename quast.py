#!/usr/bin/env python3

############################################################################
# Copyright (c) 2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import snakemake
import os
import sys
import yaml
import argparse
import logging
from src import qutils
from src import qconfig


logger = logging.getLogger("quast")
config_name = "config.yaml"
config_template = os.path.join(os.path.dirname(__file__), config_name)
snakemake_file_general = os.path.join(os.path.dirname(__file__), "quast.smk")
snakemake_file_no_ref = os.path.join(os.path.dirname(__file__), "quast_no_ref.smk")


def process_path(path):
    return os.path.abspath(os.path.expanduser(path))


def proper_quote(arg):
    import os
    if os.name == 'nt':
        from subprocess import list2cmdline
        return list2cmdline([arg])
    else:
        import shlex
        return shlex.quote(arg)


def quote_init(self, quote_func=proper_quote, *args, **kwargs):
    self.quote_func = quote_func
    super(type(self), self).__init__(*args, **kwargs)


def init_logger(args):
    import sys

    mlogger = logging.getLogger()
    mlogger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    ch.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    mlogger.addHandler(ch)


def parse_command_line(description="Snakemake-based QUAST: Quality Assessment Tool for Genome Assemblies"):
    parser = argparse.ArgumentParser(description=description, add_help=False)

    positional = parser.add_argument_group('Positional argument (required)')
    positional.add_argument("contigs_fpaths",
                            nargs="+",
                            type=str,
                            metavar="ASSEMBLY_FILE",
                            # action=,  # TODO: do some pre-check/preprocessing in the callback, e.g. make the path absolute
                            help="Assembly file(s) in the FASTA format.")

    main_options = parser.add_argument_group('Main options')
    main_options.add_argument("--output-dir", "-o",
                              type=str,
                              metavar="DIR",
                              # action=,  # TODO: do some pre-check/preprocessing in the callback, e.g. make the path absolute
                              help="Output directory.")
    main_options.add_argument("--reference", "-r",
                              type=str,
                              metavar="FILE",
                              default=None,
                              # action=,  # TODO: do some pre-check/preprocessing in the callback, e.g. make the path absolute
                              help="Reference genome file in the FASTA format.")
    main_options.add_argument("--features", "-g",
                              type=str,
                              metavar="FILE",
                              action="append",
                              # action=,  # TODO: do some pre-check/preprocessing in the callback, e.g. make the path absolute
                              help="File with genomic feature coordinates in the reference (GFF, BED, NCBI or TXT). "
                                   "Could be specified in the '[type:]<filepath>' format where 'type' limits "
                                   "the search to features of a specific type only, e.g., 'CDS'.")

    main_options.add_argument("--min-contig", "-m",
                              type=int,
                              metavar="N",
                              default=500,
                              # action=,  # TODO: check it is >= 0
                              help="Lower threshold for contig/scaffold length.")
    main_options.add_argument("--threads", "-t", "--cores", "--jobs", "-j",
                              action="store",
                              dest="threads",
                              const=snakemake.utils.available_cpu_count(),
                              nargs="?",
                              metavar="N",
                              type=int,
                              help="Use at most N cores in parallel (default: 1). "
                                   "If N is omitted, the limit is set to the number of "
                                   "available cores.")

    snakemake_log_options = parser.add_argument_group('Advanced: Snakemake logging options')
    snakemake_log_options.add_argument("--verbose",
                                       action="store_true",
                                       help="Print debugging output.")
    snakemake_log_options.add_argument("--stats",
                                       metavar="FILE",
                                       help="Write pipeline execution statistics to the JSON file.")
    snakemake_log_options.add_argument("--nocolor",
                                       action="store_true",
                                       help="Do not use a colored output.")
    snakemake_log_options.add_argument("--quiet", "-q",
                                       action="store_true",
                                       help="Do not output any progress or rule information.")
    snakemake_log_options.add_argument("--printshellcmds", "-p",
                                       action="store_true",
                                       help="Print out the shell commands that will be executed.")

    snakemake_run_options = parser.add_argument_group('Advanced: Snakemake execution options')
    snakemake_run_options.add_argument("--forceall", "-F",
                                       action="store_true",
                                       help="Force the execution of the pipeline regardless of already created output.")
    snakemake_run_options.add_argument("--forcerun", "-R",
                                       nargs="*", help=argparse.SUPPRESS)  # TODO: add explicit help message
    snakemake_run_options.add_argument("--target",
                                       nargs="*",
                                       default=None,
                                       help=argparse.SUPPRESS)  # TODO: add explicit help message

    other_options = parser.add_argument_group('Other options')
    other_options.add_argument("--help", "-h",
                               action='help',
                               default=argparse.SUPPRESS,
                               help="Print help message and exit.")
    other_options.add_argument("--version", "-v",
                               action="version",
                               version="6.0.0",  # TODO: store version somewhere else (version.py? config.yaml?) and read from there
                               default=argparse.SUPPRESS,
                               help="Print version and exit.")

    args = parser.parse_args()

    # TODO: various checks and assertions, e.g., existence of contigs_fpaths
    return args


def prepare_config(args):
    with open(config_template) as cfg:
        try:
            config = yaml.safe_load(cfg)
        except yaml.YAMLError as exc:
            logger.error(exc)
            return 1

    config["output_dir"] = process_path(args.output_dir)

    # TODO: implement more advanced assembly labels options: user-provided and "all_labels_from_dirs"
    contigs_labels = qutils.process_labels(args.contigs_fpaths)
    config["samples"] = dict(zip(contigs_labels, map(process_path, args.contigs_fpaths)))

    if args.reference is not None:
        config["reference"] = process_path(args.reference)

    config["features"] = []
    config["features_files"] = []
    if args.features:
        for value in args.features:
            if ':' in value:
                feature, fpath = value.split(':')
            else:
                feature, fpath = qconfig.ALL_FEATURES_TYPE, value  # special case -- read all features
            config["features"].append(feature)
            config["features_files"].append(process_path(fpath))

    config["min_contig"] = args.min_contig

    with open(os.path.join(args.output_dir, config_name), 'w') as dst:
        yaml.dump(config, dst)

    return config


def main():
    args = parse_command_line()

    init_logger(args)  # TODO: organise logging similar to the conventional QUAST

    if args.threads is None:  # TODO: use the "default number of threads" strategy from the conventional QUAST
        args.threads = 1

    logger.info("QUAST started!")

    command_line = " ".join(sys.argv)
    logger.info("Command line: " + command_line)

    snakemake.utils.QuotedFormatter.__init__ = quote_init
    snakemake.utils.makedirs(args.output_dir)  # TODO: if output_dir is not specified, use "quast_results/XXX" as in the conventional QUAST

    logger.debug("Args:")
    logger.debug(str(args))

    config = prepare_config(args)

    rc = snakemake.snakemake(snakefile=snakemake_file_no_ref if args.reference is None else snakemake_file_general,
                             workdir=os.path.dirname(__file__),
                             cores=args.threads,
                             forceall=args.forceall,
                             verbose=args.verbose,
                             stats=args.stats,
                             nocolor=args.nocolor,
                             quiet=args.quiet,
                             targets=args.target,
                             forcerun=args.forcerun,
                             printshellcmds=args.printshellcmds,
                             config=config)

    return 0 if rc else 1


if __name__ == "__main__":
    rc = main()
    sys.exit(rc)