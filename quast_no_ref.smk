import os
from os.path import join, isdir

from scripts.gene_finding.run_busco import get_lineage
from src.qutils import get_path_to_program, get_dir_for_download

configfile: "config.yaml"

corrected_dirpath = join(config['output_dir'], config['CORRECTED_DIRPATH'])


gene_pred_dirpath = join(config['output_dir'], 'gene_prediction')
tmp_gene_pred_dirpath = join(gene_pred_dirpath, 'tmp')
jobs_threads = max(1, config['threads']//len(config['samples']))

gene_pred_output = list()
if config['gene_prediction']:
    gene_pred_output = expand(join(gene_pred_dirpath, "{sample}.gff"), sample=config['samples'])
    if not isdir(tmp_gene_pred_dirpath):
        os.makedirs(tmp_gene_pred_dirpath)
else:
    tmp_gene_pred_dirpath = None

busco_dirpath = join(config['output_dir'], 'busco')
lineage = get_lineage(is_prokaryote=config['is_prokaryote'],is_fungus=config['is_fungus'])
if config['busco'] and get_path_to_program('busco'):
    busco_output = expand(join(busco_dirpath, "{sample}/short_summary.specific." + lineage + ".{sample}.txt"), sample=config['samples'])
else:
    busco_output = list()
    busco_dirpath = None

rule all:
    input:
        join(config['output_dir'], "report.txt")

def get_input_fastas(wildcards):
    return config["samples"][wildcards.sample]

rule correct_contigs:
    input:
        get_input_fastas
    log:
        out=join(corrected_dirpath,"{sample}.log"),
        err=join(corrected_dirpath,"{sample}.err")
    output:
        join(corrected_dirpath,"{sample}.fasta")
    shell:
        "python -m scripts.preprocessing.correct_contig {input} {output} {config[min_contig]} >{log.out} 2>{log.err}"

rule gene_prediction:
    input:
        contig=join(corrected_dirpath,"{sample}.fasta"),
    conda:
        "envs/basic.yaml"
    log:
        out = join(gene_pred_dirpath, '{sample}.log'),
        err = join(gene_pred_dirpath, '{sample}.err')
    params:
        label="{sample}",
        tmp_dir=tmp_gene_pred_dirpath,
        output_dir=gene_pred_dirpath,
        tool = 'prodigal' if config['is_prokaryote'] else 'glimmerhmm'
    output:
        join(gene_pred_dirpath, '{sample}.gff')
    shell:
        "python -m scripts.gene_finding.gene_prediction {params.output_dir} {input.contig} {params.label} "
        "{params.tmp_dir} {params.tool} >{log.out} 2>{log.err}"

if busco_dirpath:
    db_dirpath = join(get_dir_for_download('busco', 'busco_db', [lineage]), lineage)
    rule download_busco:
        conda:
            "envs/busco.yaml"
        log:
            out=join(busco_dirpath,'busco.log'),
            err=join(busco_dirpath,'busco.err')
        params:
            output_dir=busco_dirpath
        output:
            directory(db_dirpath)
        shell:
            "python -m scripts.gene_finding.busco_download >{log.out} 2>{log.err}"

    rule run_busco:
        input:
            database=db_dirpath,
            contig=join(corrected_dirpath,"{sample}.fasta"),
        conda:
            "envs/busco.yaml"
        log:
            out=join(busco_dirpath,'{sample}.log'),
            err=join(busco_dirpath,'{sample}.err')
        params:
            label="{sample}",
            output_dir=busco_dirpath
        output:
            join(busco_dirpath, "{sample}/short_summary.specific." + lineage + ".{sample}.txt")
        shell:
            "python -m scripts.gene_finding.run_busco {params.label} {input.contig} {params.output_dir} "
            "{input.database} {jobs_threads} >{log.out} 2>{log.err}"

rule save_stats:
    input:
        contigs=expand(join(corrected_dirpath, "{sample}.fasta"), sample=config['samples']),
        gene_pred_output=gene_pred_output,
        busco_output=busco_output
    conda:
        "envs/basic.yaml"
    log:
        out=join(config['output_dir'], 'quast.log'),
        err=join(config['output_dir'], 'quast.err')
    params:
        busco_dirpath=busco_dirpath,
        lineage=lineage,
        tmp_gene_pred_dirpath=tmp_gene_pred_dirpath
    output:
        join(config['output_dir'], "report.txt")
    shell:
        "python -m scripts.make_reports -m {config[min_contig]} "
        "-o {config[output_dir]} --gene_pred_dirpath {params.tmp_gene_pred_dirpath} "
        "--busco_dirpath {params.busco_dirpath} --lineage {params.lineage} --contigs_fpaths {input.contigs} >{log.out} 2>{log.err}"

