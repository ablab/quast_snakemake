import os
from os.path import join, isdir

configfile: "config.yaml"

corrected_dirpath = join(config['output_dir'], config['CORRECTED_DIRPATH'])

glimmer_dirpath = join(config['output_dir'], 'gene_prediction')
tmp_glimmer_dirpath = join(glimmer_dirpath, 'tmp')

glimmer_output = list()
if config['gene_prediction']:
    glimmer_output = expand(join(glimmer_dirpath, "{sample}_glimmer.gff"), sample=config['samples'])
    if not isdir(glimmer_output):
        os.makedirs(glimmer_output)


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

rule glimmer:
    input:
        contig=join(corrected_dirpath,"{sample}.fasta"),
    conda:
        "envs/basic.yaml"
    log:
        out = join(glimmer_dirpath, '{sample}.log'),
        err = join(glimmer_dirpath, '{sample}.err')
    params:
        label="{sample}",
        tmp_dir=tmp_glimmer_dirpath,
        output_dir=glimmer_dirpath,
    output:
        join(glimmer_dirpath, '{sample}_glimmer.gff')
    shell:
        "python -m scripts.gene_finding.glimmer {params.output_dir} {input.contig} {params.label} "
        "{params.tmp_dir} >{log.out} 2>{log.err}"

rule save_stats:
    input:
        contigs=expand(join(corrected_dirpath, "{sample}.fasta"), sample=config['samples']),
        glimmer_output=glimmer_output
    conda:
        "envs/basic.yaml"
    log:
        out=join(config['output_dir'], 'quast.log'),
        err=join(config['output_dir'], 'quast.err')
    params:
        tmp_glimmer_dirpath=tmp_glimmer_dirpath
    output:
        join(config['output_dir'], "report.txt")
    shell:
        "python -m scripts.make_reports -m {config[min_contig]} "
        "-o {config[output_dir]} --glimmer_dirpath {params.tmp_glimmer_dirpath} "
        "--contigs_fpaths {input.contigs} >{log.out} 2>{log.err}"

