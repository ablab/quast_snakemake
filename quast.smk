import os
from os.path import join, isdir

configfile: "config.yaml"

contig_analyzer_dirpath = join(config['output_dir'], 'contig_analyzer')
genome_analyzer_dirpath = join(config['output_dir'], 'genome_analyzer')
minimap_dirpath = join(contig_analyzer_dirpath, 'minimap_output')
icarus_dirpath = join(contig_analyzer_dirpath, 'contigs_reports')
aux_dirpath = join(contig_analyzer_dirpath, 'aux')

corrected_dirpath = join(config['output_dir'], config['CORRECTED_DIRPATH'])
corrected_reference = config['reference'].split('/')[-1]

if not isdir(minimap_dirpath):
    os.makedirs(minimap_dirpath)

if not isdir(icarus_dirpath):
    os.makedirs(icarus_dirpath)

if not isdir(aux_dirpath):
    os.makedirs(aux_dirpath)

rule all:
    input:
        join(config['output_dir'], "report.txt")

def get_input_fastas(wildcards):
    return config["samples"][wildcards.sample]


rule correct_reference:
    input:
        fasta=config['reference'],
    log:
        out=join(corrected_dirpath, 'quast.log'),
        err=join(corrected_dirpath, 'quast.err')
    output:
        join(corrected_dirpath, corrected_reference)
    shell:
        "python -m scripts.preprocessing.correct_ref {input.fasta} {output} >{log.out} 2>{log.err}"

rule reference_stats:
    input:
        join(corrected_dirpath, corrected_reference)
    output:
        join(corrected_dirpath, corrected_reference + ".csv")
    shell:
        "python -m scripts.preprocessing.save_ref_stats {input} {output}"

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

rule contig_aligner:
    input:
        contig=join(corrected_dirpath,"{sample}.fasta"),
        reference=join(corrected_dirpath, corrected_reference),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv")
    conda:
        "envs/basic.yaml"
    log:
        out = join(contig_analyzer_dirpath, '{sample}_stdout.log'),
        err = join(contig_analyzer_dirpath, '{sample}_stderr.err')
    params:
        label="{sample}",
        output_dir=contig_analyzer_dirpath
    output:
        join(contig_analyzer_dirpath, 'contigs_report_{sample}.stdout')
    shell:
        #TODO: fix threads
        "python -m scripts.alignment.contig_aligner {input.reference} {params.label} {input.contig} {params.output_dir} "
        "{config[is_prokariote]} {input.reference_csv} {config[threads]} >{log.out} 2>{log.err}"

rule prepare_genome_analyzer:
    input:
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
        operons_fpaths=expand(config['operons']),
        features_fpaths=expand("{feature}", feature=config['features']),
    log:
        out = join(genome_analyzer_dirpath, 'quast.log'),
        err = join(genome_analyzer_dirpath, 'quast.err')
    params:
        labels=expand("{sample}", sample=config['samples']),
        output_dir=genome_analyzer_dirpath,
        features=expand("{feature}", feature=config['features_type']),
    output:
        join(genome_analyzer_dirpath, 'genome_info.txt'),
        join(genome_analyzer_dirpath, 'gene.csv'),
        join(genome_analyzer_dirpath, 'operon.csv')
    shell:
        "python -m scripts.gene_finding.prepare_genome_analyzer "
        "--output_dir {params.output_dir} --reference {input.reference_csv} --labels {params.labels} "
        "--features {params.features} --features_fpaths {input.features_fpaths} --operons {input.operons_fpaths} "
        ">{log.out} 2>{log.err}"

rule genome_analyzer:
    input:
        contig=join(corrected_dirpath,"{sample}.fasta"),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
        gene_container=join(genome_analyzer_dirpath, 'gene.csv'),
        operon_container=join(genome_analyzer_dirpath, 'operon.csv'),
    log:
        out = join(genome_analyzer_dirpath, '{sample}.log'),
        err = join(genome_analyzer_dirpath, '{sample}.err')
    params:
        label="{sample}",
        output_dir=genome_analyzer_dirpath,
        coords_dir=minimap_dirpath
    output:
        join(genome_analyzer_dirpath, '{sample}_info.txt')
    shell:
        "python -m scripts.gene_finding.genome_analyzer {params.output_dir} {input.reference_csv} {input.contig} {params.label} "
        "{params.coords_dir} {input.gene_container} {input.operon_container} >{log.out} 2>{log.err}"

rule save_stats:
    input:
        contigs=expand(join(corrected_dirpath, "{sample}.fasta"), sample=config['samples']),
        reference=join(corrected_dirpath, corrected_reference),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
        contig_stdout=expand(join(config['output_dir'], "contig_analyzer/contigs_report_{sample}.stdout"), sample=config['samples']),
        genome_info=expand(join(genome_analyzer_dirpath, '{sample}_info.txt'), sample=config['samples'])
    conda:
        "envs/basic.yaml"
    log:
        out=join(config['output_dir'], 'quast.log'),
        err=join(config['output_dir'], 'quast.err')
    params:
        contig_analyzer_dirpath=contig_analyzer_dirpath,
        genome_analyzer_dirpath=genome_analyzer_dirpath
    output:
        join(config['output_dir'], "report.txt")
    shell:
        "python -m scripts.make_reports {config[min_contig]} {input.reference_csv} {input.reference} {config[output_dir]} "
        "{params.contig_analyzer_dirpath} {params.genome_analyzer_dirpath} {input.contigs}  >{log.out} 2>{log.err}"

