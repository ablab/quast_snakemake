import os
from os.path import join, isdir

configfile: "config.yaml"

corrected_dirpath = join(config['output_dir'], config['CORRECTED_DIRPATH'])
corrected_reference = config['reference'].split('/')[-1]

contig_analyzer_dirpath = join(config['output_dir'], 'contig_analyzer')
minimap_dirpath = join(contig_analyzer_dirpath, 'minimap_output')
icarus_dirpath = join(contig_analyzer_dirpath, 'contigs_reports')
aux_dirpath = join(contig_analyzer_dirpath, 'aux')

genome_analyzer_dirpath = join(config['output_dir'], 'genome_analyzer')
glimmer_dirpath = join(config['output_dir'], 'gene_prediction')
tmp_glimmer_dirpath = join(glimmer_dirpath, 'tmp')

for dir in [minimap_dirpath, icarus_dirpath, aux_dirpath, tmp_glimmer_dirpath]:
    if not isdir(dir):
        os.makedirs(dir)

glimmer_output = list()
if config['gene_prediction']:
    glimmer_output = expand(join(glimmer_dirpath, "{sample}_glimmer.gff"), sample=config['samples'])

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
        features_fpaths=config['features'],
    log:
        out = join(genome_analyzer_dirpath, 'quast.log'),
        err = join(genome_analyzer_dirpath, 'quast.err')
    params:
        labels=config['samples'],
        features=config['features_type'],
        output_dir=genome_analyzer_dirpath,
    output:
        info=join(genome_analyzer_dirpath, 'genome_info.txt'),
        containers=expand(join(genome_analyzer_dirpath, "{feature}.csv"), feature=config['features_type']),
    shell:
        "python -m scripts.gene_finding.prepare_genome_analyzer "
        "--output_dir {params.output_dir} --reference {input.reference_csv} --labels {params.labels} "
        "--features {params.features} --features_fpaths {input.features_fpaths} "
        ">{log.out} 2>{log.err}"

rule genome_analyzer:
    input:
        contig=join(corrected_dirpath,"{sample}.fasta"),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
        containers=expand(join(genome_analyzer_dirpath, "{feature}.csv"), feature=config['features_type']),
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
        "{params.coords_dir} {input.containers} >{log.out} 2>{log.err}"

rule glimmer:
    input:
        contig=join(corrected_dirpath,"{sample}.fasta"),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
    conda:
        "envs/basic.yaml"
    log:
        out = join(glimmer_dirpath, '{sample}.log'),
        err = join(glimmer_dirpath, '{sample}.err')
    params:
        label="{sample}",
        tmp_dir=tmp_glimmer_dirpath,
        output_dir=glimmer_dirpath,
        coords_dir=minimap_dirpath
    output:
        join(glimmer_dirpath, '{sample}_glimmer.gff')
    shell:
        "python -m scripts.gene_finding.glimmer {params.output_dir} {input.contig} {params.label} "
        "{params.tmp_dir} >{log.out} 2>{log.err}"

rule save_stats:
    input:
        contigs=expand(join(corrected_dirpath, "{sample}.fasta"), sample=config['samples']),
        reference=join(corrected_dirpath, corrected_reference),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
        contig_stdout=expand(join(config['output_dir'], "contig_analyzer/contigs_report_{sample}.stdout"), sample=config['samples']),
        genome_info=expand(join(genome_analyzer_dirpath, '{sample}_info.txt'), sample=config['samples']),
        features=expand(join(genome_analyzer_dirpath, "{feature}.csv"), feature=config['features_type']),
        glimmer_output=glimmer_output
    conda:
        "envs/basic.yaml"
    log:
        out=join(config['output_dir'], 'quast.log'),
        err=join(config['output_dir'], 'quast.err')
    params:
        contig_analyzer_dirpath=contig_analyzer_dirpath,
        genome_analyzer_dirpath=genome_analyzer_dirpath,
        tmp_glimmer_dirpath=tmp_glimmer_dirpath
    output:
        join(config['output_dir'], "report.txt")
    shell:
        "python -m scripts.make_reports -m {config[min_contig]} --csv {input.reference_csv} -r {input.reference} "
        "-o {config[output_dir]} --contig_analyzer_dirpath {params.contig_analyzer_dirpath} "
        "--genome_analyzer_dirpath {params.genome_analyzer_dirpath} --glimmer_dirpath {params.tmp_glimmer_dirpath} "
        "--contigs_fpaths {input.contigs} --features {input.features} >{log.out} 2>{log.err}"

