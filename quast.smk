import os
from os.path import join, isdir

from scripts.gene_finding.run_busco import get_lineage
from src.qutils import get_path_to_program, get_dir_for_download

# configfile: "config.yaml"

corrected_dirpath = join(config['output_dir'], config['CORRECTED_DIRPATH'])
corrected_reference = config['reference'].split('/')[-1]

contig_analyzer_dirpath = join(config['output_dir'], 'contig_analyzer')
minimap_dirpath = join(contig_analyzer_dirpath, 'minimap_output')
icarus_dirpath = join(contig_analyzer_dirpath, 'contigs_reports')
aux_dirpath = join(contig_analyzer_dirpath, 'aux')
jobs_threads = max(1, config['threads']//len(config['samples']))

aligned_stats_dirpath = join(config['output_dir'], 'aligned_stats')

genome_analyzer_dirpath = join(config['output_dir'], 'genome_analyzer')

gene_pred_dirpath = join(config['output_dir'], 'gene_prediction')
tmp_gene_pred_dirpath = join(gene_pred_dirpath, 'tmp')

gene_pred_output = list()
if config['gene_prediction']:
    gene_pred_output = expand(join(gene_pred_dirpath, "{sample}.gff"), sample=config['samples'])
else:
    tmp_gene_pred_dirpath = None

if config['features']:
    features_input = expand(join(genome_analyzer_dirpath, "{feature}.csv"),feature=config['features'])
else:
    features_input = list()
    config['features_files'] = list()

busco_dirpath = join(config['output_dir'], 'busco')
lineage = get_lineage(is_prokaryote=config['is_prokaryote'],is_fungus=config['is_fungus'])
if config['busco'] and get_path_to_program('busco'):
    busco_output = expand(join(busco_dirpath, "{sample}/short_summary.specific." + lineage + ".{sample}.txt"), sample=config['samples'])
else:
    busco_output = list()

kmer_analyzer_dirpath = join(config['output_dir'], 'kmer_stats')
tmp_kmer_analyzer_dirpath = join(kmer_analyzer_dirpath, 'tmp')
if config['kmer_analysis'] and get_path_to_program('kmc'):
    kmer_output = expand(join(kmer_analyzer_dirpath, "{sample}.stat"), sample=config['samples'])
    os.makedirs(tmp_kmer_analyzer_dirpath)
else:
    kmer_output = list()

reads_analyzer_dirpath = join(config['output_dir'], 'reads_analyzer')
if config['reads_files']:
    reads_analyzer_output = expand(join(reads_analyzer_dirpath, "{sample}.stat"), sample=config['samples'])
    reads_analyzer_ref_output = join(reads_analyzer_dirpath, "reference.cov")
else:
    reads_analyzer_output = list()
    reads_analyzer_ref_output = list()

for d in [minimap_dirpath, icarus_dirpath, aux_dirpath, tmp_gene_pred_dirpath]:
    if d and not isdir(d):
        os.makedirs(d)

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
        out=join(corrected_dirpath, "{sample}.log"),
        err=join(corrected_dirpath, "{sample}.err")
    output:
        join(corrected_dirpath, "{sample}.fasta")
    shell:
        "python -m scripts.preprocessing.correct_contig {input} {output} {config[min_contig]} >{log.out} 2>{log.err}"

rule contig_aligner:
    input:
        contig=join(corrected_dirpath, "{sample}.fasta"),
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
        "python -m scripts.alignment.contig_aligner {input.reference} {params.label} {input.contig} {params.output_dir} "
        "{config[is_prokaryote]} {input.reference_csv} {jobs_threads} >{log.out} 2>{log.err}"

rule prepare_genome_analyzer:
    input:
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
        features_fpaths=config['features_files'],
    log:
        out = join(genome_analyzer_dirpath, 'quast.log'),
        err = join(genome_analyzer_dirpath, 'quast.err')
    params:
        labels=config['samples'],
        output_dir=genome_analyzer_dirpath,
        features=config['features'],
        features_option='--features' if features_input else '',
        features_paths_option='--features_fpaths' if features_input else '',
    output:
        info=join(genome_analyzer_dirpath, 'genome_info.txt'),
        features=features_input
    shell:
        "python -m scripts.gene_finding.prepare_genome_analyzer "
        "--output_dir {params.output_dir} --reference {input.reference_csv} --labels {params.labels} "
        "{params.features_option} {params.features} {params.features_paths_option} {input.features_fpaths}"
        ">{log.out} 2>{log.err}"

rule genome_analyzer:
    input:
        contig=join(corrected_dirpath, "{sample}.fasta"),
        contig_stdout=expand(join(config['output_dir'], "contig_analyzer/contigs_report_{sample}.stdout"), sample=config['samples']),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
    log:
        out = join(genome_analyzer_dirpath, '{sample}.log'),
        err = join(genome_analyzer_dirpath, '{sample}.err')
    params:
        label="{sample}",
        features_input=features_input,
        output_dir=genome_analyzer_dirpath,
        coords_dir=minimap_dirpath
    output:
        join(genome_analyzer_dirpath, '{sample}_info.txt')
    shell:
        "python -m scripts.gene_finding.genome_analyzer {params.output_dir} {input.reference_csv} {input.contig} {params.label} "
        "{params.coords_dir} {params.features_input} >{log.out} 2>{log.err}"

rule gene_prediction:
    input:
        contig=join(corrected_dirpath, "{sample}.fasta"),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
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

db_dirpath = join(get_dir_for_download('busco', 'busco_db', [lineage]), lineage)
rule download_busco:
    conda:
        "envs/busco.yaml"
    log:
        out=join(busco_dirpath, 'busco.log'),
        err=join(busco_dirpath, 'busco.err')
    params:
        output_dir=busco_dirpath
    output:
        directory(db_dirpath)
    shell:
        "python -m scripts.gene_finding.busco_download >{log.out} 2>{log.err}"

rule run_busco:
    input:
        database=db_dirpath,
        contig=join(corrected_dirpath, "{sample}.fasta"),
    conda:
        "envs/busco.yaml"
    log:
        out=join(busco_dirpath, '{sample}.log'),
        err=join(busco_dirpath, '{sample}.err')
    params:
        label="{sample}",
        output_dir=busco_dirpath
    output:
        join(busco_dirpath, "{sample}/short_summary.specific." + lineage + ".{sample}.txt")
    shell:
        "python -m scripts.gene_finding.run_busco {params.label} {input.contig} {params.output_dir} "
        "{input.database} {jobs_threads} >{log.out} 2>{log.err}"

rule kmer_analysis_ref:
    input:
        reference=join(corrected_dirpath, corrected_reference),
    conda:
        "envs/basic.yaml"
    log:
        out=join(kmer_analyzer_dirpath, 'kmc.log'),
        err=join(kmer_analyzer_dirpath, 'kmc.err')
    params:
        output_dir=kmer_analyzer_dirpath,
        tmp_output_dir=tmp_kmer_analyzer_dirpath
    output:
        join(tmp_kmer_analyzer_dirpath, 'kmc.downsampled.txt')
    shell:
        "python -m scripts.large_assembly_analysis.run_kmc_ref {params.output_dir} "
        "{params.tmp_output_dir} {input.reference} {config[threads]} >{log.out} 2>{log.err}"

rule kmer_analysis:
    input:
        contig=join(corrected_dirpath, "{sample}.fasta"),
        reference_csv=join(corrected_dirpath,corrected_reference + ".csv"),
        downsampled_kmers_fpath=join(tmp_kmer_analyzer_dirpath, 'kmc.downsampled.txt')
    conda:
        "envs/basic.yaml"
    log:
        out=join(kmer_analyzer_dirpath, '{sample}.log'),
        err=join(kmer_analyzer_dirpath, '{sample}.err')
    params:
        label="{sample}",
        ref_kmc_out_path=join(tmp_kmer_analyzer_dirpath, 'reference.kmc'),
        output_dir=kmer_analyzer_dirpath,
        tmp_output_dir=tmp_kmer_analyzer_dirpath
    output:
        join(kmer_analyzer_dirpath, '{sample}.stat')
    shell:
        "python -m scripts.large_assembly_analysis.run_kmc {params.output_dir} {params.tmp_output_dir} "
        "{params.ref_kmc_out_path} {input.reference_csv} {config[is_prokaryote]} {input.downsampled_kmers_fpath} "
        "{input.contig} {params.label} "
        "{jobs_threads} >{log.out} 2>{log.err}"

rule align_reads_contigs:
    input:
        contigs=join(corrected_dirpath, "{sample}.fasta"),
    conda:
        "envs/basic.yaml"
    log:
        out=join(reads_analyzer_dirpath, '{sample}.log'),
        err=join(reads_analyzer_dirpath, '{sample}.err')
    params:
        reads_analyzer_dirpath=reads_analyzer_dirpath,
        reads_fpaths=config['reads_files'],
        reads_types=config['reads_types'],
        reads_option='--reads_fpaths',
        reads_types_option='--reads_types',
    output:
        join(reads_analyzer_dirpath, '{sample}.bam')
    shell:
        "python -m scripts.read_mapping.align_reads "
        "--reference {input.contigs} --output_dir {reads_analyzer_dirpath} --threads {jobs_threads} "
        "{params.reads_types_option} {params.reads_types} {params.reads_option} {params.reads_fpaths} >{log.out} 2>{log.err}"

rule align_reads_ref:
    input:
        reference=join(corrected_dirpath, corrected_reference),
    conda:
        "envs/basic.yaml"
    log:
        out=join(reads_analyzer_dirpath, 'reference.log'),
        err=join(reads_analyzer_dirpath, 'reference.err')
    params:
        reads_analyzer_dirpath=reads_analyzer_dirpath,
        reads_types=config['reads_types'],
        reads_fpaths=config['reads_files'],
        reads_option='--reads_fpaths',
        reads_types_option='--reads_types',
    output:
        join(reads_analyzer_dirpath, 'reference.bam')
    shell:
        "python -m scripts.read_mapping.align_reads "
        "--reference {input.reference} --output_dir {reads_analyzer_dirpath} --threads {config[threads]} "
        "{params.reads_types_option} {params.reads_types} {params.reads_option} {params.reads_fpaths} >{log.out} 2>{log.err}"

rule calculate_coverage:
    input:
        contigs=join(corrected_dirpath, "{sample}.fasta"),
        bam_file=join(reads_analyzer_dirpath, "{sample}.bam"),
    conda:
        "envs/basic.yaml"
    log:
        out=join(reads_analyzer_dirpath, '{sample}.log'),
        err=join(reads_analyzer_dirpath, '{sample}.err')
    params:
        reads_analyzer_dirpath=reads_analyzer_dirpath,
        reads_types=config['reads_types'],
        reads_option='--reads_fpaths',
        reads_types_option='--reads_types',
    output:
        join(reads_analyzer_dirpath, '{sample}.stat')
    shell:
        "python -m scripts.read_mapping.analyze_reads {params.reads_analyzer_dirpath} {input.contigs} "
        "{input.bam_file} {jobs_threads} >{log.out} 2>{log.err}"

rule analyze_ref:
    input:
        reference=join(corrected_dirpath, corrected_reference),
        bam_file=join(reads_analyzer_dirpath, "reference.bam"),
    conda:
        "envs/basic.yaml"
    log:
        out=join(reads_analyzer_dirpath, 'reference.log'),
        err=join(reads_analyzer_dirpath, 'reference.err')
    params:
        reads_analyzer_dirpath=reads_analyzer_dirpath,
        reads_types=config['reads_types'],
        reads_option='--reads_fpaths',
        reads_types_option='--reads_types',
    output:
        reads_analyzer_ref_output
    shell:
        "python -m scripts.read_mapping.process_reference "
        "{input.reference} {input.bam_file} {reads_analyzer_dirpath} {config[threads]} {config[search_sv]} "
        ">{log.out} 2>{log.err}"

rule save_stats:
    input:
        contigs=expand(join(corrected_dirpath, "{sample}.fasta"), sample=config['samples']),
        reference=join(corrected_dirpath, corrected_reference),
        reference_csv=join(corrected_dirpath, corrected_reference + ".csv"),
        contig_stdout=expand(join(config['output_dir'], "contig_analyzer/contigs_report_{sample}.stdout"), sample=config['samples']),
        genome_info=expand(join(genome_analyzer_dirpath, '{sample}_info.txt'), sample=config['samples']),
        features=features_input,
        gene_pred_output=gene_pred_output,
        reads_analyzer_output=reads_analyzer_output,
        reads_analyzer_ref_output=reads_analyzer_ref_output,
        busco_output=busco_output,
        kmer_output=kmer_output
    conda:
        "envs/basic.yaml"
    log:
        out=join(config['output_dir'], 'quast.log'),
        err=join(config['output_dir'], 'quast.err')
    params:
        contig_analyzer_dirpath=contig_analyzer_dirpath,
        genome_analyzer_dirpath=genome_analyzer_dirpath,
        features_option='--features' if features_input else '',
        tmp_gene_pred_dirpath=tmp_gene_pred_dirpath if gene_pred_output else None,
        kmer_analyzer_dirpath=kmer_analyzer_dirpath if kmer_output else None,
        reads_analyzer_dirpath=reads_analyzer_dirpath if reads_analyzer_output else None,
        busco_dirpath=busco_dirpath if busco_output else None,
        lineage=lineage,
    output:
        join(config['output_dir'], "report.txt")
    shell:
        "python -m scripts.make_reports -m {config[min_contig]} --csv {input.reference_csv} -r {input.reference} "
        "-o {config[output_dir]} --contig_analyzer_dirpath {params.contig_analyzer_dirpath} "
        "{params.features_option} {input.features} --genome_analyzer_dirpath {params.genome_analyzer_dirpath} "
        "--gene_pred_dirpath {params.tmp_gene_pred_dirpath} "
        "--kmer_analyzer_dirpath {params.kmer_analyzer_dirpath} "
        "--reads_analyzer_dirpath {params.reads_analyzer_dirpath} "
        "--busco_dirpath {params.busco_dirpath} --lineage {params.lineage} --contigs_fpaths {input.contigs} >{log.out} 2>{log.err}"

