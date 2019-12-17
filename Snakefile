#!/usr/bin/env python3

import multiprocessing
import pandas


def find_input_reads(wildcards):
    return {
        'r1': sorted(set(sample_data['r1_path'].tolist())),
        'r2': sorted(set(sample_data['r2_path'].tolist()))
    }


###########
# GLOBALS #
###########

ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
sample_csv = 'data/sample_data.csv'

# software
honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.3'
    '@e7e37748bde42ab8d6ad8dffecd5ca008089276c')
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'

########
# MAIN #
########

# get a list of samples
sample_data = pandas.read_csv(sample_csv,
                              index_col='sample')
all_samples = sorted(set(sample_data.index))

#########
# RULES #
#########

rule target:
    input:
        'output/010_genotypes/calls.vcf.gz'

rule genotype:
    input:
        find_input_reads,        # just to make sure combine gets run first
        csv = sample_csv,
        ref = honeybee_ref
    output:
        cutoffs = 'output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/calls.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        ref = 'output/010_genotypes/015_ref/ref.fasta',
        fai = 'output/010_genotypes/015_ref/ref.fasta.fai',
    params:
        wd = 'output/010_genotypes',
        ploidy = '2'
    log:
        'output/logs/genotype.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        honeybee_genotype_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--ploidy {params.ploidy} '
        '--threads {threads} '
        '--restart_times 1 '
        '&> {log}'


# combine reads for BB indivs
rule combine_reads:
    input:
        r1_1 = 'data/bb_lane1/{sample}_{run}/{sample}_{run}_R1.fq.gz',
        r1_2 = 'data/bb_lane1/{sample}_{run}/{sample}_{run}_R2.fq.gz',
        r2_1 = 'data/bb_lane2/{sample}_{run}/{sample}_{run}_R1.fq.gz',
        r2_2 = 'data/bb_lane2/{sample}_{run}/{sample}_{run}_R2.fq.gz'
    output:
        r1 = temp('output/000_tmp/reads/{sample}_{run}_R1.fq.gz'),
        r2 = temp('output/000_tmp/reads/{sample}_{run}_R2.fq.gz')
    singularity:
        samtools
    shell:
        'cat {input.r1_1} {input.r2_1} > {output.r1} & '
        'cat {input.r1_2} {input.r2_2} > {output.r2} & '
        'wait'

