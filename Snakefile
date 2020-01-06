#!/usr/bin/env python3

import multiprocessing
import pandas


def find_input_reads(wildcards):
    return {
        'r1': sorted(set(sample_data['r1_path'].tolist())),
        'r2': sorted(set(sample_data['r2_path'].tolist()))
    }


def get_min_cutoff(wildcards):
    cutoff_file = checkpoints.genotype.get(**wildcards).output['cutoffs']
    cutoffs = pandas.read_csv(cutoff_file,
                              header=None,
                              index_col=0)
    return cutoffs.loc['min_depth', 1]


def get_max_cutoff(wildcards):
    cutoff_file = checkpoints.genotype.get(**wildcards).output['cutoffs']
    cutoffs = pandas.read_csv(cutoff_file,
                              header=None,
                              index_col=0)
    return cutoffs.loc['max_depth', 1]


###########
# GLOBALS #
###########

ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
sample_csv = 'data/sample_data.csv'

# software
honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.4'
    '@d9ac3b038cc4244bb2d5a0863a18168570fab218')
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
        'output/010_genotypes/calls.vcf.gz',
        'output/020_filtered-genotypes/filtered.vcf'


rule filter:
    input:
        cutoffs = 'output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/calls.vcf.gz'
    output:
        'output/020_filtered-genotypes/filtered.vcf'
    params:
        min_depth = get_min_cutoff,
        max_depth = get_max_cutoff,
        maf = 0.1,
        max_missing = 0.9,
        qual = 30
    log:
        'output/logs/filter.log'
    singularity:
        honeybee_genotype_pipeline
    shell:
        'vcftools '
        '--gzvcf {input.vcf} '
        '--maf {params.maf} '
        '--max-missing {params.max_missing} '
        '--minQ {params.qual} '
        '--min-meanDP {params.min_depth} '
        '--max-meanDP {params.max_depth} '
        '--minDP {params.min_depth} '
        '--maxDP {params.max_depth} '
        '--recode '
        '--stdout '
        '> {output} '
        '2> {log}'


checkpoint genotype:
    input:
        unpack(find_input_reads), # just to make sure combine gets run first
        csv = sample_csv,
        ref = ref
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

