#!/usr/bin/env python3

import peppy
import pathlib2

#############
# FUNCTIONS #
#############

def get_peptide_dbs(wildcards):
    input_keys = ['peptide_db']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

##needed to get BUSCO running in new folder
def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_species = pep.sample_table['sample_name']

bbduk_adapters = '/adapters.fa'
bbduk_ref = '/phix174_ill.ref.fa.gz'

##############
# CONTAINERS #
##############

busco_container = 'docker://ezlabgva/busco:v4.0.2_cv1'
tidyverse_container = 'docker://rocker/tidyverse'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

#########
# RULES #
#########

rule target:
    input:
       ##GC content plots
        expand('output/gc_depth/{hic_species}/hic_gc_boxplot.pdf', hic_species=['Mh', 'IR']),
        expand('output/gc_depth/{not_hic_species}/busco_gc_boxplot.pdf', not_hic_species=['MO', 'FR']),
        expand('output/gc_stats/{species}/kruskal_res.txt', species=['MO', 'FR', 'Mh']),
        expand('output/depth_stats/{species}/kruskal_res.txt', species=['MO', 'FR', 'Mh']),
        expand('output/depth_analysis/{species}_boxplot_y_zoom.pdf', species=['MO', 'FR', 'Mh'])

#######################
## GC vs depth table ##
#######################

rule depth_stat_test:
    input:
        gc_depth_table = 'output/gc_depth/{species}/gc_vs_depth_table.csv'
    output:
        shapiro_res = 'output/depth_stats/{species}/gc_shapiro_res.txt',
        summary_stats = 'output/depth_stats/{species}/summary_stats.txt',
        kruskal_res = 'output/depth_stats/{species}/kruskal_res.txt',
        wilcox_res = 'output/depth_stats/{species}/wilcox_res.txt'
    log:
        'output/logs/gc_depth_stats/{species}_depth_stats.log'
    script:
        'src/depth_stat_test.R'

rule GC_stat_test:
    input:
        gc_depth_table = 'output/gc_depth/{species}/gc_vs_depth_table.csv'
    output:
        shapiro_res = 'output/gc_stats/{species}/gc_shapiro_res.txt',
        summary_stats = 'output/gc_stats/{species}/summary_stats.txt',
        kruskal_res = 'output/gc_stats/{species}/kruskal_res.txt',
        wilcox_res = 'output/gc_stats/{species}/wilcox_res.txt'
    log:
        'output/logs/gc_depth_stats/{species}_gc_stats.log'
    script:
        'src/GC_stat_test.R'

rule plot_depth:
    input:
        gc_table = 'output/gc_depth/{species}/gc_table.csv',
        coverage_file = 'output/samtools_coverage/{species}/coverage.out'
    output:
        depth_boxplot = 'output/gc_depth/{species}/depth_boxplot.pdf',
        gc_depth_table = 'output/gc_depth/{species}/gc_vs_depth_table.csv'
    log:
        'output/logs/gc_depth/{species}_plot_depth.log'
    script:
        'src/plot_depth.R'

#################################
## sequencing depth of contigs ##
#################################

rule samtools_coverage:
    input:
        bam = 'output/samtools/{species}/sorted.bam'
    output:
        coverage_out = 'output/samtools_coverage/{species}/coverage.out'
    log:
        'output/logs/samtools_coverage_{species}.log'
    threads:
        20
    shell:
        'samtools coverage '
        '{input.bam} '
        '-o {output.coverage_out} '
        '2> {log}'

##don't work because MO and FR have so many contigs (3000+), for Mh shows much the same as scatterplot
rule depth_boxplot:
    input:
        samtools_depth = 'output/samtools_depth/{species}/depth.out',
        gc_depth_table = 'output/gc_depth/{species}/gc_vs_depth_table.csv'
    output:
        boxplot = 'output/depth_analysis/{species}_boxplot.jpeg',
        boxplot_y_zoom = 'output/depth_analysis/{species}_boxplot_y_zoom.pdf'
    threads:
        20
    log:
        'output/logs/boxplots/{species}_depth_boxplot.log'
    script:
        'src/depth_boxplot.R'

##include -a option - to print all positions even if depth = 0
rule samtools_depth:
    input:
        sorted_bam = 'output/samtools/{species}/sorted.bam'
    output:
        depth_out = 'output/samtools_depth/{species}/depth.out'
    log:
        'output/logs/samtools_depth_{species}.log'
    threads:
        20
    shell:
        'samtools depth '
        '{input.sorted_bam} '
        '> {output.depth_out} '
        '-a '
        '2> {log}'

rule samtools_index:
    input:
        bam = 'output/samtools/{species}/sorted.bam'
    output:
        index = 'output/samtools/{species}/sorted.bam.bai'
    log:
        'output/logs/samtools_index_{species}.log'
    threads:
        20
    shell:
        'samtools index '
        '{input.bam} '
        '2> {log}'

rule samtools_sort:
    input:
        sam = 'output/bwa/{species}/bwa_mem.sam'
    output:
        sorted_bam = 'output/samtools/{species}/sorted.bam'
    log:
        'output/logs/samtools_sort_{species}.log'
    threads:
        20
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '2> {log}'

##bwa-mem to map for read depth of scaffolds
rule bwa_mem:
    input:
        index = 'output/bwa/{species}/index.bwt',
        trimr1 = 'output/bbduk_trim_dna/{species}/{species}_trimr1.fq.gz',
        trimr2 = 'output/bbduk_trim_dna/{species}/{species}_trimr2.fq.gz'
    output:
        sam = 'output/bwa/{species}/bwa_mem.sam'
    params:
        index_dir = 'output/bwa/{species}/index'
    threads:
        50
    log:
        'output/logs/bwa_mem_{species}.log'
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.index_dir} '
        '{input.trimr1} '
        '{input.trimr2} '
        '> {output.sam} '
        '2> {log}'

rule bwa_index:
    input:
        genome = 'data/final_genomes/{species}.fa'
    output:
        index = 'output/bwa/{species}/index.bwt'
    params:
        outdir = 'output/bwa/{species}/index'
    threads:
        20
    log:
        'output/logs/bwa_index_{species}.log'
    shell:
        'bwa index '
        '{input.genome} '
        '-p {params.outdir} '
        '2> {log} '

##trim and decontaminate DNA reads to map onto genome
rule bbduk_trim_dna:
    input:
        filr1 = 'output/bbduk_trim_dna/{species}/{species}_filr1.fq.gz',
        filr2 = 'output/bbduk_trim_dna/{species}/{species}_filr2.fq.gz'
    output:
        trimr1 = 'output/bbduk_trim_dna/{species}/{species}_trimr1.fq.gz',
        trimr2 = 'output/bbduk_trim_dna/{species}/{species}_trimr2.fq.gz',
        t_stats = 'output/bbduk_trim_dna/{species}/trim-stats.txt'
    log:
        trim = 'output/logs/bbduk_trim_dna/{species}_trim.log'
    params:
        adapters = bbduk_adapters
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in1={input.filr1} '
        'in2={input.filr2} '
        'int=f '
        'out1={output.trimr1} '
        'out2={output.trimr2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '

rule bbduk_filter_dna:  
    input:
        r1 = 'data/dna_reads/{species}/{species}_1.fastq.gz',
        r2 = 'data/dna_reads/{species}/{species}_2.fastq.gz'
    output:
        filr1 = 'output/bbduk_trim_dna/{species}/{species}_filr1.fq.gz',
        filr2 = 'output/bbduk_trim_dna/{species}/{species}_filr2.fq.gz',
        f_stats = 'output/bbduk_trim_dna/{species}/filter-stats.txt'
    log:
        filter = 'output/logs/bbduk_trim_dna/{species}_filter.log'
    params:
        ref = bbduk_ref
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in1={input.r1} '
        'in2={input.r2} '
        'out1={output.filr1} '
        'out2={output.filr2} '
        'ref={params.ref} '
        'hdist=1 '
        'stats={output.f_stats} '       
        '2> {log.filter} ' 

###########################
## GC content of contigs ##
###########################

##only run busco if Hi-C scaffold list not available
ruleorder: plot_GC_content_hic > plot_GC_content_not_hic

rule plot_GC_content_not_hic:
    input:
        gc = 'output/bb_stats/{not_hic_species}/gc.txt',
        gc_hist_file = 'output/bb_stats/{not_hic_species}/gc_hist.out',
        viral_contig_list = 'data/viral_contig_lists/{not_hic_species}.csv',
        busco_results_file = 'output/busco/{not_hic_species}/run_hymenoptera_odb10/full_table.tsv'
    output:
        busco_only = 'output/gc_depth/{not_hic_species}/busco_gc_boxplot.pdf',
        gc_histogram = 'output/gc_depth/{not_hic_species}/gc_histogram.pdf',
        gc_table = 'output/gc_depth/{not_hic_species}/gc_table.csv'
    log:
        'output/logs/gc_depth/{not_hic_species}_plot_gc.log'
    script:
        'src/plot_GC_content_not_hic.R'

rule plot_GC_content_hic:
    input:
        gc = 'output/bb_stats/{hic_species}/gc.txt',
        gc_hist = 'output/bb_stats/{hic_species}/gc_hist.out',
        viral_contig_list = 'data/viral_contig_lists/{hic_species}.csv',
        hic_scaffold_list = 'data/hi-c_genomes/{hic_species}_hic_scaffold_ids.txt'
    output:
        hic_only = 'output/gc_depth/{hic_species}/hic_gc_boxplot.pdf',
        gc_histogram = 'output/gc_depth/{hic_species}/gc_histogram.pdf',
        gc_table = 'output/gc_depth/{hic_species}/gc_table.csv'
    log:
        'output/logs/gc_depth/{hic_species}_plot_gc.log'
    script:
        'src/plot_GC_content_hic.R'

rule bb_stats:
    input:
        genome = 'data/final_genomes/{species}.fa'
    output:
        gc = 'output/bb_stats/{species}/gc.txt',
        stats = 'output/bb_stats/{species}/bb_stats.out',
        gc_hist = 'output/bb_stats/{species}/gc_hist.out'
    log:
        'output/logs/bbstats/{species}.log'
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input.genome} '
        'out={output.stats} '
        'gc={output.gc} '
        'gcformat=4 '
        'gchist={output.gc_hist} '

########################################
## BUSCO for ID of eukaryotic contigs ##
########################################

##for two assemblies without hi-c need busco to compare busco containing to viral
rule busco_non_hic_genomes:
    input:
        genome = 'data/final_genomes/{not_hic_species}.fa'
    output:
        'output/busco/{not_hic_species}/run_hymenoptera_odb10/full_table.tsv'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                            'busco_{not_hic_species}.log'))
    params:
        wd = 'output/busco',
        genome = lambda wildcards, input: resolve_path(input.genome)
    singularity:
        busco_container
    threads:
        20
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--force '
        '--in {params.genome} '
        '--out {wildcards.not_hic_species} '
        '--lineage hymenoptera_odb10 '
        '--augustus_species nasonia '
        '--cpu {threads} '
        '--mode genome '
        '-f '
        '&> {log} '
