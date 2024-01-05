import pandas as pd
import yaml
from pathlib import Path
import re
import os
import sys
from utils import utils

BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"

# big picture variables
OUTPUT = config['output_path']

rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/geneTable.csv',
        OUTPUT + 'fastqc/report.html',
        OUTPUT + 'alignments/alignments.bamstats',
        OUTPUT + 'alignments/alignments.sorted.bam.bai',
        OUTPUT + 'alignments/alignments.tagged.bam',
        OUTPUT + "counts/counts.txt",
        OUTPUT + "scanpy/anndata.h5ad",
     
    
rule getReference:
    input:
        refgenome=config['ref_path'],
    output:
        OUTPUT + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"


rule prepReference:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        ref=OUTPUT + 'references/reference.fa',
        flag=touch(OUTPUT + 'reference.done')
    shell:
        "cat {input} | gzip -d > {output.ref}"
        
            
rule get_gene_table:
    input:
        annotations=config['gtf_path'],
    output:
        OUTPUT + "references/geneTable.csv"
    shell:
        "python scripts/getGeneTable.py {input} {output}"


rule minimap2_index:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        OUTPUT + 'references/reference.mmi'
    shell:
        "minimap2 -d {output} {input.refgenome}"


rule fastqc:
    input:
        fastq=config['fastq'],
    output:
        html=OUTPUT + "fastqc/report.html",
        zip=OUTPUT + "fastqc/report_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/report.log"
    threads:
        8
    wrapper:
        "v1.29.0/bio/fastqc"


rule minimap2_align:
   input:
       fastq=config['fastq'],
       refgenome=OUTPUT + 'references/reference.fa.gz',
       refindex=OUTPUT + 'references/reference.mmi',
   output:        
       OUTPUT + 'alignments/alignments.sam'
   params:
       args=config['minimap2_args'],
       threads=36
   shell:
       "minimap2 {params.args} -t {params.threads} {input.refgenome} {input.fastq} > {output}"


rule make_bam:
    input:
        OUTPUT + 'alignments/alignments.sam'
    output:       
        OUTPUT + 'alignments/alignments.bam'
    shell:
        "samtools view -Sb {input} > {output}"


rule samtools_sort:
    input:
         OUTPUT + 'alignments/alignments.bam'
    output:
         OUTPUT + 'alignments/alignments.sorted.bam'
    shell:
        "samtools sort -T {input} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        OUTPUT + 'alignments/alignments.sorted.bam'
    output:
        OUTPUT + 'alignments/alignments.sorted.bam.bai'
    shell:
        "samtools index {input}"


rule bamtools_stats:
    input:
        OUTPUT + 'alignments/alignments.bam'
    output:
        OUTPUT + 'alignments/alignments.bamstats'
    params:
        "-insert" # optional summarize insert size data
    log:
        "logs/bamtools/stats/bamstats.log"
    wrapper:
        "v2.1.1/bio/bamtools/stats"

rule tag_bam:
    input:
        OUTPUT + 'alignments/alignments.sorted.bam'
    output:
        OUTPUT + 'alignments/alignments.tagged.bam'
    shell:
        "python scripts/tag_bam.py {input} {output}"


rule htseq_count:
    input:
        bam=OUTPUT + 'alignments/alignments.tagged.bam',
        annotations=config['gtf_path'],
    output:
        OUTPUT + "counts/counts.txt"
    params:
        d=int(config['umi_distance'])
    shell:
        "htseq-count-barcodes --nonunique all {input.bam} {input.annotations} > {output}"


rule make_andata:
    input:
        counts=OUTPUT + "counts/counts.txt",
        genes=OUTPUT + "references/geneTable.csv",
    output:
        OUTPUT + "scanpy/anndata.h5ad",
    shell:
        """python scripts/make_andata.py {input.counts} {input.genes} {output}"""
        
