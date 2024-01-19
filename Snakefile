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
print("\nOUTPUT PATH:")
print(OUTPUT)

# load in fastq path
input_path = os.path.abspath(config['inputs'])
input_df = pd.read_csv(input_path, comment="#")
samples = input_df['sample_id'].to_list()

# get input names 
input_names = utils.get_input_names(input_df, OUTPUT)

print("\nINPUT FILES:")
for x in input_names:
    print(x)


################ RULE FILES ################
include: "rules/reference.smk"
include: "rules/process_reads.smk"

rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/annotations.gtf',
        OUTPUT + 'references/geneTable.csv',
        OUTPUT + 'reports/seqkit_stats/raw_report.txt',
        OUTPUT + 'reports/seqkit_stats/demultiplexed_report.txt',
        expand(f"{OUTPUT}fastq/{{sid}}.raw.fastq.gz", sid=samples),
        expand(f"{OUTPUT}demultiplex/{{sid}}.done", sid=samples),
        expand(f"{OUTPUT}reports/fastqc/{{sid}}.report.html", sid=samples),
        expand(f"{OUTPUT}mapping/{{sid}}.bam.bai", sid=samples),
        expand(f"{OUTPUT}mapping/{{sid}}.tagged.bam", sid=samples),
        expand(f"{OUTPUT}reports/bamstats/{{sid}}.bamstats", sid=samples),
        OUTPUT + 'merged/merged.bam.bai',
        OUTPUT + 'merged/merged.stats',
        OUTPUT + 'merged/merged.bamstats',
        OUTPUT + 'counts/counts.txt',
        OUTPUT + 'scanpy/anndata.h5ad',
        OUTPUT + 'scanpy/anndata.processed.h5ad',


rule htseq_count:
    input:
        bam=OUTPUT + 'merged/merged.bam',
        annotations=config['gtf_path'],
    output:
        OUTPUT + "counts/counts.txt"
    conda:
        "envs/htseq_count.yml"
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


rule process_anndata:
    input:
        anndata=OUTPUT + "scanpy/anndata.h5ad",
    output:
        OUTPUT + "scanpy/anndata.processed.h5ad",
    params:
        config=str(BASE_DIR) + "/config/config.yaml"
    shell:
        """python scripts/process_anndata.py {input.anndata} {params.config} {output}"""
