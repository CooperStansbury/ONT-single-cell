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
include: "rules/demultiplex.smk"
include: "rules/flair.smk"
include: "rules/core.smk"
include: "rules/anndata.smk"
include: "rules/v5tags.smk"
include: "rules/nanocount.smk"



rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/transcripts.fa',
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
        expand(f"{OUTPUT}counts/individual/{{sid}}.counts.txt", sid=samples),
        expand(f"{OUTPUT}nanocount/mapping/{{sid}}.bam", sid=samples),
        expand(f"{OUTPUT}nanocount/{{sid}}.tx_counts.tsv", sid=samples),
        OUTPUT + 'merged/merged.bam.bai',
        OUTPUT + 'merged/merged.stats',
        OUTPUT + 'merged/merged.bamstats',
        OUTPUT + 'counts/counts.txt',
        OUTPUT + 'scanpy/anndata.h5ad',
        OUTPUT + 'scanpy/anndata.processed.h5ad',
        OUTPUT + 'nanocount/merged/merged_tx_counts.tsv',
        OUTPUT + 'v5_tagged/v5_tagged.fastq.gz',



rule archive:
    input:
        expand(f"{OUTPUT}nanocount/mapping/{{sid}}.bam", sid=samples),
        expand(f"{OUTPUT}nanocount/{{sid}}.tx_counts.tsv", sid=samples),
        expand(f"{OUTPUT}flair/align/{{sid}}.bam", sid=samples),
        expand(f"{OUTPUT}flair/correct/{{sid}}_all_corrected.bed", sid=samples),
        expand(f"{OUTPUT}flair/collapse/{{sid}}.isoforms.bed", sid=samples),
        expand(f"{OUTPUT}velocyto/{{sid}}.done", sid=samples),





