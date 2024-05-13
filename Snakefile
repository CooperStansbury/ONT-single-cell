from datetime import datetime
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
[print(x) for x in samples]

# timestamp
now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
test_file = OUTPUT + "test/" + now + ".txt"


################ RULE FILES ################
include: "rules/reference.smk"
include: "rules/demultiplex.smk"
include: "rules/core.smk"
include: "rules/anndata.smk"
include: "rules/v5tags.smk"
include: "rules/isoquant.smk"



rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/transcripts.fa',
        OUTPUT + 'references/annotations.gtf',
        OUTPUT + 'references/geneTable.csv',
        OUTPUT + 'reports/seqkit_stats/raw_report.txt',
        OUTPUT + 'reports/seqkit_stats/demultiplexed_report.txt',
        expand(OUTPUT + "fastq/{sid}.raw.fastq.gz", sid=samples),
        expand(OUTPUT + "demultiplex/{sid}.done", sid=samples),
        expand(OUTPUT + "reports/fastqc/{sid}.report.html", sid=samples),
        expand(OUTPUT + "mapping/{sid}.bam.bai", sid=samples),
        expand(OUTPUT + "mapping/{sid}.tagged.bam", sid=samples),
        expand(OUTPUT + "reports/bamstats/{sid}.bamstats", sid=samples),
        expand(OUTPUT + "individual_counts/{sid}.counts.txt", sid=samples),
        OUTPUT + 'merged/merged.bam.bai',
        OUTPUT + 'merged/merged.stats',
        OUTPUT + 'merged/merged.bamstats',
        OUTPUT + 'merged/merged.counts.txt',
        OUTPUT + 'scanpy/raw.anndata.h5ad',
        OUTPUT + "scanpy/processed.anndata.h5ad",
        OUTPUT + "scanpy/clustered.anndata.h5ad",
        OUTPUT + 'v5_tagged/v5_tagged.fastq.gz',
        OUTPUT + "isoquant/annotations.db",
        test_file,

rule test:
    output:
        touch(test_file),
    





