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

weng_ids = [
    "young2_HSC", 
    "young1_all_t2", 
    "old2_BMMC_HSPC",
    "young2_all", 
    "old1_BMMC_HSPC",
    "young1_all_t1",
]


geneformer_ids = [
    "iHSC",
    "pellin",
    "tabula_sapiens",
    "young2_HSC", 
    "young1_all_t2", 
    "old2_BMMC_HSPC",
    "young2_all", 
    "old1_BMMC_HSPC",
    "young1_all_t1",
]


################ RULE FILES ################
include: "rules/reference.smk"
include: "rules/demultiplex.smk"
include: "rules/core.smk"
include: "rules/anndata.smk"
include: "rules/v5tags.smk"
include: "rules/geneformer.smk"
include: "rules/isoquant.smk"
include: "rules/rna_velocity.smk"



rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/transcripts.fa',
        OUTPUT + 'references/annotations.gtf',
        OUTPUT + 'references/geneTable.csv',
        expand(OUTPUT + "fastq/{sid}.raw.fastq.gz", sid=samples),
        expand(OUTPUT + "demultiplex/{sid}.done", sid=samples),
        expand(OUTPUT + "reports/fastqc/{sid}.report.html", sid=samples),
        expand(OUTPUT + "mapping/{sid}.bam.bai", sid=samples),
        expand(OUTPUT + "mapping/{sid}.tagged.bam.bai", sid=samples),
        expand(OUTPUT + "reports/bamstats/{sid}.bamstats", sid=samples),
        expand(OUTPUT + "individual_counts/{sid}.counts.txt", sid=samples),
        expand(OUTPUT + "v5_tagged/{sid}.tagged.csv", sid=samples),
        OUTPUT + 'v5_tagged/factors.mmi',
        OUTPUT + 'v5_tagged/all_reads_factor_mapped.bam',
        OUTPUT + 'reports/seqkit_stats/raw_report.txt',
        OUTPUT + 'reports/seqkit_stats/demultiplexed_report.txt',
        OUTPUT + "v5_tagged/read_ids.txt",
        OUTPUT + 'v5_tagged/v5_result.factor_table.csv',
        OUTPUT + 'merged/merged.bam.bai',
        OUTPUT + 'merged/merged.stats',
        OUTPUT + 'merged/merged.bamstats',
        OUTPUT + 'merged/merged.counts.txt',
        OUTPUT + 'scanpy/raw.anndata.h5ad',
        OUTPUT + "scanpy/processed.anndata.h5ad",
        OUTPUT + "scanpy/clustered.anndata.h5ad",
        OUTPUT + 'v5_tagged/v5_result.table.csv',
        OUTPUT + 'v5_tagged/all_reads_merged.fastq',
        expand(OUTPUT + "isoquant/{sid}.done", sid=samples),
        OUTPUT + "isoquant/annotations.db",
        OUTPUT + "isoquant_prepared/gene_counts.csv",
        OUTPUT + "isoquant_prepared/transcript_counts.csv",
        OUTPUT + "isoquant_prepared/isoforms.csv",
        OUTPUT + 'velocyto/merged.tagged.bam',
        OUTPUT + 'velocyto/run_velocyto.done',



rule geneformer:
    input:
        expand(OUTPUT + "geneformer_adata/{pid}.anndata.h5ad", pid=weng_ids), 
        OUTPUT + "geneformer_adata/processed.anndata.h5ad",
        OUTPUT + "geneformer_adata/processed.anndata.lt.h5ad",
        OUTPUT + "geneformer_inputs/iHSC.dataset",
        OUTPUT + "geneformer_adata/iHSC.anndata.h5ad",
        OUTPUT + "geneformer_adata/pellin.anndata.h5ad",
        OUTPUT + "geneformer_adata/tabula_sapiens.anndata.h5ad",
        
rule test:
    output:
        touch(test_file),