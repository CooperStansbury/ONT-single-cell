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

# structure file names and get Id lists
fastqPath = os.path.abspath(config['fastq'])
fastqDf = utils.getFastq(fastqPath)
fids = fastqDf['fileId'].to_list()



rule all:
    input:
        OUTPUT + 'references/reference.fa',
        OUTPUT + 'references/geneTable.csv',
     
        
        
rule getReference:
    input:
        refgenome=config['ref_path'],
    output:
        OUTPUT + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"
        
    
rule getFastq:
    input:
        fastq=fastqDf['filePath'].to_list(),
    output:
        expand(f"{OUTPUT}fastq/{{fileId}}.raw.fastq.gz", fileId=fids)
    run:
        from shutil import copyfile
        for i, refPath in enumerate(input.fastq):

            outPath = output[i]
            copyfile(refPath, outPath)
            
         
        
rule get_gene_table:
    input:
        annotations=config['gtf_path'],
    output:
        OUTPUT + "references/geneTable.csv"
    shell:
        "python scripts/getGeneTable.py {input} {output}"
        
