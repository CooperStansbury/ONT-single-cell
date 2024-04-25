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
        expand(f"{OUTPUT}velocyto/{{sid}}.done", sid=samples),


rule all_flair:
    input:
        expand(f"{OUTPUT}flair/align/{{sid}}.bam", sid=samples),
        expand(f"{OUTPUT}flair/correct/{{sid}}_all_corrected.bed", sid=samples),
        expand(f"{OUTPUT}flair/collapse/{{sid}}.isoforms.bed", sid=samples),


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


rule htseq_count_individual:
    input:
        bam=OUTPUT + 'mapping/{sid}.tagged.bam',
        annotations=config['gtf_path'],
    output:
        OUTPUT + "counts/individual/{sid}.counts.txt"
    conda:
        "envs/htseq_count.yml"
    params:
        d=int(config['umi_distance'])
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        "htseq-count-barcodes --nonunique all {input.bam} {input.annotations} > {output}"


rule make_anndata:
    input:
        counts=OUTPUT + "counts/counts.txt",
        gtf=OUTPUT + "references/annotations.gtf",
    output:
        OUTPUT + "scanpy/anndata.h5ad",
    shell:
        """python scripts/make_anndata.py {input.counts} {input.gtf} {output}"""


rule process_anndata:
    input:
        anndata=OUTPUT + "scanpy/anndata.h5ad",
    output:
        OUTPUT + "scanpy/anndata.processed.h5ad",
    params:
        config=str(BASE_DIR) + "/config/config.yaml"
    shell:
        """python scripts/process_anndata.py {input.anndata} {params.config} {output}"""


rule prepare_nanocount:
    input:
        flag=OUTPUT + "demultiplex/{sid}.done",
        ref=OUTPUT + 'references/transcripts.fa',
    output:        
        bam=OUTPUT + 'nanocount/mapping/{sid}.bam',
    params:
        fastq=OUTPUT + "demultiplex/{sid}.matched_reads.fastq.gz",
    threads:
        int(config['threads'])
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    log:
        OUTPUT + "nanocount/mapping/{sid}.log",
    shell:
        """minimap2 -ax map-ont -N 10 -t {threads} \
        {input.ref} {params.fastq} | samtools view \
        -bh > {output.bam} """


rule run_nanocount:
    input:
        OUTPUT + 'nanocount/mapping/{sid}.bam',
    output:
        counts=OUTPUT + 'nanocount/{sid}.tx_counts.tsv',
        bam=OUTPUT + 'nanocount/{sid}.filtered.bam',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """NanoCount -i {input} -b {output.bam} \
        --extra_tx_info -o {output.counts} """


rule merge_nanocount_bam:
    input:
        expand(f"{OUTPUT}nanocount/mapping/{{sid}}.bam", sid=samples),
    output:
        OUTPUT + 'nanocount/merged/merged.bam',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        int(config['threads'])
    shell:
        """samtools merge -@ {threads} {output} {input} """


rule run_nanocount_merged:
    input:
        OUTPUT + 'nanocount/merged/merged.bam',
    output:
        counts=OUTPUT + 'nanocount/merged/merged_tx_counts.tsv',
        bam=OUTPUT + 'nanocount/merged/merged.filtered.bam',
    shell:
        """NanoCount -i {input} -b {output.bam} \
        --extra_tx_info -o {output.counts} """


rule run_velocyto:
    input:
        bam=OUTPUT + 'mapping/{sid}.tagged.bam',
        gtf=config['gtf_path'],
    output:
        flag=touch(OUTPUT + 'velocyto/{sid}.done'),
    conda:
        "envs/velocyto.yml"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 2
    params:
        prefix=lambda wildcards: OUTPUT + "velocyto/" + wildcards.sid + "/",
    shell:
        """velocyto run --multimap \
        --samtools-threads {threads} \
        -o {params.prefix} {input.bam} {input.gtf} """


rule get_v5_tagged_reads:
    input:
        bam=OUTPUT + "merged/merged.bam",
        gtf=config['gtf_path'],
    output:
        OUTPUT + "v5_tagged/read_ids.txt"
    shell:
        """python scripts/extract_TF_reads.py {input.bam} {input.gtf} {output}"""


rule extract_read_v5_tag:
    input:
        ids=OUTPUT + "v5_tagged/read_ids.txt",
        fastq=OUTPUT + "demultiplex/{sid}.matched_reads.fastq.gz",
    output:
        OUTPUT + "v5_tagged/{sid}.v5_tagged.fastq.gz"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit grep -n -f {input.ids} -j {threads} -o {output} {input.fastq}"""


rule merge_v5_tags:
    input:
        expand(f"{OUTPUT}v5_tagged/{{sid}}.v5_tagged.fastq.gz", sid=samples)
    output:
        OUTPUT + "v5_tagged/v5_tagged.fastq.gz"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """cat {input} > {output} """


rule flair_align:
    input:
        fastq=OUTPUT + "demultiplex/{sid}.matched_reads.fastq.gz",
        ref=OUTPUT + 'references/reference.fa',
        refindex=OUTPUT + 'references/reference.mmi',
    output:
        bam=OUTPUT + 'flair/align/{sid}.bam',
        bai=OUTPUT + 'flair/align/{sid}.bam.bai',
        bed=OUTPUT + 'flair/align/{sid}.bed',
    conda:
        "envs/flair.yml"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    params:
        prefix=lambda wildcards: OUTPUT + "flair/align/" + wildcards.sid, 
    threads:
        int(config['threads']) // 4
    shell:
        """flair align -g {input.ref} -r {input.fastq} \
        --threads {threads} --output {params.prefix} """


rule flair_correct:
    input:
        bed=OUTPUT + 'flair/align/{sid}.bed',
        gtf=OUTPUT + 'references/annotations.gtf',
        ref=OUTPUT + 'references/reference.fa',
    output:
        bed=OUTPUT + 'flair/correct/{sid}_all_corrected.bed',
    conda:
        "envs/flair.yml"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    params:
        prefix=lambda wildcards: OUTPUT + "flair/correct/" + wildcards.sid, 
    threads:
        int(config['threads']) // 4
    shell:
        """flair correct -q {input.bed} -f {input.gtf} \
        -g {input.ref} --threads {threads} --output {params.prefix} """


rule flair_collapse:
    input:
        fastq=OUTPUT + "demultiplex/{sid}.matched_reads.fastq.gz",
        bed=OUTPUT + 'flair/correct/{sid}_all_corrected.bed',
        gtf=OUTPUT + 'references/annotations.gtf',
        ref=OUTPUT + 'references/reference.fa',
    output:
        bam=OUTPUT + 'flair/collapse/{sid}.isoforms.bed',
    conda:
        "envs/flair.yml"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    params:
        prefix=lambda wildcards: OUTPUT + "flair/collapse/" + wildcards.sid + ".", 
        tmp=OUTPUT + "flair"
    threads:
        int(config['threads']) // 4
    shell:
        """flair collapse --threads {threads} \
        -g {input.ref} -q {input.bed} -r {input.fastq} \
        --isoformtss --annotation_reliant generate \
        --generate_map --trust_ends \
        --temp_dir {params.tmp} --output {params.prefix} """
