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