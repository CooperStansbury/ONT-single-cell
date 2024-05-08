
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