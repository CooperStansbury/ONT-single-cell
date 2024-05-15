rule align_reads:
    input:
        flag=OUTPUT + "demultiplex/{sid}.done",
        ref=OUTPUT + 'references/reference.fa.gz',
        refindex=OUTPUT + 'references/reference.mmi',
    output:        
        bam=OUTPUT + 'mapping/{sid}.bam',
    params:
        args=config['minimap2_args'],
        fastq=OUTPUT + "demultiplex/{sid}.matched_reads.fastq.gz",
    threads:
        int(config['threads'])
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    log:
        OUTPUT + "mapping/{sid}.log",
    shell:
        """minimap2 {params.args} -t {threads} \
        {input.ref} {params.fastq} | samtools sort \
        -@ {threads} -O bam -o {output.bam} """


rule samtools_index:
    input:
        OUTPUT + 'mapping/{sid}.bam'
    output:
        OUTPUT + 'mapping/{sid}.bam.bai'
    conda:
        "../envs/samtools.yml"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        int(config['threads'])
    shell:
        """samtools index -@ {threads} {input}"""


rule bamtools_stats:
    input:
        OUTPUT + 'mapping/{sid}.bam'
    output:
        OUTPUT + 'reports/bamstats/{sid}.bamstats'
    params:
        "-insert" # optional summarize insert size data
    log:
        OUTPUT + "reports/bamstats/{sid}.log"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        int(config['threads']) // 4
    wrapper:
        "v2.1.1/bio/bamtools/stats"


rule tag_bam:
    input:
        OUTPUT + 'mapping/{sid}.bam'
    output:
        bam=OUTPUT + 'mapping/{sid}.tagged.bam',
        records=OUTPUT + 'mapping/{sid}.records.csv',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/tag_bam.py {input} {output.bam} {output.records}"""


rule samtools_index_tagged:
    input:
        OUTPUT + 'mapping/{sid}.tagged.bam'
    output:
        OUTPUT + 'mapping/{sid}.tagged.bam.bai'
    conda:
        "../envs/samtools.yml"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        int(config['threads'])
    shell:
        """samtools index -@ {threads} {input}"""


rule merge_bam:
    input:
        expand(OUTPUT + "mapping/{sid}.tagged.bam", sid=samples),
    output:
        OUTPUT + 'merged/merged.bam',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    conda:
        "../envs/samtools.yml"
    threads:
        int(config['threads'])
    shell:
        """samtools merge -@ {threads} {output} {input} """


rule samtools_index_merged:
    input:
        OUTPUT + 'merged/merged.bam',
    output:
        OUTPUT + 'merged/merged.bam.bai'
    conda:
        "../envs/samtools.yml"
    threads:
        int(config['threads'])
    shell:
        """samtools index -@ {threads} {input}"""


rule samtools_stats_merged:
    input:
        OUTPUT + 'merged/merged.bam',
    output:
        OUTPUT + 'merged/merged.stats'
    conda:
        "../envs/samtools.yml"
    threads:
        int(config['threads']) // 4
    shell:
        """samtools stats -@ {threads} {input} > {output}"""


rule bamtools_stats_merged:
    input:
        OUTPUT + 'merged/merged.bam',
    output:
       OUTPUT + 'merged/merged.bamstats'
    params:
        "-insert" # optional summarize insert size data
    log:
        OUTPUT + "merged/merged.bamstats.log"
    threads:
        int(config['threads']) // 4
    wrapper:
        "v2.1.1/bio/bamtools/stats"


rule htseq_count:
    input:
        bam=OUTPUT + 'merged/merged.bam',
        annotations=config['gtf_path'],
    output:
        OUTPUT + "merged/merged.counts.txt"
    conda:
        "../envs/htseq_count.yml"
    params:
        d=int(config['umi_distance'])
    shell:
        "htseq-count-barcodes --nonunique all {input.bam} {input.annotations} > {output}"


rule htseq_count_individual:
    input:
        bam=OUTPUT + 'mapping/{sid}.tagged.bam',
        annotations=config['gtf_path'],
    output:
        OUTPUT + "individual_counts/{sid}.counts.txt"
    conda:
        "../envs/htseq_count.yml"
    params:
        d=int(config['umi_distance'])
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        "htseq-count-barcodes --nonunique all {input.bam} {input.annotations} > {output}"