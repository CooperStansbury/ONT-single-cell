rule get_fastq:
    input:
        fastq=input_df['file_path'].to_list()
    output:
        input_names
    run:
        from shutil import copyfile
        for i, refPath in enumerate(input.fastq):

            outPath = output[i]
            copyfile(refPath, outPath)


rule raw_report:
    input:
        expand(f"{OUTPUT}fastq/{{sid}}.raw.fastq.gz", sid=samples),
    output:
        OUTPUT + "reports/seqkit_stats/raw_report.txt",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule demultiplex:
    input:
        fastq=OUTPUT + "fastq/{sid}.raw.fastq.gz",
        whitelist=config['barcode_whitelist'],
    output:
        touch(OUTPUT + "demultiplex/{sid}.done"),
        OUTPUT + 'demultiplex/{sid}.emtpy_bc_list.csv',
        OUTPUT + 'demultiplex/{sid}.knee_plot.png',
        OUTPUT + 'demultiplex/{sid}.matched_reads.fastq.gz',
        OUTPUT + 'demultiplex/{sid}.putative_bc.csv',
        OUTPUT + 'demultiplex/{sid}.summary.txt',
        OUTPUT + 'demultiplex/{sid}.whitelist.csv',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] 
    params:
        expected=config['expected_cells'],
        output_prefix=lambda wildcards: OUTPUT + "demultiplex/" + wildcards.sid + ".", 
    log:
        OUTPUT + "demultiplex/{sid}.log",
    shell:
        """blaze --expect-cells {params.expected} \
        --output-prefix {params.output_prefix} --threads {threads} \
        --full-bc-whitelist {input.whitelist} {input.fastq} """


rule demultiplexed_report:
    input:
        flags=expand(f"{OUTPUT}demultiplex/{{sid}}.done", sid=samples),
    output:
        OUTPUT + "reports/seqkit_stats/demultiplexed_report.txt",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    params:
        files=expand(f"{OUTPUT}demultiplex/{{sid}}.matched_reads.fastq.gz", sid=samples),
    shell:
        """seqkit stats -a -b -j {threads} {params.files} -o {output}"""


rule fastqc:
    input:
        OUTPUT + "demultiplex/{sid}.matched_reads.fastq.gz",
    output:
        html=OUTPUT + "reports/fastqc/{sid}.report.html",
        zip=OUTPUT + "reports/fastqc/{sid}.report.zip"
    params: "--quiet"
    log:
        OUTPUT + "reports/fastqc/{sid}.log"
    threads:
        config['threads'] // 4
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    wrapper:
        "v1.29.0/bio/fastqc"


rule minimap2_align:
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


rule merge_bam:
    input:
        expand(f"{OUTPUT}mapping/{{sid}}.tagged.bam", sid=samples),
    output:
        OUTPUT + 'merged/merged.bam',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        int(config['threads'])
    shell:
        """samtools merge -@ {threads} {output} {input} """


rule samtools_index_merged:
    input:
        OUTPUT + 'merged/merged.bam',
    output:
        OUTPUT + 'merged/merged.bam.bai'
    threads:
        int(config['threads'])
    shell:
        """samtools index -@ {threads} {input}"""


rule samtools_stats_merged:
    input:
        OUTPUT + 'merged/merged.bam',
    output:
        OUTPUT + 'merged/merged.stats'
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