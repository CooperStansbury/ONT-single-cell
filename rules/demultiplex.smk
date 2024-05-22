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
        expand(OUTPUT + "fastq/{sid}.raw.fastq.gz", sid=samples),
    output:
        OUTPUT + "reports/seqkit_stats/raw_report.txt",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    conda:
        "../envs/seqkit.yml"
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
    conda:
        "blaze"
    log:
        OUTPUT + "demultiplex/{sid}.log",
    shell:
        """blaze --expect-cells {params.expected} \
        --output-prefix {params.output_prefix} --threads {threads} \
        --full-bc-whitelist {input.whitelist} {input.fastq} """


rule demultiplexed_report:
    input:
        flags=expand(OUTPUT + "demultiplex/{sid}.done", sid=samples),
    output:
        OUTPUT + "reports/seqkit_stats/demultiplexed_report.txt",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    conda:
        "../envs/seqkit.yml"
    params:
        files=expand(OUTPUT + "demultiplex/{sid}.matched_reads.fastq.gz", sid=samples),
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
