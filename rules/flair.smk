
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