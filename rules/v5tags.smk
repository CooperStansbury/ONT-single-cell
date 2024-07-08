rule get_v5_tagged_reads:
    input:
        bam=OUTPUT + "merged/merged.bam",
        gtf=config['gtf_path'],
        index=OUTPUT + 'merged/merged.bam.bai',
    output:
        read_id=OUTPUT + "v5_tagged/read_ids.txt",
        read_map=OUTPUT + "v5_tagged/read_map.csv",
    conda:
        "../envs/genomic_query.yml"
    shell:
        """python scripts/extract_transcription_factors.py {input.bam} \
        {input.gtf} {output.read_id} {output.read_map}"""


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
    conda:
        "../envs/seqkit.yml"
    shell:
        """seqkit grep -n -f {input.ids} -j {threads} -o {output} {input.fastq}"""
        
        
   
rule merge_v5_reads:
    input:
        expand(OUTPUT + "v5_tagged/{sid}.v5_tagged.fastq.gz", sid=samples),
    output:
        OUTPUT + 'v5_tagged/all_reads_merged.fastq',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """cat {input} > {output}"""
        
        
rule get_reference_factors:
    input:
        '/home/cstansbu/git_repositories/ONT-single-cell/config/reprogramming_factors.fasta'
    output:
        OUTPUT + 'v5_tagged/factors.fasta'
    shell:
        "cp {input} {output}"
        

rule minimap2_index_factors:
    input:
        OUTPUT + 'v5_tagged/factors.fasta'
    output:
        OUTPUT + 'v5_tagged/factors.mmi'
    threads:
        config['threads']
    conda:
        "aligner"
    shell:
        "minimap2 -t {threads} -d {output} {input}"

        
rule align_reads_factors:
    input:
        fastq=OUTPUT + 'v5_tagged/all_reads_merged.fastq',
        ref=OUTPUT + 'v5_tagged/factors.fasta',
        refindex=OUTPUT + 'v5_tagged/factors.mmi',
    output:        
        bam=OUTPUT + 'v5_tagged/all_reads_factor_mapped.bam',
    params:
        args=config['minimap2_args'],
    threads:
        int(config['threads'])
    log:
        OUTPUT + "v5_tagged/factor_mapping.log",
    conda:
        "aligner"
    shell:
        """minimap2 -ax splice -uf --MD -t {threads} \
        {input.ref} {input.fastq} | samtools sort \
        -@ {threads} -O bam -o {output.bam} """
        
        
rule locate_factors:
    input:
        fastq=OUTPUT + "v5_tagged/{sid}.v5_tagged.fastq.gz",
        codes=OUTPUT + 'v5_tagged/factors.fasta',
    output:
        OUTPUT + 'v5_tagged/{sid}.factor_tagged.csv'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    conda:
        "../envs/seqkit.yml"
    params:
        mismatches=10
    shell:
        """seqkit locate -m {params.mismatches} -f {input.codes} {input.fastq} -j {threads} > {output} """
        
        
rule compile_factors:
    input:
        read_map=OUTPUT + "v5_tagged/read_map.csv",
        tags=expand(OUTPUT + "v5_tagged/{sid}.factor_tagged.csv", sid=samples),
    output:
        OUTPUT + 'v5_tagged/v5_result.factor_table.csv'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/compile_factor_tags.py {output} {input.tags}"""
        
                
rule locate_tags:
    input:
        fastq=OUTPUT + "v5_tagged/{sid}.v5_tagged.fastq.gz",
        codes="/home/cstansbu/git_repositories/ONT-single-cell/config/v5tag.fasta"
    output:
        OUTPUT + 'v5_tagged/{sid}.tagged.csv'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 4
    conda:
        "../envs/seqkit.yml"
    params:
        mismatches=2
    shell:
        """seqkit locate -m {params.mismatches} -f {input.codes} {input.fastq} -j {threads} > {output} """
        
        
        
        
rule compile_tags:
    input:
        read_map=OUTPUT + "v5_tagged/read_map.csv",
        tags=expand(OUTPUT + "v5_tagged/{sid}.tagged.csv", sid=samples),
    output:
        OUTPUT + 'v5_tagged/v5_result.table.csv'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/compile_v5_results.py {input.read_map} {output} {input.tags}"""