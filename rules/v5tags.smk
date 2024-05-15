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