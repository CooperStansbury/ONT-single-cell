rule build_db:
    input:
        OUTPUT + 'references/annotations.gtf',
    output:
        OUTPUT + "isoquant/annotations.db",
    conda:
        "../envs/isoquant.yml"
    shell:
        """python scripts/build_isoquant_db.py {input} {output}"""


rule run_isoquant:
    input:
        ref=OUTPUT + 'references/reference.fa',
        db=OUTPUT + "isoquant/annotations.db",
        bam=OUTPUT + "mapping/{sid}.tagged.bam",
        bam_index=OUTPUT + "mapping/{sid}.tagged.bam.bai",
    output:
        directory(OUTPUT + "isoquant/{sid}"),
        touch(OUTPUT + "isoquant/{sid}.done"),
    params:
        outdir=OUTPUT + "isoquant",
        prefix=lambda wildcards: wildcards.sid
    conda:
        "../envs/isoquant.yml"
    threads:
        config['threads']
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """isoquant.py --reference {input.ref} --genedb {input.db} \
        --threads {threads} --count_exons --prefix {params.prefix} \
        --bam {input.bam} --data_type 'nanopore' \
        --bam_tags 'CB,UB,RD' \
        --complete_genedb -o {params.outdir}"""
        
        
rule merge_gene_counts:
    input:
        expand(OUTPUT + "isoquant/{sid}/{sid}.gene_counts.tsv", sid=samples),
    output:
        OUTPUT + "isoquant_prepared/gene_counts.csv"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/merge_isoquant.py {output} {input}"""
        
        
rule merge_transcript_counts:
    input:
        expand(OUTPUT + "isoquant/{sid}/{sid}.transcript_counts.tsv", sid=samples),
    output:
        OUTPUT + "isoquant_prepared/transcript_counts.csv"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/merge_isoquant.py {output} {input}"""
        
        
rule make_prepped_isoquant:
    input:
        gene_table=OUTPUT + 'references/geneTable.csv',
        gene_path=OUTPUT + "isoquant_prepared/gene_counts.csv",
        trx_path=OUTPUT + "isoquant_prepared/transcript_counts.csv",
    output:
        OUTPUT + "isoquant_prepared/isoforms.csv"
    conda:
        "bioinf"
    shell:
        """python scripts/process_isoquant.py {input.gene_table} {input.gene_path} {input.trx_path} {output}"""
    
        