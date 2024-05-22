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