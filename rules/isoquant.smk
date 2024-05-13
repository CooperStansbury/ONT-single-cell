
rule run_isoquant:
    input:
        ref=OUTPUT + 'references/reference.fa',
        gtf=OUTPUT + 'references/annotations.gtf',
        bam=OUTPUT + 'merged/merged.bam',
        bam_index=OUTPUT + 'merged/merged.bam.bai',
    output:
        log=OUTPUT + "isoquant/isoquant.log",
        db=OUTPUT + "isoquant/annotations.db",
    params:
        prefix=OUTPUT + "isoquant"
    conda:
        "../envs/isoquant.yml"
    shell:
        """isoquant.py --reference {input.ref} --genedb {input.gtf} --bam {input.bam} --data_type 'nanopore' --complete_genedb -o {params.prefix}"""

