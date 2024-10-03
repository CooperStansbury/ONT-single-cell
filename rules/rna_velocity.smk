

rule prep_velocyto:
    input:
        OUTPUT + 'merged/merged.bam',
    output:
        bam=OUTPUT + 'velocyto/merged.tagged.bam',
        bcfile=OUTPUT + 'velocyto/bcfile.tsv',
    conda:
        'velocyto'
    shell:
        """python scripts/prep_velocyto.py {input} {output.bam} {output.bcfile}"""



rule run_velocyto:
    input:
        bam=OUTPUT + 'velocyto/merged.tagged.bam',
        gtf=OUTPUT + 'references/annotations.gtf',
        bcfile=OUTPUT + 'velocyto/bcfile.tsv',
    output:
        flag=touch(OUTPUT + 'velocyto/run_velocyto.done'),
    conda:
        'velocyto'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 2
    params:
        outdir=lambda wildcards: OUTPUT + "velocyto/" 
    shell:
        """velocyto run {input.bam} {input.gtf} \
        --bcfile {input.bcfile} --samtools-threads {threads} \
        --outputfolder {params.outdir}"""
