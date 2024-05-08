




rule run_velocyto:
    input:
        bam=OUTPUT + 'mapping/{sid}.tagged.bam',
        gtf=config['gtf_path'],
    output:
        flag=touch(OUTPUT + 'velocyto/{sid}.done'),
    conda:
        "envs/velocyto.yml"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    threads:
        config['threads'] // 2
    params:
        prefix=lambda wildcards: OUTPUT + "velocyto/" + wildcards.sid + "/",
    shell:
        """velocyto run --multimap \
        --samtools-threads {threads} \
        -o {params.prefix} {input.bam} {input.gtf} """
