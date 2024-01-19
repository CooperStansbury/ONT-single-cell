rule get_reference:
    input:
        refgenome=config['ref_path'],
    output:
        OUTPUT + 'references/reference.fa.gz'
    shell:
        "cp {input} {output}"


rule get_annotations:
    input:
        refgenome=config['gtf_path'],
    output:
        OUTPUT + 'references/annotations.gtf'
    shell:
        "cp {input} {output}"


rule prep_reference:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        ref=OUTPUT + 'references/reference.fa',
        flag=touch(OUTPUT + 'references/reference.done')
    shell:
        "cat {input} | gzip -d > {output.ref}"
        
            
rule get_gene_table:
    input:
        annotations=config['gtf_path'],
    output:
        OUTPUT + "references/geneTable.csv"
    shell:
        "python scripts/getGeneTable.py {input} {output}"


rule minimap2_index:
    input:
        refgenome=OUTPUT + 'references/reference.fa.gz'
    output:
        OUTPUT + 'references/reference.mmi'
    threads:
        config['threads']
    shell:
        "minimap2 -t {threads} -d {output} {input.refgenome}"

