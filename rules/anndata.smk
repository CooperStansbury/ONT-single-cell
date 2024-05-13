rule make_anndata:
    input:
        counts=OUTPUT + 'merged/merged.counts.txt',
        gtf=OUTPUT + "references/annotations.gtf",
    conda:
        "../envs/scanpy.yml"
    output:
        OUTPUT + "scanpy/raw.anndata.h5ad",
    shell:
        """python scripts/make_anndata.py {input.counts} {input.gtf} {output}"""


rule process_anndata:
    input:
        OUTPUT + "scanpy/raw.anndata.h5ad",
    output:
        OUTPUT + "scanpy/processed.anndata.h5ad",
    conda:
        "../envs/scanpy.yml"
    params:
        annotations="/home/cstansbu/git_repositories/ONT-single-cell/config/gene_annotations/"
    shell:
        """python scripts/process_anndata.py {input} {output} {params.annotations}"""
        
        
        
rule cluster_anndata:
    input:
        OUTPUT + "scanpy/processed.anndata.h5ad",
    output:
        OUTPUT + "scanpy/clustered.anndata.h5ad",
    conda:
        "../envs/scanpy.yml"
    shell:
        """python scripts/cluster_cells.py {input} {output}"""