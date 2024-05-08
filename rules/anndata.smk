rule make_anndata:
    input:
        counts=OUTPUT + "counts/counts.txt",
        gtf=OUTPUT + "references/annotations.gtf",
    output:
        OUTPUT + "scanpy/anndata.h5ad",
    shell:
        """python scripts/make_anndata.py {input.counts} {input.gtf} {output}"""


rule process_anndata:
    input:
        anndata=OUTPUT + "scanpy/anndata.h5ad",
    output:
        OUTPUT + "scanpy/anndata.processed.h5ad",
    params:
        config=str(BASE_DIR) + "/config/config.yaml"
    shell:
        """python scripts/process_anndata.py {input.anndata} {params.config} {output}"""