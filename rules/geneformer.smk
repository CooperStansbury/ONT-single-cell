rule ihsc_to_geneformer:
    input:
        adata=OUTPUT + "scanpy/clustered.anndata.h5ad",
        gene_table=OUTPUT + 'references/geneTable.csv',
    output:
        OUTPUT + "geneformer_adata/iHSC.anndata.h5ad"
    conda:
        "geneformer"
    params:
        gt_key="gene_name",
        gene_id="gene_name",
        gene_column_type="gene_name",
        gene_index="gene_name",
        counts_layer="raw_counts",
    shell:
        """python scripts/to_geneformer_adata.py {input.adata} \
        {input.gene_table} {params.gt_key} {params.gene_id} \
        {params.gene_column_type} {params.gene_index} \
        {params.counts_layer} {output} 
        """
        
        
rule pellin_to_geneformer:
    input:
        adata="/nfs/turbo/umms-indikar/shared/projects/HSC/data/pellin_2019/pellin.anndata.h5ad",
        gene_table=OUTPUT + 'references/geneTable.csv',
    output:
        OUTPUT + "geneformer_adata/pellin.anndata.h5ad"
    conda:
        "geneformer"
    params:
        gt_key="gene_name",
        gene_id="gene_name",
        gene_column_type="gene_name",
        gene_index="gene_name",
        counts_layer="raw_counts",
    shell:
        """python scripts/to_geneformer_adata.py {input.adata} \
        {input.gene_table} {params.gt_key} {params.gene_id} \
        {params.gene_column_type} {params.gene_index} \
        {params.counts_layer} {output} 
        """
    
    
rule weng_to_geneformer:
    input:
        adata="/nfs/turbo/umms-indikar/shared/projects/HSC/data/weng_2024/scanpy_objects/{pid}.h5ad",
        gene_table=OUTPUT + 'references/geneTable.csv',
    output:
        OUTPUT + "geneformer_adata/{pid}.anndata.h5ad"
    conda:
        "geneformer"
    params:
        gt_key="gene_name",
        gene_id="gene_name",
        gene_column_type="gene_name",
        gene_index="gene_name",
        counts_layer="raw_counts",
    shell:
        """python scripts/to_geneformer_adata.py {input.adata} \
        {input.gene_table} {params.gt_key} {params.gene_id} \
        {params.gene_column_type} {params.gene_index} \
        {params.counts_layer} {output} 
        """
        
        
rule tabula_to_geneformer:
    input:
        adata="/nfs/turbo/umms-indikar/shared/projects/HSC/data/tabula_sapiens/tabula_sapiens_filtered.h5ad",
        gene_table=OUTPUT + 'references/geneTable.csv',
    output:
        OUTPUT + "geneformer_adata/tabula_sapiens.anndata.h5ad"
    conda:
        "geneformer"
    params:
        gt_key="gene_name",
        gene_id="gene_name",
        gene_column_type="gene_name",
        gene_index="gene_name",
        counts_layer="raw_counts",
    shell:
        """python scripts/to_geneformer_adata.py {input.adata} \
        {input.gene_table} {params.gt_key} {params.gene_id} \
        {params.gene_column_type} {params.gene_index} \
        {params.counts_layer} {output} 
        """
        
rule merge_geneformer_inputs:
    input:
        gene_table=OUTPUT + 'references/geneTable.csv',
        adatas=expand(OUTPUT + "geneformer_adata/{gid}.anndata.h5ad", gid=geneformer_ids)
    output:
        OUTPUT + "geneformer_adata/merged.anndata.h5ad"
    conda:
        "geneformer"
    shell:
        """python scripts/merge_geneformer_adata.py {input.gene_table} {output} {input.adatas}"""
        
        
rule process_merged_adata_gene:
    input:
        OUTPUT + "geneformer_adata/merged.anndata.h5ad"
    output:
        OUTPUT + "geneformer_adata/processed.anndata.h5ad"
    conda:
        "../envs/scanpy.yml"
    params:
        annotations="/home/cstansbu/git_repositories/ONT-single-cell/config/gene_annotations/"
    shell:
        """python scripts/process_merged_adata.py {input} {output} {params.annotations}"""
    
        
rule tokenize:
    input:
        adata= OUTPUT + "geneformer_adata/processed.anndata.h5ad",
        gene_lengths="/home/cstansbu/git_repositories/Geneformer/geneformer/gene_median_dictionary.pkl",
        tokens="/home/cstansbu/git_repositories/Geneformer/geneformer/token_dictionary.pkl",
    output:
        directory(OUTPUT + "geneformer_inputs/iHSC.dataset")
    conda:
        "geneformer"
    shell:
        """python scripts/tokenize_data.py {input.adata} {output} {input.gene_lengths} {input.tokens} """
    