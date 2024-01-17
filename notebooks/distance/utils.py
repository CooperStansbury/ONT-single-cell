import pandas as pd
import numpy as np
import anndata as ad
import os
import sys
from collections import Counter
import gget
import scipy
from sklearn import preprocessing
from sklearn.feature_extraction.text import TfidfTransformer

import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns


def get_stacked_cell_data(adata, genes='marker'):
    """Get the genes """
    if genes is None:
        genes = adata.var_names.copy()
    elif genes == 'marker':
        genes = adata.var.copy()
        genes = genes[genes['is_fb_marker'] | genes['is_hsc_marker']]
        genes = genes.index.to_list()
    
    """Get the cells"""
    our_idx = adata[adata.obs['set'] == 'our_data'].obs.index.to_list()
    fb_idx = adata[adata.obs['celltype'] == 'FB'].obs.index.to_list()
    hsc_idx = adata[adata.obs['celltype'] == 'HSC'].obs.index.to_list()
    
    # # get the numer of fibroblast signatures to sample
    n_fb_to_sample = len(our_idx) - len(hsc_idx)
    
    fb_idx = np.random.choice(fb_idx, n_fb_to_sample, replace=False)
    fb_idx = list(fb_idx)
    
    idx = our_idx + fb_idx + hsc_idx
    
    """ extract the data """
    pdf = adata[idx, genes]
    return pdf