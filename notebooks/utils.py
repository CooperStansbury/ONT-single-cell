import os
import pandas as pd
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
import io

def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.05, 
                point['y'],
                str(point['val']),
                fontsize=4,
                fontweight='bold')

def ncolor(n, cmap='viridis'):
    cmap = matplotlib.cm.get_cmap(cmap)
    arr = np.linspace(0, 1, n)
    return [matplotlib.colors.rgb2hex(cmap(x)) for x in arr] 


def get_colors(data, cmap):
    """A function to return seaborn colormap
    dict from a colum """
    color_list = sns.palettes.color_palette(cmap,
                                            data.nunique(), 
                                            as_cmap=False)
    return color_list


def parseKEGG(pathId):
    genes = []
    results = REST.kegg_get(pathId).read()
    current_section = None
    for line in results.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section

        if current_section == "GENE":
            linesplit = line[12:].split("; ")
            gene_identifiers = linesplit[0]
            gene_id, gene_symbol = gene_identifiers.split()
    
            if not gene_symbol in genes:
                genes.append(gene_symbol)
    return genes

def getPathname(pathId):
    """A function to return the legg pathname"""
    result = REST.kegg_list(pathId).read()
    return result.split("\t")[1].split("-")[0].strip()


def makeColorbar(cmap, width, hieght, title, orientation, tickLabels):
    a = np.array([[0,1]])
    plt.figure(figsize=(width, hieght))
    img = plt.imshow(a, cmap=cmap)
    plt.gca().set_visible(False)
    cax = plt.axes([0.1, 0.2, 0.8, 0.6])
    ticks = np.linspace(0,1 , len(tickLabels))
    cbar = plt.colorbar(orientation=orientation, 
                        cax=cax, 
                        label=title,
                        ticks=ticks)

    if orientation == 'vertical':
        cbar.ax.set_yticklabels(tickLabels)
    else:
        cbar.ax.set_xticklabels(tickLabels)

        
def _normalize_data(X, counts, after=None, copy=False):
    X = X.copy() if copy else X
    if issubclass(X.dtype.type, (int, np.integer)):
        X = X.astype(np.float32)  # TODO: Check if float64 should be used
    else:
        counts_greater_than_zero = counts[counts > 0]

    after = np.median(counts_greater_than_zero, axis=0) if after is None else after
    counts += counts == 0
    counts = counts / after
    if scipy.sparse.issparse(X):
        sparsefuncs.inplace_row_scale(X, 1 / counts)
    elif isinstance(counts, np.ndarray):
        np.divide(X, counts[:, None], out=X)
    else:
        X = np.divide(X, counts[:, None])  # dask does not support kwarg "out"
    return X


def normalize(df, target_sum=1):
    """A function to normalize spots """
    index = df.index
    columns = df.columns
    X = df.to_numpy().copy()
    counts_per_cell = X.sum(1)
    counts_per_cell = np.ravel(counts_per_cell)
    cell_subset = counts_per_cell > 0
    Xnorm = _normalize_data(X, counts_per_cell, target_sum)
    
    ndf = pd.DataFrame(Xnorm, columns=columns, index=index)
    return ndf