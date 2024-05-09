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

def min_max(v):
    return (v - v.min()) / (v.max() - v.min())

def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.05, 
                point['y'],
                str(point['val']),
                fontsize=4,
                fontweight='bold')

def labeler(x, y, label, ax, 
            kws={'fontsize': 4, 'fontweight':'bold'}):
    ax.text(x, y, str(label), **kws)
    
        
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


def get_pangloa(fpath, celltypes, cell_names):
    """A fucntion to get marker genes from panglao"""
    pdf = pd.read_csv(fpath, sep='\t')
    pdf = pdf[pdf['cell type'].isin(celltypes)]
    pdf = pdf[['official gene symbol', 'cell type']].drop_duplicates()
    pdf.columns = ['gene_name', 'cell_type']
    pdf['values'] = 1
    pdf = pd.pivot_table(pdf, 
                         columns='cell_type', 
                         index='gene_name',
                         values='values')
    pdf = pdf.fillna(0)
    pdf.columns = cell_names
    pdf = pdf.reset_index(drop=False)
    return pdf
    

def get_marker_genes(fpath, adata, celltypes, cell_names):
    """A function to load marker genes for given cell
    types """
    pdf = get_pangloa(fpath, celltypes, cell_names)
    var = adata.var.copy()
    var = var[var['gene_name'].isin(pdf['gene_name'].to_list())]
    var = var.reset_index(drop=False)
    var = pd.merge(var, pdf, 
                   how='left',
                   left_on='gene_name',
                   right_on='gene_name',)
    var = var.set_index('gene_id')
    return var


def get_genes_by_cell_type(pdf, cell_type, ui_upper=None):
    """Retrieves genes based on cell type and optional ubiquitousness index filtering.

    Args:
        pdf (pandas.DataFrame): The DataFrame containing gene data.
        cell_type (str): The cell type to filter for.
        ui_upper (float, optional): Maximum ubiquitousness index for filtering. 
                                    Defaults to None (no filtering).

    Returns:
        list: A list of official gene symbols matching the criteria.
    """

    genes = pdf[pdf['cell type'] == cell_type]

    if ui_upper is not None:
        genes = genes[genes['ubiquitousness index'] < ui_upper]

    return genes['official gene symbol'].to_list()

