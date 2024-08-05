import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
import scipy 
import seaborn as sns 
import scanpy as sc
from matplotlib.colors import ListedColormap

CUSTOM_COLORS = [
    "#A30015", "#FF6978", "#FCB9B2", "#993300", "#D95F02", "#FDB462", 
    "#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A", "#FFFF99", "#6A3D9A", 
    "#CAB2D6", "#000000", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", 
    "#FDB462", "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F", "#FFFFB3", 
    "#666666", "#006400", "#00FFFF", "#B3DE69", "#00008B", "#808080", 
    "#FF0000", "#4B0082", "#EEE685", "#708090", "#008080"
]

def plot_gene_percent(adata, gene='GFI1B', threshold=0.0):
    """Plots the percentage of cells expressing a gene across different cell types.

    Args:
      adata: AnnData object containing the gene expression data.
      gene: The name of the gene to analyze.
      threshold: The expression threshold to consider a cell positive for the gene.
    """

    pdf = adata.obs.copy()
    pdf[gene] = adata[:, gene].to_df(layer='raw_counts')
    pdf[gene] = np.where(pdf[gene] > threshold, 1, 0)

    pdf = pdf.groupby('cell_label').agg(
      total=(gene, 'count'),
      positive=(gene, 'sum'),
    )

    pdf['percent'] = pdf['positive'] / pdf['total']

    sns.barplot(data=pdf, y='cell_label', x='percent', ec='k')
    sns.despine()
    plt.xlim([0, 1])
    plt.title(gene)
    plt.ylabel("")
    plt.xlabel(f"Percent {gene}+")
    plt.show()


def create_custom_colormap(colors):
    """Creates a custom Matplotlib colormap from a list of colors.

    This function takes a list of color specifications (e.g., hex codes, RGB tuples)
    and constructs a discrete colormap that can be used in Matplotlib plots.

    Args:
        colors: A list of color specifications. Each element can be:
            - A hex code string (e.g., "#FF0000" for red)
            - An RGB tuple (e.g., (1.0, 0.0, 0.0) for red)
            - Any other format accepted by Matplotlib's `ListedColormap`.

    Returns:
        matplotlib.colors.ListedColormap: The custom colormap object.

    Example:
        ```
        cmap = create_custom_colormap(["#FF0000", "#00FF00", "#0000FF"])
        plt.scatter(x, y, c=labels, cmap=cmap)
        ```
    """
    return ListedColormap(colors)



def make_colorbar(cmap='viridis', 
                  width=0.2,
                  height=2.5, 
                  title='', 
                  orientation='vertical', 
                  tick_labels=[0, 1]):
    """
    Creates and displays a standalone colorbar using Matplotlib.

    Args:
        cmap (str or matplotlib.colors.Colormap): The colormap to use for the colorbar.
        width (float): The width of the colorbar figure in inches.
        height (float): The height of the colorbar figure in inches.
        title (str): The title to display above or next to the colorbar.
        orientation (str): The orientation of the colorbar ('vertical' or 'horizontal').
        tick_labels (list of str): The labels to display at each tick on the colorbar.

    Returns:
        None: This function displays the colorbar directly using Matplotlib.

    Raises:
        ValueError: If the `orientation` is not 'vertical' or 'horizontal'.
    """
    
    a = np.array([[0, 1]])  # Dummy data for the image
    plt.figure(figsize=(width, height))
    img = plt.imshow(a, cmap=cmap)
    plt.gca().set_visible(False)  # Hide the axes of the image
    cax = plt.axes([0.1, 0.2, 0.8, 0.6])  # Define the colorbar position

    ticks = np.linspace(0, 1, len(tick_labels)) 
    cbar = plt.colorbar(
        orientation=orientation,
        cax=cax,
        label=title,
        ticks=ticks
    )

    if orientation == 'vertical':
        cbar.ax.set_yticklabels(tick_labels)
    elif orientation == 'horizontal':
        cbar.ax.set_xticklabels(tick_labels)


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
        

def generate_color_mapping(data, palette_name):
    """Generates a color mapping for unique values in a pandas Series.

    Args:
        data: A pandas Series containing categorical values.
        palette_name: The name of the seaborn color palette to use.

    Returns:
        A dictionary mapping unique values in the Series to corresponding colors.
    """

    unique_values = data.unique()
    color_palette = sns.color_palette(palette_name, n_colors=len(unique_values))
    return dict(zip(unique_values, color_palette))


def get_n_colors(n, cmap_name='viridis'):
    """Generates n evenly spaced hex color codes from a Matplotlib colormap.

    Args:
        n: The number of colors to generate.
        cmap_name: The name of the Matplotlib colormap (default: 'viridis').

    Returns:
        A list of n hex color codes (e.g., '#RRGGBB').
    """

    cmap = plt.get_cmap(cmap_name)
    return [plt.colors.rgb2hex(cmap(i / (n - 1))) for i in range(n)]


def label_point(x, y, val, ax, offset=0.05, fontsize=4):
    """Annotates points in a plot with their values.

    Args:
        x: x-coordinates of the points.
        y: y-coordinates of the points.
        val: Values to display as labels.
        ax: The Matplotlib Axes object to draw on.
        offset: Distance between the point and the label.
        fontsize: Font size of the labels.
    """

    for x_coord, y_coord, label in zip(x, y, val):
        ax.text(x_coord + offset, y_coord, str(label), fontsize=fontsize, fontweight='bold')


def plot_basin(df, x='UMAP 1', y='UMAP 2', cmap="Reds", bins=100, pthresh=.1, levels=4):
    """
    Plots a basin based on UMAP dimensionality reduction, combining scatter, 2D histogram, and density contours with automatic color selection.

    Args:
        df (pandas.DataFrame): The dataframe containing the UMAP coordinates.
        x (str): The column name for the x-axis ('UMAP 1' by default).
        y (str): The column name for the y-axis ('UMAP 2' by default).
        cmap (str): The colormap for the 2D histogram ("Reds" by default).
        bins (int): The number of bins for the 2D histogram (100 by default).
        pthresh (float): The percentile threshold for the 2D histogram (0.1 by default).
        levels (int): The number of contour levels for the density plot (4 by default).
    """

    # 2D Histogram: Shows density distribution
    sns.histplot(data=df, x=x, y=y, bins=bins, pthresh=pthresh, cmap=cmap)

    # Determine Color for Contours based on colormap
    cmap_obj = plt.cm.get_cmap(cmap)  # Get the colormap object
    contour_color = cmap_obj(0.8)  # Choose a color towards the end of the colormap
    
    # 3. Density Contours: Smoothly outlines areas of high density
    sns.kdeplot(data=df, x=x, y=y, levels=levels, color=contour_color, linewidths=0.5)


def plot_labels(data, label_col, ax=None, offset=(5, 2), text_kwargs=None):
    """
    Creates a scatterplot with labels directly on the provided axes (or a new one),
    without plotting points. Assumes UMAP coordinates for x and y axes.
    Scales plot limits to encompass the data and existing axes limits.

    Args:
        data (pd.DataFrame): DataFrame with 'UMAP 1', 'UMAP 2', and `label_col` columns.
        label_col (str): The column name for the labels.
        ax (plt.Axes, optional): The axes to plot on. If None, a new figure and axes will be created.
        offset (tuple, optional): Offset (x, y) for label positioning (default is (5, 2)).
        text_kwargs (dict, optional): Additional keyword arguments to pass to `ax.text`.
    """

    if ax is None:
        _, ax = plt.subplots()

    # Default text properties (can be overridden by text_kwargs)
    default_text_kwargs = {
        'color': 'w',
        'path_effects': [pe.withStroke(linewidth=1.5, foreground="k", alpha=0.95)],
        'verticalalignment': 'center',
        'horizontalalignment': 'center',
        'fontsize': '10',
    }
    # Combine default and user-provided kwargs
    if text_kwargs is None:
        text_kwargs = {}
    text_kwargs = {**default_text_kwargs, **text_kwargs}  # Merge dictionaries

    for _, row in data.iterrows():
        ax.text(row['UMAP 1'], row['UMAP 2'], row[label_col], **text_kwargs)

    # Adjust x and y limits to fit data and/or existing plot elements
    ax.set_xlim(
        min(data['UMAP 1'].min(), ax.get_xlim()[0]),  
        max(data['UMAP 1'].max(), ax.get_xlim()[1])   
    )
    ax.set_ylim(
        min(data['UMAP 2'].min(), ax.get_ylim()[0]),
        max(data['UMAP 2'].max(), ax.get_ylim()[1])
    )

    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    plt.show()


def create_standalone_colorbar(scatter_result, figsize=(0.15, 1.5), labelsize=8):
    """Creates a standalone colorbar figure from a matplotlib scatter plot result.

    Args:
        scatter_result: The return value of plt.scatter().
        figsize (tuple, optional): Size of the figure (width, height). Defaults to (0.15, 1.5).
        labelsize (int, optional): Font size of the colorbar tick labels. Defaults to 8.

    Returns:
        A matplotlib figure object containing the colorbar.
    """
    
    # Get colormap directly or from name
    cmap = scatter_result.cmap
    if isinstance(cmap, str):
        cmap = plt.cm.get_cmap(cmap)

    # Create figure and axes for colorbar with specified size
    fig, ax = plt.subplots(figsize=figsize)

    # Create ScalarMappable for colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=scatter_result.norm)
    sm.set_array([])  # No actual data needed, just colormapping

    # Add colorbar and customize tick labels
    cbar = plt.colorbar(sm, cax=ax)
    cbar.ax.tick_params(labelsize=labelsize)

    return fig


def plot_umap_scatter(
    adata, 
    color="CD34", 
    cmap="viridis", 
    vmin=None, 
    vmax=None, 
    label=True, 
    title=None,
    colorbar=True,
    **kwargs
):
    """
    Creates a scatterplot of UMAP data from an AnnData object with color mapping and a colorbar.

    Args:
        adata (anndata.AnnData): The AnnData object containing the UMAP embedding and color data.
        color (str, optional): Column name in `adata.var` or `adata.obs` for color mapping (default: "CD34").
        cmap (str or matplotlib.colors.Colormap, optional): Colormap to use (default: "viridis").
        vmin (float, optional): Minimum value for colormap normalization (default: None).
        vmax (float, optional): Maximum value for colormap normalization (default: None).
        label (bool, optional): Whether to add labels to the axes (default: True).
        colorbar (bool, optional): Whether to add a colorbar (default: True).
        **kwargs: Additional keyword arguments passed to `plt.scatter`.

    Returns:
        matplotlib.axes.Axes: The Axes object containing the plot.

    Raises:
        KeyError: If `color` is not found in `adata.var` or `adata.obs`.
    """
    
    
    if color in adata.var_names: # Color is a gene expression value
        df = adata.to_df()
        expression = df[color].values
    elif color in adata.obs:
        expression = adata.obs[color].values
    else:
        raise KeyError(f"Color key '{color}' not found in adata.var or adata.obs")

    # Sort points for smooth colorbar
    order = np.argsort(expression)

    # Determine vmin and vmax if not provided
    if vmin is None:
        vmin = np.min(expression)
    if vmax is None:
        vmax = np.max(expression)


    fig, ax = plt.subplots()
    
    scatter = ax.scatter(
        adata.obsm['X_umap'][:, 0][order],
        adata.obsm['X_umap'][:, 1][order],
        c=expression[order],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        **kwargs
    )
    
    ax.set_aspect('auto')

    if colorbar:
        create_standalone_colorbar(scatter)
    
    # Simplified axis label setting
    x_label, y_label = ("UMAP 1", "UMAP 2") if label else ("", "")
    ax.set(xlabel=x_label, ylabel=y_label, xticks=[], yticks=[])
    
    if not title is None:
        ax.set_title(title)
    return 