import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import seaborn as sns 
import scanpy as sc


def plot_umap_scatter(
    adata, x= "UMAP 1", 
    y="UMAP 2", color="CD34", 
    cmap="viridis",
    vmin=None,
    vmax=None,
    **kwargs
) -> plt.Axes:
    """
    Creates a scatterplot of UMAP data from an AnnData object with color mapping and a colorbar.

    Args:
        adata (anndata.AnnData): The AnnData object containing the UMAP embedding and color data.
        x (str, optional): Column name in `adata.obsm` for x-coordinates (default: "UMAP 1").
        y (str, optional): Column name in `adata.obsm` for y-coordinates (default: "UMAP 2").
        color (str, optional): Column name in `adata.var` or `adata.obs` for color mapping (default: "CD34").
        cmap (str or matplotlib.colors.Colormap, optional): Colormap to use (default: "viridis").
        vmin (float, optional): Minimum value for colormap normalization (default: None).
        vmax (float, optional): Maximum value for colormap normalization (default: None).
        **kwargs: Additional keyword arguments passed to `plt.scatter`.

    Returns:
        matplotlib.axes.Axes: The Axes object containing the plot.

    Raises:
        KeyError: If `color` is not found in `adata.var` or `adata.obs`.
    """

    if color in adata.var_names:  # Color is a gene expression value
        df = adata.to_df()
        df.columns = adata.var['gene_name'].values  # Use gene names as column labels
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

    # Create scatterplot with normalized colormap
    fig, ax = plt.subplots()
    scatter = ax.scatter(
        adata.obs[x].iloc[order],
        adata.obs[y].iloc[order],
        c=expression[order],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        **kwargs
    )

    # Add colorbar
    cbar = plt.colorbar(scatter, shrink=0.4, ax=ax)

    # Additional formatting
    sns.despine(ax=ax)
    ax.set(
        xlabel=x,
        ylabel=y,
        xticks=[],
        yticks=[]
    )

    return ax



def plot_diffusing_umap_average(df, ax=None, color="blue", base_alpha=0.8, scale=1.0, margin=0.1, num_circles=20):
    """
    Plots the average UMAP coordinates as a faded circle, scaled by variability, and adjusts plot limits.

    Args:
        df (pandas.DataFrame): The dataframe containing UMAP coordinates.
        ax (matplotlib.axes.Axes, optional): An existing axes object to plot on.
        color (str): The color of the circle.
        base_alpha (float): The base transparency for the innermost circle.
        scale (float): A scaling factor to adjust the circle size.
        margin (float): The fraction of the plot range to add as a margin around the circle.
        num_circles (int): The number of circles to create for the fading effect.
    """

    if ax is None:
        fig, ax = plt.subplots()

    # Calculate average and standard deviation of UMAP coordinates
    avg_umap1 = df['UMAP 1'].mean()
    avg_umap2 = df['UMAP 2'].mean()
    std_umap1 = df['UMAP 1'].std()
    std_umap2 = df['UMAP 2'].std()

    # Combined standard deviation
    combined_std = np.sqrt(std_umap1**2 + std_umap2**2)

    # Determine the step size for increasing the radius and decreasing alpha
    radius_step = combined_std * scale / num_circles
    alpha_step = base_alpha / num_circles

    # Create multiple circles with fading transparency
    for i in range(num_circles):
        radius = radius_step * (i + 1)
        alpha = base_alpha - alpha_step * i
        circle = plt.Circle((avg_umap1, avg_umap2), radius, color=color, alpha=alpha)
        ax.add_patch(circle)

    return ax