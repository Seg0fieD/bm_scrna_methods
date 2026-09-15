"""PAGA topology of the annotated cell types on each diffusion map embedding."""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

import numpy as np
import pandas as pd 
import scanpy as sc
from anndata import AnnData

from trajectory_utils import (
    CELL_TYPE_KEY, EMBEDDINGS, N_NEIGHBORS, FIGSIZE, SEED,
    attach_embeddings, embedding_key, load_annotated, save_figure, save_table, set_seed
)

PLOT_THRESHOLD  = 0.10
EDGE_CMAP       = "cividis"
HEATMAP_CMAP    = "magma"
HEATMAP_FIGSIZE = (11, 9)
BACKGROUND      = "#12141a"
FOREGROUND      = "#e6e6e6"

def ensure_categorical(adata: AnnData) -> None:
    """Cast the cell type column to categorical dtype, as required for PAGA grouping"""
    if not isinstance(adata.obs[CELL_TYPE_KEY].dtype, pd.CategoricalDtype):
        adata.obs[CELL_TYPE_KEY] = adata.obs[CELL_TYPE_KEY].astype("category")

def cell_counts(adata :AnnData) -> pd.DataFrame:
    """Return the number of cells per cel type, in category order."""
    categories = adata.obs[CELL_TYPE_KEY].cat.categories
    counts = adata.obs[CELL_TYPE_KEY].value_counts().reindex(categories)
    return pd.DataFrame({"cells": counts})

def compute_paga(adata: AnnData, embedding: str) -> None:
    """Build the neighbour graph on the named embedding and run PAGA over the cell types"""
    sc.pp.neighbors(adata, 
                    n_neighbors  = N_NEIGHBORS,
                    use_rep      = embedding_key(embedding),
                    random_state = SEED,
                    )
    sc.tl.paga(adata, groups = CELL_TYPE_KEY)
    adata.uns["paga"].pop("pos", None)
    print(f"PAGA computed on {embedding}")

def connectivity_frame(adata: AnnData, key: str) -> pd.DataFrame:
    """Return a PAGA connectivity matrix as a dense frame labelled by cell type."""
    categories = adata.obs[CELL_TYPE_KEY].cat.categories
    matrix = adata.uns["paga"][key].toarray()
    return pd.DataFrame(matrix, index = categories, columns = categories)

def edge_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Return the non-zero edges of a symmetric connectivity matrix, strongest first."""
    matrix = frame.to_numpy()
    rows, cols = np.triu_indices_from(matrix, k = 1)
    table = pd.DataFrame(
        {
            "source"       : frame.index[rows],
            "target"       : frame.columns[cols],
            "connectivity" : matrix[rows, cols],
        }
    ) 
    table = table[table["connectivity"] > 0]
    return table.sort_values("connectivity", ascending = False, ignore_index = True)

def recolor_edges(axis: Axes) -> None:
    """
        Recolor the drawn edges by line width, so connection strength reads as color as well as thickness
    """
    palette = plt.get_cmap(EDGE_CMAP)
    for collection in axis.collections:
        if not isinstance(collection, LineCollection):
            continue
        widths = np.asarray(collection.get_linewidth(), dtype = float)
        if widths.size == 0:
            continue
        norm = Normalize(vmin = widths.min(), vmax = widths.max())
        collection.set_color(palette(0.30 + 0.70 * norm(widths)))

def style_dark(figure: Figure, axis: Axes) -> None:
    """Apply a dark canvas with ligt text, ticks and spines"""
    figure.patch.set_facecolor(BACKGROUND)
    axis.set_facecolor(BACKGROUND)
    axis.title.set_color(FOREGROUND)
    axis.tick_params(colors = FOREGROUND)
    for spine in axis.spines.values():
        spine.set_color(FOREGROUND)

def plot_graph(adata: AnnData, embedding: str) -> None:
    """Draw the PAGA graph colored by cell type and write it to the figure directory."""
    figure, axis = plt.subplots(figsize = FIGSIZE)
    sc.pl.paga( adata, color = CELL_TYPE_KEY, 
                threshold = PLOT_THRESHOLD,  layout = "fr",
                random_state = SEED, fontsize = 9, fontoutline = 1, 
                ax = axis, node_size_scale = 1.5, edge_width_scale = 0.7,
                min_edge_width = 0.3, max_edge_width = 6.0, frameon = False,
                show = False )
    recolor_edges(axis)
    axis.set_title(f"PAGA topology on {embedding}")
    style_dark(figure, axis)
    save_figure(figure, f"paga_topology_{embedding}")

def plot_heatmap(frame: pd.DataFrame, embedding: str) -> None: 
    """Draw the connectivity matrix as a heatmap on a fixed zero to one scale."""
    labels       = list(frame.index)
    figure, axis = plt.subplots(figsize = HEATMAP_FIGSIZE)
    image = axis.imshow(frame.to_numpy(), cmap = HEATMAP_CMAP, vmin = 0.0, vmax = 1.0)
    axis.set_xticks(range(len(labels)), labels, rotation = 90, fontsize = 8)    
    axis.set_yticks(range(len(labels)), labels, fontsize = 8)        
    axis.set_title(f"PAGA connectivity on {embedding}")
    bar = figure.colorbar(image, ax = axis, shrink = 0.8, label = "connectivity")
    # bar.ax.yaxis.label.set_color(FOREGROUND)
    # bar.ax.tick_params(colors = FOREGROUND)
    # style_dark(figure, axis)
    save_figure(figure, f"paga_topology_{embedding}_heatmap")

def main() -> None:
    """Compute and store the PAGA topology for every embedding under comparison."""
    set_seed()
    adata = load_annotated()
    attach_embeddings(adata)
    ensure_categorical(adata)
    save_table(cell_counts(adata), "paga_topology_cell_counts")


    for embedding in EMBEDDINGS:
        compute_paga(adata, embedding)
        full = connectivity_frame(adata, "connectivities")
        tree = connectivity_frame(adata, "connectivities_tree")
        save_table(full, f"paga_topology_{embedding}_connectivities")
        save_table(tree, f"paga_topology_{embedding}_tree")
        save_table(edge_table(full), f"paga_topology_{embedding}_edges", index = False)
        plot_graph(adata, embedding)
        plot_heatmap(full, embedding)

if __name__ == "__main__":
    main()
