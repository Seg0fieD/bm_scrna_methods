"""PAGA topology of the annotated cell types on each diffusion map embedding."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import scanpy as sc
from anndata import AnnData


from trajectory_utils import (
    CELL_TYPE_KEY, EMBEDDINGS, N_NEIGHBORS, FIGSIZE, SEED,
    attach_embeddings, embedding_key, load_annotated, save_figure, save_table, set_seed
)

PLOT_THRESHOLD = 0.05

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


def plot_graph(adata: AnnData, embedding: str) -> None:
    """Draw the PAGA graph colored by cell type and write it to the figure directory."""
    figure, axis = plt.subplots(figsize = FIGSIZE)
    sc.pl.paga(adata, color = CELL_TYPE_KEY, 
               threshold = PLOT_THRESHOLD,  layout = "fr",
               random_state = SEED, fontsize = 8, ax = axis,
                node_size_scale = 1.5, edge_width_scale = 0.7,
                 show = False )
    axis.set_title(f"PAGA topology on {embedding}")
    save_figure(figure, f"paga_topology_{embedding}")


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

if __name__ == "__main__":
    main()
