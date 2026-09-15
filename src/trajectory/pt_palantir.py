"""Palantir pseudotime and fate probablities across the embedding and root grid."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import palantir 
import pandas as pd
from anndata import AnnData
from scipy.stats import spearmanr

from trajectory_utils import (
    CELL_TYPE_KEY, SEED, 
    attach_embeddings, embedding_key, load_annotated, run_grid, run_tag,
    save_figure, save_matrix, save_table, select_root, set_seed, 
    style_bar, style_dark,
)

BASE_EMBEDDING   = {"diffmap_pca" : "pca", "diffmap_scvi" : "scvi"}
N_DIFF_COMPS     = 10
NUM_WAYPOINTS    = 1200
PALANTIR_KNN     = 30
PSEUDOTIME_CMAP  = "plasma"
GRID_FIGSIZE     = (16, 12)
CELLTYPE_FIGSIZE = (9, 11)

def multiscale_space(adata: AnnData, embedding: str) -> str:
    """ 
        Build the Palantir diffusion map of the matching base embedding and
        return its ``obsm`` key.
    """
    suffix = f"_{embedding}"
    palantir.utils.run_diffusion_maps(
        adata, n_components = N_DIFF_COMPS, seed = SEED,
        pca_key     = embedding_key(BASE_EMBEDDING[embedding]),
        kernel_key  = f"DM_Kernel{suffix}",
        sim_key     = f"DM_Similarity{suffix}",
        eigval_key  = f"DM_EigenValues{suffix}",
        eigvec_key  = f"DM_EigenVectors{suffix}",
    )
    out_key = f"DM_EigenVector_multiscaled{suffix}"
    palantir.utils.determine_multiscale_space(
        adata, n_eigs = N_DIFF_COMPS,
        eigval_key  = f"DM_EigenValues{suffix}",
        eigvec_key  = f"DM_EigenVectors{suffix}",
        out_key     = out_key,
    )
    print(f"diffusion map built for {embedding}")
    return out_key

def run_palantir_once(adata: AnnData, embedding: str, root_method: str, space_key: str):
    """Run Palantir for one embedding and root selection, returing its result object and root barcode. """
    root = select_root(adata, root_method, embedding)
    tag  = run_tag(embedding, root_method)
    result = palantir.core.run_palantir(
        adata, 
        early_cell      = root,
        knn             = PALANTIR_KNN, 
        num_waypoints   = NUM_WAYPOINTS,
        seed            = SEED, 
        eigvec_key      = space_key,
        pseudo_time_key = f"pseudotime_{tag}",
        entropy_key     = f"entropy_{tag}",
        fate_prob_key   = f"fate_{tag}",
        waypoints_key   = f"waypoints_{tag}"
    )
    print(f"{tag}: {result.branch_probs.shape[1]} terminal states")
    return result, root

def run_correlation(frame: pd.DataFrame) -> pd.DataFrame:
    """Return the Spearman correlation of pseudotime between every pair of runs."""
    matrix, _ = spearmanr(frame.to_numpy())
    return pd.DataFrame(matrix, index = frame.columns, columns = frame.columns)

def plot_umap_grid(adata: AnnData, frame: pd.DataFrame) -> None:
    """Draw pseudotime on the UMAP for every run, on one shared zero to one scale."""
    coords = adata.obsm["X_umap"]
    figure, axes = plt.subplots(2, 2, figsize = GRID_FIGSIZE)
    for axis, tag in zip(axes.flat, frame.columns):
        points = axis.scatter(coords[:, 0], coords[:, 1], c = frame[tag],
                              cmap = PSEUDOTIME_CMAP, s = 2 , linewidths = 0,
                              vmin = 0.0, vmax = 1.0
                              )
        axis.set_title(tag)
        axis.set_xticks([])
        axis.set_yticks([])
        style_dark(figure, axis)
    bar = figure.colorbar(points, ax = axes, shrink = 0.6, label = "pseudotime")
    style_bar(bar)
    save_figure(figure, "pt_palantir_umap_grid")

def plot_celltype_pseudotime(means: pd.DataFrame) -> None:
    """Draw the mean pseudotime of each cell type per run, ordered by the primary run."""
    ordered = means.sort_values(means.columns[0])
    figure, axis = plt.subplots(figsize = CELLTYPE_FIGSIZE)
    image = axis.imshow(ordered.to_numpy(), cmap = PSEUDOTIME_CMAP,
                        vmin = 0.0, vmax = 1.0, aspect = "auto")
    axis.set_xticks(
        range(ordered.shape[1]), ordered.columns, rotation=45, ha="right", fontsize=8
    )
    axis.set_yticks(range(ordered.shape[0]), ordered.index, fontsize=8)
    for row, column in np.ndindex(ordered.shape):
        axis.text(
            column,
            row,
            f"{ordered.iat[row, column]:.2f}",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
        )
    axis.set_title("Mean Palantir pseudotime")
    bar = figure.colorbar(image, ax=axis, shrink=0.8, label="pseudotime")
    style_bar(bar)
    style_dark(figure, axis)
    save_figure(figure, "pt_palantir_pseudotime_by_celltype")

def main() -> None:
    """Run Palantir over every embedding and root pair, then store the tables and figures."""
    set_seed()
    adata = load_annotated()
    attach_embeddings(adata)
    attach_embeddings(adata, tuple(BASE_EMBEDDING.values()))
    labels = adata.obs[CELL_TYPE_KEY]

    spaces, pseudotime, entropy, roots, terminals = {}, {}, {}, [], []

    for embedding, root_method in run_grid():
        if embedding not in spaces:
            spaces[embedding] = multiscale_space(adata, embedding)
        tag = run_tag(embedding, root_method)
        result, root = run_palantir_once(
            adata, embedding, root_method, spaces[embedding]
        )

        pseudotime[tag] = result.pseudotime.reindex(adata.obs_names)
        entropy[tag] = result.entropy.reindex(adata.obs_names)
        roots.append({"run": tag, "root": root, "cell_type": labels.loc[root]})

        fates = result.branch_probs.reindex(adata.obs_names)
        names = [f"{labels.loc[cell]} | {cell}" for cell in fates.columns]
        save_matrix(fates.to_numpy(), f"pt_palantir_fate_{tag}", names)
        terminals.extend(
            {"run": tag, "terminal": cell, "cell_type": labels.loc[cell]}
            for cell in fates.columns
        )

    frame = pd.DataFrame(pseudotime)
    means = frame.groupby(labels, observed=True).mean()

    save_table(frame, "pt_palantir_pseudotime")
    save_table(pd.DataFrame(entropy), "pt_palantir_entropy")
    save_table(pd.DataFrame(roots), "pt_palantir_roots", index=False)
    save_table(pd.DataFrame(terminals), "pt_palantir_terminal_states", index=False)
    save_table(means, "pt_palantir_pseudotime_by_celltype")
    save_table(run_correlation(frame), "pt_palantir_run_correlation")

    plot_umap_grid(adata, frame)
    plot_celltype_pseudotime(means)


if __name__ == "__main__":
    main()