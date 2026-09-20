"""Cell-Rank cross-check of the trajectory: marostates, terminal states and fate probablities"""
from __future__ import annotations

import cellrank as cr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from trajectory_utils import (
    CELL_TYPE_KEY,  MAIN_EMBEDDING, MAIN_ROOT_METHODS, RESULT_DIR, SEED,
    attach_embeddings, embedding_key, load_annotated, run_tag, save_figure,
    save_matrix, save_table, set_seed, style_bar, style_dark,
)

BASE_EMBEDDING      = "pca"
GRAPH_NEIGHBORS     = 30
STATE_COUNTS        = (8, 5, 4)
SCHUR_METHOD        = "krylov"
PSEUDOTIME_WEIGHT   = 0.8
FATE_CMAP           = "plasma"
UMAP_CMAP           = "summer"
GRID_COLUMNS        = 3
PANEL_SIZE          = (6, 5)
HEATMAP_FIGSIZE     = (8, 9)


def attach_pseudotime(adata: AnnData) -> str:
    """Attach the accepted Palantir pseudotime to ``obs`` and return its key."""
    tag = run_tag(MAIN_EMBEDDING, MAIN_ROOT_METHODS)
    frame = pd.read_csv(RESULT_DIR / "pt_palantir_pseudotime.csv", index_col = 0)
    key = f"pseudotime_{tag}"
    adata.obs[key] = frame[tag].reindex(adata.obs_names).to_numpy()
    print(f"attached pseudotime from {tag}")
    return key


def build_kernel(adata: AnnData, time_key: str):
    """Return a transition matrix combining pseudotime direction with graph connectivity."""
    sc.pp.neighbors(
        adata,
        n_neighbors  = GRAPH_NEIGHBORS,
        use_rep      = embedding_key(BASE_EMBEDDING),
        random_state = SEED,
    )
    directed = cr.kernels.PseudotimeKernel(
        adata, time_key = time_key
    ).compute_transition_matrix()
    undirected = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    return PSEUDOTIME_WEIGHT * directed + (1 - PSEUDOTIME_WEIGHT) * undirected


def fit_estimator(kernel, n_states: int):
    """Fit the estimator at the given macrostate resolution and compute fate probabilities."""
    estimator = cr.estimators.GPCCA(kernel)
    estimator.compute_schur(n_components=n_states + 2, method=SCHUR_METHOD) # type: ignore
    estimator.compute_macrostates(n_states=n_states, cluster_key=CELL_TYPE_KEY)
    estimator.predict_terminal_states()
    estimator.compute_fate_probabilities() # type: ignore
    return estimator


def composition(states: pd.Series, labels: pd.Series) -> pd.DataFrame:
    """Return the cell type composition of each state, aligned on cell barcodes."""
    joined = pd.DataFrame({"state": states, "cell_type": labels})
    return pd.crosstab(joined["state"], joined["cell_type"])


def fate_frame(estimator, index: pd.Index) -> pd.DataFrame:
    """Return the fate probabilities as a frame labelled by lineage."""
    fates = estimator.fate_probabilities
    return pd.DataFrame(np.asarray(fates), index = index, columns = list(fates.names))


def plot_fate_grid(adata: AnnData, frame: pd.DataFrame, tag: str) -> None:
    """Draw every fate probability on the UMAP, on one shared zero to one scale."""
    coords = adata.obsm["X_umap"]
    rows = int(np.ceil(frame.shape[1] / GRID_COLUMNS))
    figure, axes = plt.subplots(
        rows,
        GRID_COLUMNS,
        figsize=(PANEL_SIZE[0] * GRID_COLUMNS, PANEL_SIZE[1] * rows),
        squeeze=False,
    )
    for axis, name in zip(axes.flat, frame.columns):
        points = axis.scatter(
            coords[:, 0],
            coords[:, 1],
            c          = frame[name],
            cmap       = UMAP_CMAP,
            s          = 2,
            linewidths = 0,
            vmin       = 0.0,
            vmax       = 1.0,
        )
        axis.set_title(name)
        axis.set_xticks([])
        axis.set_yticks([])
        style_dark(figure, axis)
    for axis in axes.flat[frame.shape[1] :]:
        axis.set_visible(False)
    bar = figure.colorbar(points, ax = axes, shrink = 0.6, label = "fate probability")
    style_bar(bar)
    save_figure(figure, f"pt_cellrank_fate_umap_{tag}")


def plot_fate_heatmap(means: pd.DataFrame, tag: str) -> None:
    """Draw the mean fate probability of each cell type per lineage, on a fixed scale."""
    figure, axis = plt.subplots(figsize=HEATMAP_FIGSIZE)
    image = axis.imshow(
        means.to_numpy(), cmap = FATE_CMAP, vmin = 0.0, vmax = 1.0, aspect = "auto"
    )
    axis.set_xticks(
        range(means.shape[1]), means.columns, rotation = 45, ha = "right", fontsize = 8
    )
    axis.set_yticks(range(means.shape[0]), means.index, fontsize=8)
    for row, column in np.ndindex(means.shape):
        axis.text(
            column,
            row,
            f"{means.iat[row, column]:.2f}",
            ha       = "center",
            va       = "center",
            fontsize = 7,
            color    = "white",
        )
    axis.set_title(f"Mean CellRank fate probability, {tag}")
    bar = figure.colorbar(image, ax = axis, shrink = 0.8, label = "fate probability")
    style_bar(bar)
    style_dark(figure, axis)
    save_figure(figure, f"pt_cellrank_fate_by_celltype_{tag}")


def main() -> None:
    """Run the CellRank cross-check at each macrostate resolution and store its results."""
    set_seed()
    adata = load_annotated()
    attach_embeddings(adata, (BASE_EMBEDDING,))
    labels = adata.obs[CELL_TYPE_KEY]
    kernel = build_kernel(adata, attach_pseudotime(adata))

    for n_states in STATE_COUNTS:
        tag = f"k{n_states}"
        estimator = fit_estimator(kernel, n_states)

        save_table(
            composition(estimator.macrostates, labels), f"pt_cellrank_macrostates_{tag}"
        )
        save_table(
            composition(estimator.terminal_states, labels),
            f"pt_cellrank_terminal_states_{tag}",
        )

        frame = fate_frame(estimator, adata.obs_names)
        means = frame.groupby(labels, observed=True).mean()
        save_matrix(frame.to_numpy(), f"pt_cellrank_fate_{tag}", list(frame.columns))
        save_table(means, f"pt_cellrank_fate_by_celltype_{tag}")

        plot_fate_grid(adata, frame, tag)
        plot_fate_heatmap(means, tag)

if __name__ == "__main__":
    main()
    