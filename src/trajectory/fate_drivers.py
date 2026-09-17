"""Genes whose expression tracks each lineage's Fate probability"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.stats import t 

t_distribution = t

from trajectory_utils import (
    RESULT_DIR, load_annotated, save_figure, save_table, set_seed,
    style_bar , style_dark,
)


FATE_TAG        = "k8"
MIN_CELLS       = 50
TOP_N           = 50
HEATMAP_TOP     = 10
DRIVER_VLIM     = 0.6
DRIVER_CMAP     = "RdBu_r"
EXPRESSION_CMAP = "summer"
HEATMAP_FIGSIZE = (12, 19)
GRID_COLUMNS    = 3
PANEL_SIZE      = (7, 6)

def load_fates(index: pd.Index) -> pd.DataFrame:
    """Read the sorted fate probablity matrix and its lineage names"""
    matrix = np.load(RESULT_DIR / f"pt_cellrank_fate_{FATE_TAG}.npy")
    names = pd.read_csv(RESULT_DIR / f"pt_cellrank_fate_{FATE_TAG}_features.csv")["feature"]
    print(f"loaded fate probablities {matrix.shape} for {FATE_TAG}")
    return pd.DataFrame(matrix, index = index, columns = list(names))

def expression_matrix(adata: AnnData) -> sp.csc_matrix:
    """Return the log-normalised expression matrix in column-sparse form."""
    return sp.csc_matrix(adata.X)

def correlate(matrix: sp.csc_matrix, fates: pd.DataFrame) -> np.ndarray:
    """Return the pearson correlation of every gene with every fate probability"""
    cells = matrix.shape[0]
    values = fates.to_numpy()
    standardised = (values - values.mean(axis = 0)) / values.std(axis = 0)
    mean = np.asarray(matrix.mean(axis = 0)).ravel()
    mean_square = np.asarray(matrix.multiply(matrix).mean(axis = 0)).ravel()
    deviation = np.sqrt(np.maximum(mean_square - mean**2, 0.0))
    products = matrix.T @ standardised
    scale = cells * deviation[:, None]
    return np.divide(products, scale, out = np.zeros_like(products), where = scale > 0)

def significance(correlations: np.ndarray, cells: int) -> np.ndarray:
    """Return two-sided p-values for Pearson correlations at the given sample size"""
    statistic = correlations * np.sqrt((cells - 2) / np.maximum(1 - correlations ** 2, 1e-12))
    return (2 * t_distribution.sf(np.abs(statistic), cells - 2))

def top_drivers(correlations: pd.DataFrame, pvalues: pd.DataFrame, fraction: pd.Series) -> pd.DataFrame:
    """Return the most strongly correlated genes per lineage, strongest first."""
    blocks = []
    for lineage in correlations.columns:
        ranked = correlations[lineage].sort_values(ascending = False).head(TOP_N)
        blocks.append(pd.DataFrame({
            "lineage"      : lineage,
            "gene"         : ranked.index,
            "correlation"  : ranked.to_numpy(),
            "p_value"      : pvalues.loc[ranked.index, lineage].to_numpy(),
            "pct_cells"    : fraction.reindex(ranked.index).to_numpy(),
        }))
    return pd.concat(blocks, ignore_index = True)

def heatmap_genes(correlations: pd.DataFrame) -> list[str]:
    """Return the strongest genes per lineage, in lineage order and without repeats."""
    chosen: list[str] = []
    for lineage in correlations.columns:
        ranked = correlations[lineage].sort_values(ascending = False).head(HEATMAP_TOP)
        for gene in ranked.index:
            if gene not in chosen:
                chosen.append(gene)
    return chosen

def plot_driver_heatmap(correlations: pd.DataFrame) -> None:
    """Draw the correlation of the leading driver genes with every lineage."""
    genes   = heatmap_genes(correlations)
    frame   = correlations.loc[genes]
    figure, axis = plt.subplots(figsize = HEATMAP_FIGSIZE)
    image = axis.imshow(frame.to_numpy(), cmap = DRIVER_CMAP,
                        vmin = -DRIVER_VLIM, vmax = DRIVER_VLIM, aspect = "auto")
    axis.set_xticks(range(frame.shape[1]), frame.columns,
                          rotation = 45, ha = "right", fontsize = 6)
    axis.set_yticks(range(frame.shape[0]), frame.index, fontsize = 6)     
    axis.set_title("Fate probablity correlatoin of the leading driver genes")
    bar = figure.colorbar(image, ax = axis, shrink = 0.5, label = "Pearson correlation")
    style_bar(bar)
    style_dark(figure, axis)
    save_figure(figure, "fate_drivers_heatmap")

def plot_driver_umaps(adata: AnnData, top: pd.DataFrame) -> None: 
    """Draw the expression of the leading driver gene of every lineage on the UMAP"""
    leading = top.groupby("lineage", observed = True).head(1)
    coords = adata.obsm["X_umap"]
    rows = int(np.ceil(len(leading) / GRID_COLUMNS))
    figure, axes = plt.subplots(rows, GRID_COLUMNS,
                                figsize = (PANEL_SIZE[0] * GRID_COLUMNS, PANEL_SIZE[1] * rows),
                                squeeze = False)
    for axis, row in zip(axes.flat, leading.itertuples()):
        expression = np.asarray(adata[:, row.gene].X.todense()).ravel()
        points  = axis.scatter(coords[:, 0], coords[:, 1], c = expression,
                               cmap = EXPRESSION_CMAP, s = 2, linewidths = 0)
        axis.set_title(f"{row.lineage} - {row.gene} (r = {row.correlation:.2f})")
        axis.set_xticks([])
        axis.set_yticks([])
        bar = figure.colorbar(points, ax = axis, shrink = 0.7)
        style_bar(bar)
        style_dark(figure, axis)

    for axis in axes.flat[len(leading): ]:
        axis.set_visible(False)
    save_figure(figure, "fate_drives_umap")


def main() -> None:
    """Correlate gene expression with fate probablity and stor the driver tabels and features"""
    set_seed()
    adata = load_annotated()
    fates = load_fates(adata.obs_names)

    matrix   = expression_matrix(adata)
    detected = np.asarray((matrix > 0).sum(axis = 0)).ravel() # type: ignore
    keep     = detected >= MIN_CELLS
    genes    = adata.var_names[keep]
    print(f"testing {keep.sum()} of {adata.n_vars} genes")

    values       = correlate(matrix[:, keep], fates)
    correlations = pd.DataFrame(values, index = genes, columns = fates.columns)
    pvalues      = pd.DataFrame(significance(values, adata.n_obs), index = genes, columns = fates.columns)
    fraction     = pd.Series(detected[keep] / adata.n_obs, index = genes)


    save_table(correlations, "fate_drivers_correlations")
    save_table(top_drivers(correlations, pvalues, fraction), "fate_drivers_top", index = False)

    plot_driver_heatmap(correlations)
    plot_driver_umaps(adata, top_drivers(correlations, pvalues, fraction))



if __name__ == "__main__":
    main()



