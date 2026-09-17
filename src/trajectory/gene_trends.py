"""Driver gene expression and cell type composition along pseudotime, per lineage"""

from __future__ import annotations

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


from trajectory_utils import(
    CELL_TYPE_KEY, MAIN_EMBEDDING, MAIN_ROOT_METHODS, RESULT_DIR, 
    load_annotated, run_tag, save_figure, save_table, set_seed, style_bar, style_dark,
)

FATE_TAG         = "k8"
FATE_THRESHOLD   = 0.60
MIN_CELLS        = 200
MIN_FRACTION     = 0.02
N_BINS           = 40
TOP_GENES        = 30
TREND_VLIM       = 2.0
TREND_CMAP       = "RdBu_r"
COMPOSITION_CMAP = "magma"
FIGURE_SIZE      = (12, 11)

## put in utils 
def slug(name: str) -> str: 
    """Return a file-safe form of a cell type name."""
    return name.lower().replace(" ", "_").replace("-","_")

def load_pseudotime(index: pd.Index) -> pd.Series:
    """Read the accepted Palantir pseudotime."""
    tag   = run_tag(MAIN_EMBEDDING, MAIN_ROOT_METHODS)
    frame = pd.read_csv(RESULT_DIR/ "pt_palantir_pseudotime.csv", index_col = 0)
    print(f"Using pseudotime from {tag}")
    return frame[tag].reindex(index)

def load_fates(index: pd.Index) -> pd.DataFrame:
    """Read the stored fate probablity matrix and its lineage names"""
    matrix = np.load(RESULT_DIR / f"pt_cellrank_fate_{FATE_TAG}.npy")
    names  = pd.read_csv(RESULT_DIR / f"pt_cellrank_fate_{FATE_TAG}_features.csv")["feature"]
    return pd.DataFrame(matrix, index = index, columns = list(names))


def load_drivers() -> pd.DataFrame:
    """Read the ranked driver genes of every lineage."""
    return pd.read_csv(RESULT_DIR / "fate_drivers_top.csv")

def lineage_expression(adata: AnnData, cells: pd.Index, genes: list[str]) -> pd.DataFrame:
    """Return the expression of the given genes in the given cells."""
    subset = adata[cells, genes]
    values = subset.X
    dense  = np.asarray(values.todense()) if sp.issparse(values) else np.asarray(values)
    return pd.DataFrame(dense, index = subset.obs_names, columns = genes)

def bin_by_pseudotime(frame: pd.DataFrame, pseudotime: pd.Series) -> pd.DataFrame:
    """Average every column within equal-sized pseudotime bins, earliest bin first"""
    bins   = pd.qcut(pseudotime, N_BINS, labels = False, duplicates = "drop")
    binned = frame.groupby(bins, observed = True).mean().T
    binned.index.name   = "gene"
    binned.columns.name = "bin"
    return binned

def scale_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Standardise each gene across bins so trends of different magnitude are comparable."""
    values = frame.to_numpy()
    centre = values.mean(axis = 1, keepdims = True)
    spread = values.std(axis = 1, keepdims = True)
    scaled = np.divide(values - centre, spread, out = np.zeros_like(values), where = spread > 0)
    return pd.DataFrame(scaled, index = frame.index, columns = frame.columns)

def order_by_peak(frame: pd.DataFrame) -> pd.DataFrame:
    """Sort genes by the bin in which they reach their maximum."""
    peaks = frame.to_numpy().argmax(axis = 1)
    return frame.iloc[np.argsort(peaks, kind = "stable")]

def bin_composition(labels: pd.Series, pseudotime: pd.Series) -> pd.DataFrame:
    """Return the fraction of each cell type within every pseudotime bin."""
    bins = pd.qcut(pseudotime, N_BINS, labels = False, duplicates = "drop")
    joined = pd.DataFrame({"cell_type": labels.reindex(pseudotime.index), "bin": bins})
    table = pd.crosstab(joined["cell_type"], joined["bin"], normalize = "columns")
    return table[table.max(axis = 1) >= MIN_FRACTION]

def plot_lineage(trends: pd.DataFrame, composition: pd.DataFrame, lineage: str) -> None:
    """Draw the gene trends above the cell type composition on a shared pseudotime axis"""
    figure, axes = plt.subplots(2, 1, figsize = FIGURE_SIZE, sharex = True,
                                gridspec_kw = {"height_ratios": [3, 1]})
    trend_image = axes[0].imshow(trends.to_numpy(), cmap = TREND_CMAP, 
                                 vmin = -TREND_VLIM, vmax = TREND_VLIM, aspect = "auto")
    axes[0].set_yticks(range(trends.shape[0]), trends.index, fontsize = 7)
    axes[0].set_title(f"{lineage}: driver gene trends along pseudotime")
    trend_bar = figure.colorbar(trend_image, ax = axes[0], shrink = 0.8,
                                label = "scaled expression")
    style_bar(trend_bar)

    composition_image = axes[1].imshow(composition.to_numpy(), cmap = COMPOSITION_CMAP,
                                       vmin = 0.0, vmax = 1.0, aspect = "auto")
    axes[1].set_yticks(range(composition.shape[0]), composition.index, fontsize = 7)
    axes[1].set_xlabel("pseudotime bin, earliest to latest")
    composition_bar = figure.colorbar(composition_image, ax = axes[1], shrink = 0.8,
                                      label = "fraction of bin")

    style_bar(composition_bar)

    for axis in axes:
        style_dark(figure, axis)
    save_figure(figure, f"gene_trends_{slug(lineage)}")


def main() -> None:
    """Build the pseudotime-binned gene trends and composition for every lineage"""
    set_seed()
    adata       = load_annotated()
    labels      = adata.obs[CELL_TYPE_KEY]
    pseudotime  = load_pseudotime(adata.obs_names)
    fates       = load_fates(adata.obs_names)
    drivers     = load_drivers()

    trend_blocks, composition_blocks, peaks = [], [], []

    for lineage in fates.columns:
        committed = fates.index[fates[lineage] >= FATE_THRESHOLD]
        if len(committed) < MIN_CELLS:
            print(f"{lineage}: {len(committed)} committed cells, skipped")
            continue

        ordered    = pseudotime.loc[committed].sort_values()
        genes      = list(drivers.loc[drivers["lineage"] == lineage, "gene"].head(TOP_GENES))
        expression = lineage_expression(adata, ordered.index, genes)
        trends = order_by_peak(scale_rows(bin_by_pseudotime(expression, ordered)))
        composition = bin_composition(labels, ordered) # type: ignore
        print(f"{lineage}: {len(committed)} cells, {trends.shape[1]} bins")

        trend_blocks.append(trends.stack().rename("scaled").reset_index().assign(lineage = lineage))
        composition_blocks.append(
            composition.stack().rename("fraction").reset_index().assign(lineage = lineage)
        )

        peaks.append(pd.DataFrame({
            "lineage"   : lineage,
            "gene"      : trends.index,
            "peak_bin"  : trends.to_numpy().argmax(axis = 1),
        }))

        plot_lineage(trends, composition, lineage)

    save_table(pd.concat(trend_blocks, ignore_index = True), "gene_trends_binned", 
                index = False)
    save_table(pd.concat(composition_blocks, ignore_index = True),
                "gene_trends_composition", index = False)
    save_table(pd.concat(peaks, ignore_index = True), "gene_trends_peak_order", index = False)


if __name__ == "__main__":
    main()