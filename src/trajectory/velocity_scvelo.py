"""
    RNA velocity with scVelo, stochastic and dynamical models: velocity
    pseudotime, latent time, directed PAGA and a per-lane consistency check.
"""

from __future__ import annotations
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from anndata import AnnData
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from trajectory_utils import (CELL_TYPE_KEY, FOREGROUND,
                              PROJECT_ROOT, SEED, UMAP_LEGEND_FONTSIZE,
                              find_result, save_figure, save_table, set_seed,
                              style_bar, style_dark)


VELOCITY_H5AD       = PROJECT_ROOT / "data" / "interim" / "bm_velocity.h5ad"
FIT_H5AD            = PROJECT_ROOT / "data" / "interim" / "bm_velocity_fit.h5ad"
SCRIPT_NAME         = Path(__file__).stem

LANE_KEY            = "lane"
MODELS              = ("stochastic", "dynamical")
VKEYS               = {model: f"{model}_velocity" for model in MODELS}
PALANTIR_COLUMN     = "diffmap_pca_marker"
LATENT_TIME         = "latent_time"

MIN_SHARED_COUNTS   = 20
N_TOP_GENES         = 2000
N_NEIGHBORS         = 30
N_JOBS              = 4
MIN_LANE_CELLS      = 50
MIN_LIKELIHOOD      = 0.1
DEVIATION_LIMIT     = 0.3

FIGSIZE             = (16, 11)

LINEAGES            = {
    "erythroid": ("HSPC", "Erythroid early", "Erythroid late"),
    "B": ("HSPC", "Pre-B large", "Pre-B small", "B naive"),
    "myeloid": ("HSPC", "Granulocyte precursor", "Dendritic cell",
                "Monocyte classical", "Monocyte non-classical"),
    "mature lymphoid": ("CD4 T naive", "CD4 T memory", "T activated",
                        "CD8 T memory", "CD8 T cytotoxic", "NK",
                        "T proliferating"),
}

PALETTE             = plt.get_cmap("tab20")
LABEL_FONTSIZE      = 12
TITLE_FONTSIZE      = 14
TITLE_PAD           = 18
DATA_MARGIN         = 0.05
CELLTYPE_FIGSIZE    = (12, 10)
ORDERINGS_FIGSIZE   = (16, 14)
MARKERS             = {
    "palantir": ("o", "#9aa0a6", "Palantir pseudotime"),
    f"{VKEYS['stochastic']}_pseudotime": ("s", "#4fa3e0",
                                          "velocity pseudotime, stochastic"),
    f"{VKEYS['dynamical']}_pseudotime": ("D", "#f2a541",
                                         "velocity pseudotime, dynamical"),
    LATENT_TIME: ("^", "#e0607e", "latent time, dynamical"),
}


def stochastic_fit(adata: AnnData) -> None:
    """Steady-state stochastic velocity, velocity graph and pseudotime."""
    vkey = VKEYS["stochastic"]
    scv.tl.velocity(adata, mode = "stochastic", vkey = vkey)
    scv.tl.velocity_graph(adata, vkey = vkey, n_jobs = N_JOBS)
    scv.tl.velocity_pseudotime(adata, vkey = vkey)
    scv.tl.velocity_confidence(adata, vkey = vkey)


def dynamical_fit(adata: AnnData) -> None:
    """
        Full kinetic model on the stochastic velocity genes: velocity,
        velocity graph, pseudotime and latent time.
    """
    vkey = VKEYS["dynamical"]
    selected = adata.var.get(f"{VKEYS['stochastic']}_genes",
                             pd.Series(True, index = adata.var_names))
    genes = adata.var_names[selected.astype(bool)]
    scv.tl.recover_dynamics(adata, var_names = genes, n_jobs = N_JOBS)
    scv.tl.velocity(adata, mode = "dynamical", vkey = vkey)
    scv.tl.velocity_graph(adata, vkey = vkey, n_jobs = N_JOBS)
    scv.tl.velocity_pseudotime(adata, vkey = vkey)
    scv.tl.velocity_confidence(adata, vkey = vkey)
    scv.tl.latent_time(adata, vkey = vkey, min_likelihood = MIN_LIKELIHOOD)


FITTERS = {"stochastic": stochastic_fit, "dynamical": dynamical_fit}


def fitted_velocity() -> AnnData:
    """
        Velocity object with both models fitted; the saved fit is reused when
        present, otherwise the counts are filtered, smoothed and fitted.
    """
    if FIT_H5AD.exists():
        print(f"reusing {FIT_H5AD}; delete it to refit")
        return ad.read_h5ad(FIT_H5AD)
    adata = ad.read_h5ad(VELOCITY_H5AD)
    del adata.layers["ambiguous"]
    scv.pp.filter_and_normalize(adata, min_shared_counts = MIN_SHARED_COUNTS,
                                layers_normalize=["X", "spliced", "unspliced"])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = N_TOP_GENES,
                                batch_key = LANE_KEY, subset = True)
    sc.pp.neighbors(adata, n_neighbors = N_NEIGHBORS, use_rep = "X_pca",
                    random_state = SEED)
    scv.pp.moments(adata, n_neighbors = None, n_pcs = None)
    for model in MODELS:
        print(f"fitting {model} model")
        FITTERS[model](adata)
    adata.write_h5ad(FIT_H5AD, compression = "gzip")
    print(f"wrote {FIT_H5AD}")
    return adata


def cell_table(adata: AnnData) -> pd.DataFrame:
    """
        Per-cell velocity pseudotimes, latent time, velocity length and
        confidence, with cell type, lane and the Palantir pseudotime.
    """
    palantir = pd.read_csv(find_result("pt_palantir_pseudotime.csv"),
                           index_col = 0)[PALANTIR_COLUMN]
    columns = [f"{VKEYS[model]}_{suffix}" for model in MODELS
               for suffix in ("pseudotime", "length", "confidence")]
    table = adata.obs[[CELL_TYPE_KEY, LANE_KEY, *columns, LATENT_TIME]].copy()
    table["palantir"] = palantir.reindex(table.index)
    if table["palantir"].isna().any():
        raise ValueError("Palantir pseudotime missing for some cells")
    return table


def agreement(cells: pd.DataFrame,
              lineages: dict[str, tuple[str, ...]] = LINEAGES) -> pd.DataFrame:
    """
        Spearman correlation of each velocity ordering with Palantir, over all
        cells and within each lineage.
    """
    groups = {"all cells": cells, **{
        name: cells[cells[CELL_TYPE_KEY].isin(labels)]
        for name, labels in lineages.items()}}
    rows = []
    for group, frame in groups.items():
        for column in list(MARKERS)[1:]:
            rho = spearmanr(frame["palantir"], frame[column]).statistic
            rows.append({"group": group, "ordering": column,
                         "cells": len(frame), "spearman": rho})
    return pd.DataFrame(rows)


def lane_check(cells: pd.DataFrame) -> pd.DataFrame:
    """
        Per-lane medians inside each cell type and their deviation from the
        cell-type median across lanes; groups below the cell floor dropped.
    """
    metrics = [column for column in cells.columns
               if column.endswith(("_pseudotime", "_confidence"))]
    metrics.append(LATENT_TIME)
    grouped = cells.groupby([CELL_TYPE_KEY, LANE_KEY], observed = True)
    table = grouped[metrics].median()
    table["cells"] = grouped.size()
    table = table[table["cells"] >= MIN_LANE_CELLS]
    reference = table.groupby(level=CELL_TYPE_KEY, observed = True)[
        metrics].transform("median")
    deviation = (table[metrics] - reference).add_suffix("_deviation")
    return table.join(deviation).reset_index()


def paga_edges(adata: AnnData, model: str) -> pd.DataFrame:
    """
        Directed cell-type transitions from velocity-informed PAGA without a
        pseudotime prior, one row per source and target with its confidence.
    """
    scv.tl.paga(adata, groups = CELL_TYPE_KEY, vkey = VKEYS[model],
                use_time_prior = False)
    labels = adata.obs[CELL_TYPE_KEY].cat.categories
    matrix = adata.uns["paga"]["transitions_confidence"].toarray()
    targets, sources = np.nonzero(matrix)
    return pd.DataFrame({
        "model": model,
        "source": labels[sources],
        "target": labels[targets],
        "confidence": matrix[targets, sources],
    }).sort_values("confidence", ascending = False)


def gene_table(adata: AnnData) -> pd.DataFrame:
    """
        Per-gene velocity parameters of both models: steady-state ratio and
        fit of the stochastic model, kinetic rates and likelihood of the
        dynamical model.
    """
    prefixes = (VKEYS["stochastic"], "fit_")
    columns = [column for column in adata.var.columns
               if column.startswith(prefixes)]
    return adata.var[["gene_symbol", *columns]]


def set_palette(adata: AnnData) -> None:
    """Fixed tab20 colour per cell-type category, shared by every figure."""
    count = len(adata.obs[CELL_TYPE_KEY].cat.categories)
    colours = PALETTE.colors[:count]  # type: ignore
    adata.uns[f"{CELL_TYPE_KEY}_colors"] = [to_hex(colour)
                                            for colour in colours]


def dark_legend(axis: plt.Axes) -> None:  # type: ignore
    """Light legend text on the dark canvas, frame removed."""
    legend = axis.get_legend()
    if legend is None:
        return
    legend.set_frame_on(False)
    for text in legend.get_texts():
        text.set_color(FOREGROUND)


def stream_figure(adata: AnnData, model: str, prefix: str) -> None:
    """
        Velocity streamlines on the UMAP, cells coloured by cell type, with a
        margin between the data and the title.
    """
    figure, axis = plt.subplots(figsize = FIGSIZE)
    scv.pl.velocity_embedding_stream(
        adata, basis = "umap", vkey = VKEYS[model], color = CELL_TYPE_KEY,
        arrow_color = FOREGROUND, legend_loc = "right margin",
        legend_fontsize=UMAP_LEGEND_FONTSIZE, title = "", ax = axis, show = False)
    umap = adata.obsm["X_umap"]
    low, high = umap.min(axis = 0), umap.max(axis = 0)
    margin = (high - low) * DATA_MARGIN
    axis.set_xlim(low[0] - margin[0], high[0] + margin[0])
    axis.set_ylim(low[1] - margin[1], high[1] + margin[1])
    axis.set_title(f"RNA velocity, {model} model: streamlines on UMAP",
                   fontsize = TITLE_FONTSIZE, pad = TITLE_PAD)
    style_dark(figure, axis)
    dark_legend(axis)
    save_figure(figure, f"{prefix}_stream_{model}")


def paga_figure(adata: AnnData, edges: pd.DataFrame, model: str,
                prefix: str) -> None:
    """
        Velocity-directed PAGA graph at cell-type median UMAP positions,
        arrow width proportional to transition confidence.
    """
    umap = pd.DataFrame(adata.obsm["X_umap"], index = adata.obs_names,
                        columns = ("x", "y"))
    labels = adata.obs[CELL_TYPE_KEY]
    centres = umap.groupby(labels, observed = True).median()
    colours = dict(zip(labels.cat.categories,
                       adata.uns[f"{CELL_TYPE_KEY}_colors"]))
    strongest = edges["confidence"].max()
    figure, axis = plt.subplots(figsize = FIGSIZE)
    axis.scatter(umap["x"], umap["y"], s = 1, color = FOREGROUND, alpha = 0.08,
                 rasterized = True)
    widths = 0.8 + 3.7 * edges["confidence"] / strongest
    for edge, line_width in zip(edges.itertuples(), widths):
        axis.annotate("", xy = centres.loc[edge.target],
                      xytext = centres.loc[edge.source],
                      arrowprops = {"arrowstyle": "-|>", "color": FOREGROUND,
                                  "lw": line_width,
                                  "shrinkA": 9, "shrinkB": 9})
    axis.scatter(centres["x"], centres["y"], s = 160, zorder = 3,
                 color = [colours[label] for label in centres.index],
                 edgecolors = FOREGROUND, linewidths = 0.8)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_title(f"Velocity-directed PAGA, {model} model: arrows from "
                   "source to target cell type", fontsize = TITLE_FONTSIZE)
    style_dark(figure, axis)
    handles = [Line2D([], [], linestyle = "", marker = "o", markersize = 9,
                      color = colours[label], label=label)
               for label in centres.index]
    axis.legend(handles = handles, loc = "upper left", bbox_to_anchor = (1.01, 1),
                frameon = False, fontsize = LABEL_FONTSIZE - 1,
                labelcolor = FOREGROUND)
    save_figure(figure, f"{prefix}_paga_{model}")


def time_panels(adata: AnnData, cells: pd.DataFrame, prefix: str) -> None:
    """
        Two-by-two UMAP small multiples of the four orderings on one fixed
        0-1 colour scale with a shared colour bar.
    """
    umap = adata.obsm["X_umap"]
    low, high = umap.min(axis = 0), umap.max(axis = 0)
    margin = (high - low) * DATA_MARGIN
    figure, axes = plt.subplots(2, 2, figsize = ORDERINGS_FIGSIZE,
                                sharex = True, sharey = True,
                                layout = "constrained")
    figure.get_layout_engine().set(w_pad = 0.15, h_pad = 0.15, wspace = 0.04,
                                   hspace = 0.06)
    for axis, (column, (_, _, label)) in zip(axes.flat, MARKERS.items()):
        image = axis.scatter(umap[:, 0], umap[:, 1], c = cells[column], s = 2,
                             cmap = "viridis", vmin = 0, vmax = 1,
                             linewidths = 0, rasterized = True)
        axis.set_xlim(low[0] - margin[0], high[0] + margin[0])
        axis.set_ylim(low[1] - margin[1], high[1] + margin[1])
        axis.set_title(label, fontsize = TITLE_FONTSIZE, pad = 10)
        axis.set_xticks([])
        axis.set_yticks([])
        style_dark(figure, axis)
    bar = figure.colorbar(image, ax = axes, location = "right", shrink = 0.6,
                          aspect = 30, pad = 0.02)
    bar.set_label("ordering, 0 = earliest", fontsize = LABEL_FONTSIZE)
    style_bar(bar)
    figure.suptitle("Palantir pseudotime against the three velocity "
                    "orderings", color = FOREGROUND,
                    fontsize = TITLE_FONTSIZE + 2)
    save_figure(figure, f"{prefix}_orderings_umap")

def celltype_figure(cells: pd.DataFrame, prefix: str) -> None:
    """
        Median of each ordering per cell type, one row per cell type, sorted
        by Palantir median.
    """
    medians = cells.groupby(CELL_TYPE_KEY, observed = True)[
        list(MARKERS)].median().sort_values("palantir")
    rows = np.arange(len(medians))
    figure, axis = plt.subplots(figsize = CELLTYPE_FIGSIZE)
    for row in rows:
        axis.plot(medians.iloc[row].to_numpy(), [row] * len(MARKERS),
                  color = FOREGROUND, alpha = 0.25, linewidth = 1)
    for column, (marker, colour, _) in MARKERS.items():
        axis.scatter(medians[column], rows, marker = marker, color = colour,
                     s = 60, zorder = 3)
    axis.set_yticks(rows, medians.index, fontsize = LABEL_FONTSIZE)
    axis.set_xlim(0, 1)
    axis.set_xlabel("median ordering per cell type, 0 = earliest",
                    color = FOREGROUND, fontsize = LABEL_FONTSIZE)
    axis.set_title("Cell-type position under Palantir and RNA velocity",
                   fontsize = TITLE_FONTSIZE)
    style_dark(figure, axis)
    handles = [Line2D([], [], linestyle = "", marker = marker, color = colour,
                      markersize = 9, label = label)
               for marker, colour, label in MARKERS.values()]
    axis.legend(handles = handles, loc = "upper left", bbox_to_anchor = (1.01, 1),
                frameon = False, fontsize = LABEL_FONTSIZE, labelcolor = FOREGROUND)
    save_figure(figure, f"{prefix}_celltype_orderings")


def lane_figure(lanes: pd.DataFrame, metric: str, prefix: str) -> None:
    """
        Heatmap of the per-lane deviation of one metric inside each cell type,
        on a fixed symmetric scale; blank where a lane has too few cells.
    """
    matrix = lanes.pivot(index = CELL_TYPE_KEY, columns = LANE_KEY,
                         values = f"{metric}_deviation")
    figure, axis = plt.subplots(figsize = FIGSIZE)
    image = axis.imshow(matrix.to_numpy(), cmap = "RdBu_r", aspect = "auto",
                        vmin = -DEVIATION_LIMIT, vmax = DEVIATION_LIMIT)
    axis.set_xticks(range(matrix.shape[1]), matrix.columns,
                    fontsize = LABEL_FONTSIZE)
    axis.set_yticks(range(matrix.shape[0]), matrix.index,
                    fontsize = LABEL_FONTSIZE)
    axis.set_title(f"Lane deviation from the cell-type median: {metric}",
                   fontsize = TITLE_FONTSIZE)
    style_dark(figure, axis)
    bar = figure.colorbar(image, ax = axis, fraction = 0.03, pad = 0.02)
    bar.set_label("lane median minus cell-type median",
                  fontsize = LABEL_FONTSIZE)
    style_bar(bar)
    save_figure(figure, f"{prefix}_lane_{metric}")


def report_outputs(adata: AnnData, cells: pd.DataFrame, prefix: str,
                   lineages: dict[str, tuple[str, ...]]) -> pd.DataFrame:
    """
        Tables and figures shared by every velocity run: PAGA edges, cell
        orderings, Palantir agreement, lane check and gene parameters.
    """
    lanes = lane_check(cells)
    set_palette(adata)
    edge_tables = []
    for model in MODELS:
        edge_tables.append(paga_edges(adata, model))
        paga_figure(adata, edge_tables[-1], model, prefix)
        stream_figure(adata, model, prefix)
    edges = pd.concat(edge_tables)
    tables = {
        "cells": cells,
        "agreement": agreement(cells, lineages),
        "lane_check": lanes,
        "paga_edges": edges,
        "genes": gene_table(adata),
    }
    for name, table in tables.items():
        save_table(table, prefix, f"{prefix}_{name}",
                   index = name in ("cells", "genes"))
    time_panels(adata, cells, prefix)
    celltype_figure(cells, prefix)
    for metric in [column for column in lanes.columns
                   if column.endswith("_deviation")]:
        lane_figure(lanes, metric.removesuffix("_deviation"), prefix)

    pd.set_option("display.width", 250)
    print(tables["agreement"].round(3).to_string(index = False))
    print(cells.groupby(CELL_TYPE_KEY, observed = True)[list(MARKERS)]
          .median().round(3).to_string())
    print(edges.round(3).to_string(index = False))
    likely = adata.var["fit_likelihood"].dropna()
    print(f"\ndynamical fit: {likely.size} genes, median likelihood "
          f"{likely.median():.3f}, above {MIN_LIKELIHOOD}: "
          f"{(likely > MIN_LIKELIHOOD).sum()}")
    return edges


def main() -> None:
    """Fitted velocity models, their tables, figures and the lane check."""
    set_seed()
    warnings.filterwarnings("ignore", category = FutureWarning)
    scv.settings.verbosity = 2
    adata = fitted_velocity()
    report_outputs(adata, cell_table(adata), SCRIPT_NAME, LINEAGES)


if __name__ == "__main__":
    main()