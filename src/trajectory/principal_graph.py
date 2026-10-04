"""Principal tree fitted over the main diffusion embedding, with its own
pseudotime as a second opinion on the Palantir ordering."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scFates as scf
from anndata import AnnData
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr

from trajectory_utils import (
    CELL_TYPE_KEY, FOREGROUND, MAIN_EMBEDDING, MAIN_ROOT_METHODS, SEED,
    attach_embeddings, find_result, load_annotated, run_tag, save_figure,
    save_table, select_root, set_seed, style_bar, style_dark,
)

RUN      = run_tag(MAIN_EMBEDDING, MAIN_ROOT_METHODS)
UMAP_KEY = "X_umap"
REP_KEY  = "X_principal_rep"

NDIMS_REP  = 10
N_NODES    = 50
PPT_LAMBDA = 1000.0
PPT_SIGMA  = 0.1
PPT_NSTEPS = 200
MAX_TIPS   = 14

EMBED_SIZE  = (14, 11)
BOX_SIZE    = (12, 9)
SQUARE_SIZE = (10, 9)
PLOT_RIGHT  = 0.76

POINT_SIZE      = 1.5
POINT_ALPHA     = 0.5
EDGE_COLOUR     = "#e6e6e6"
TIP_COLOUR      = "#f2c14e"
FORK_COLOUR     = "#ef6f6c"
ROOT_COLOUR     = "#ffffff"
CATEGORY_CMAP   = "tab20"
PSEUDOTIME_CMAP = "viridis"
HEX_CMAP        = "magma"
HEX_GRID        = 60

TITLE_SIZE  = 17
LABEL_SIZE  = 12
TICK_SIZE   = 11
LEGEND_SIZE = 11
MARKER_SIZE = 9

GRAPH_MARKERS = (
    ("root", "*", ROOT_COLOUR, 320, MARKER_SIZE + 6, "root (stage root cell)"),
    ("tip", "^", TIP_COLOUR, 95, MARKER_SIZE + 2, "tip (lineage end)"),
    ("fork", "D", FORK_COLOUR, 130, MARKER_SIZE, "branch point"),
)


def fit_principal_tree(adata: AnnData) -> None:
    """Simple principal tree over the leading diffusion components, scaled by
    one global standard deviation so sigma and lambda carry a fixed meaning;
    a fit with more tips than ``MAX_TIPS`` is rejected."""
    matrix = np.asarray(adata.obsm[f"X_{MAIN_EMBEDDING}"], dtype = float)
    matrix = matrix[:, :NDIMS_REP]
    adata.obsm[REP_KEY] = matrix / matrix.std()
    scf.tl.tree(
        adata, method = "ppt", Nodes = N_NODES, use_rep = REP_KEY,
        ppt_lambda = PPT_LAMBDA, ppt_sigma = PPT_SIGMA,
        ppt_nsteps = PPT_NSTEPS, device = "cpu", seed = SEED,
    )
    tips = np.size(adata.uns["graph"]["tips"])
    forks = np.size(adata.uns["graph"]["forks"])
    print(f"fitted tree: {tips} tips, {forks} branch points")
    if tips > MAX_TIPS:
        raise ValueError(
            f"{tips} tips exceeds MAX_TIPS = {MAX_TIPS}; raise PPT_LAMBDA"
        )


def order_cells(adata: AnnData, barcode: str) -> None:
    """Tree rooted at the principal point holding the stage's root cell, and
    every cell projected onto it."""
    position = adata.obs_names.get_loc(barcode)
    scf.tl.root(adata, int(adata.obsm["X_R"][position].argmax()))
    try:
        scf.tl.pseudotime(adata, n_jobs = 1, n_map = 1, seed = SEED)
    except IndexError:
        # scFates fails while colouring a branch that received no milestone,
        # after every per-cell result has already been written.
        if "t" not in adata.obs:
            raise


def graph_layout(adata: AnnData) -> tuple[pd.DataFrame, np.ndarray]:
    """Principal points on the UMAP as assignment-weighted cell means, with
    pseudotime, branch, attached cell mass and tip, fork and root flags; and
    the tree edges as index pairs."""
    graph = adata.uns["graph"]
    weights = np.asarray(adata.obsm["X_R"], dtype = float)
    mass = weights.sum(axis = 0)
    xy = weights.T @ adata.obsm[UMAP_KEY] / np.where(mass > 0, mass, 1.0)[:, None]
    info = pd.DataFrame(graph["pp_info"])
    index = np.arange(len(mass))
    nodes = pd.DataFrame(
        {
            "umap_1"     : xy[:, 0],
            "umap_2"     : xy[:, 1],
            "pseudotime" : info["time"].to_numpy(),
            "segment"    : info["seg"].astype(str).to_numpy(),
            "cells"      : mass,
            "tip"        : np.isin(index, graph["tips"]),
            "fork"       : np.isin(index, graph["forks"]),
            "root"       : index == int(graph["root"]),
        }
    ).rename_axis("node")
    return nodes, np.argwhere(np.triu(graph["B"], k = 1))


def draw_tree(axis, nodes: pd.DataFrame, edges: np.ndarray) -> None:
    """Principal tree over an embedding panel: edges, principal points, and
    the root, tips and branch points marked."""
    xy = nodes[["umap_1", "umap_2"]].to_numpy()
    axis.add_collection(
        LineCollection(
            xy[edges], colors = EDGE_COLOUR, linewidths = 1.1, alpha = 0.85, # type: ignore
            zorder = 3,
        )
    )
    axis.scatter(*xy.T, s = 16, color = EDGE_COLOUR, linewidths = 0, zorder = 4)
    for flag, marker, colour, size, _, _ in GRAPH_MARKERS:
        mask = nodes[flag].to_numpy()
        axis.scatter(
            *xy[mask].T, s = size, marker = marker, color = colour,
            linewidths = 0, zorder = 6 if flag == "root" else 5,
        )


def legend_column(figure, handles: list[Line2D], title: str, top: bool) -> None:
    """Legend in the right-hand margin, anchored to its top or bottom edge."""
    legend = figure.legend(
        handles        = handles,
        title          = title,
        loc            = "upper left" if top else "lower left",
        bbox_to_anchor = (PLOT_RIGHT + 0.025, 0.95 if top else 0.06),
        frameon        = False,
        fontsize       = LEGEND_SIZE,
        title_fontsize = LEGEND_SIZE + 2,
        labelspacing   = 0.55,
        handletextpad  = 0.9,
        borderaxespad  = 0.0,
    )
    legend.get_title().set_color(FOREGROUND)
    for text in legend.get_texts():
        text.set_color(FOREGROUND)


def tree_figure(
    coords: np.ndarray, nodes: pd.DataFrame, edges: np.ndarray,
    colour_by: pd.Series, key: str, title: str, name: str,
) -> None:
    """UMAP coloured by a categorical or continuous cell annotation, with the
    principal tree drawn over it: the colour key at the top of the right-hand
    margin and the tree legend at its foot."""
    figure, axis = plt.subplots(figsize = EMBED_SIZE)
    figure.subplots_adjust(right = PLOT_RIGHT)
    style = {
        "s": POINT_SIZE, "alpha": POINT_ALPHA, "linewidths": 0,
        "rasterized": True, "zorder": 1,
    }
    if pd.api.types.is_numeric_dtype(colour_by):
        points = axis.scatter(
            *coords.T, c = colour_by.to_numpy(), cmap = PSEUDOTIME_CMAP,
            vmin = 0.0, vmax = 1.0, **style,
        )
        bar = figure.colorbar(
            points, cax = figure.add_axes((PLOT_RIGHT + 0.03, 0.42, 0.018, 0.5))
        )
        bar.set_label(key, color = FOREGROUND, fontsize = LABEL_SIZE)
        style_bar(bar)
    else:
        categories = colour_by.astype("category").cat.categories
        colours = plt.get_cmap(CATEGORY_CMAP)(
            np.linspace(0.0, 1.0, len(categories))
        )
        for colour, category in zip(colours, categories):
            mask = (colour_by == category).to_numpy()
            axis.scatter(*coords[mask].T, color = colour, **style)
        legend_column(
            figure,
            [
                Line2D(
                    [], [], marker = "o", ls = "none", color = colour,
                    markersize = MARKER_SIZE, label = category,
                )
                for colour, category in zip(colours, categories)
            ],
            key,
            top = True,
        )
    draw_tree(axis, nodes, edges)
    legend_column(
        figure,
        [
            Line2D(
                [], [], marker = marker, ls = "none", color = colour,
                markersize = size, label = label,
            )
            for _, marker, colour, _, size, label in GRAPH_MARKERS
        ]
        + [
            Line2D(
                [], [], marker = "o", color = EDGE_COLOUR,
                markersize = MARKER_SIZE - 4, label = "principal point",
            )
        ],
        "Fitted tree",
        top = False,
    )
    axis.set_title(title, color = FOREGROUND, fontsize = TITLE_SIZE, pad = 14)
    style_dark(figure, axis)
    axis.set_axis_off()
    save_figure(figure, name)


def plot_celltype_pseudotime(graph_time: pd.Series, labels: pd.Series) -> None:
    """Principal tree pseudotime distribution of every cell type, ordered by
    median."""
    order = graph_time.groupby(labels).median().sort_values().index
    figure, axis = plt.subplots(figsize = BOX_SIZE)
    parts = axis.boxplot(
        [graph_time[labels == name].to_numpy() for name in order],
        orientation = "horizontal", showfliers = False, patch_artist = True,
    )
    colours = plt.get_cmap(PSEUDOTIME_CMAP)(np.linspace(0.0, 1.0, len(order)))
    for patch, colour in zip(parts["boxes"], colours):
        patch.set(facecolor = colour, edgecolor = FOREGROUND)
    for key in ("whiskers", "caps", "medians"):
        plt.setp(parts[key], color = FOREGROUND)
    axis.set_yticks(range(1, len(order) + 1), order, fontsize = TICK_SIZE)
    axis.set_xlim(0.0, 1.0)
    axis.set_xlabel(
        "pseudotime along the tree, 0 at the root", fontsize = LABEL_SIZE
    )
    axis.set_title(
        "Principal tree pseudotime by cell type", fontsize = TITLE_SIZE
    )
    style_dark(figure, axis)
    axis.xaxis.label.set_color(FOREGROUND)
    figure.tight_layout()
    save_figure(figure, "principal_graph_pseudotime_by_celltype")


def plot_agreement(
    graph_time: pd.Series, palantir_time: pd.Series, spearman: float
) -> None:
    """Cell density of the two pseudotimes against each other on fixed unit
    axes, with the identity line for reference."""
    figure, axis = plt.subplots(figsize = SQUARE_SIZE)
    hexes = axis.hexbin(
        palantir_time, graph_time, gridsize = HEX_GRID, cmap = HEX_CMAP,
        bins = "log", extent = (0.0, 1.0, 0.0, 1.0), mincnt = 1,
    )
    axis.plot([0, 1], [0, 1], color = FOREGROUND, linewidth = 0.8, alpha = 0.6)
    axis.set(xlim = (0.0, 1.0), ylim = (0.0, 1.0))
    axis.set_xlabel("Palantir pseudotime", fontsize = LABEL_SIZE)
    axis.set_ylabel("principal tree pseudotime", fontsize = LABEL_SIZE)
    axis.set_title(
        f"Pseudotime agreement, Spearman {spearman:.2f}", fontsize = TITLE_SIZE
    )
    style_dark(figure, axis)
    axis.xaxis.label.set_color(FOREGROUND)
    axis.yaxis.label.set_color(FOREGROUND)
    bar = figure.colorbar(hexes, ax = axis, fraction = 0.04, pad = 0.02)
    bar.set_label("cells", color = FOREGROUND, fontsize = LABEL_SIZE)
    style_bar(bar)
    figure.tight_layout()
    save_figure(figure, "principal_graph_agreement")


def segment_table(
    segments: pd.Series, labels: pd.Series, graph: dict
) -> pd.DataFrame:
    """Every branch of the tree with its end points, length, cell count and
    cell type composition; branches that received no cells are kept."""
    network = pd.DataFrame(graph["pp_seg"]).astype({"n": str}).set_index("n")
    counts = pd.crosstab(segments, labels)
    table = network.rename(columns = {"d": "length"}).join(counts).fillna(0)
    table.insert(3, "cells", counts.sum(axis = 1).reindex(table.index).fillna(0))
    return table.rename_axis("segment")


def celltype_table(
    graph_time: pd.Series, palantir_time: pd.Series, labels: pd.Series
) -> pd.DataFrame:
    """Per cell type: median position under each ordering, the rank of that
    median, and the shift in rank between the two."""
    frame = pd.DataFrame({"graph": graph_time, "palantir": palantir_time})
    table = frame.groupby(labels).agg(
        cells           = ("graph", "size"),
        median_graph    = ("graph", "median"),
        median_palantir = ("palantir", "median"),
    )
    table["rank_graph"] = table["median_graph"].rank().astype(int)
    table["rank_palantir"] = table["median_palantir"].rank().astype(int)
    table["rank_shift"] = table["rank_graph"] - table["rank_palantir"]
    return table.sort_values("median_graph")


def agreement_table(
    graph_time: pd.Series, palantir_time: pd.Series
) -> pd.DataFrame:
    """Rank and linear correlation between the two pseudotimes over all
    cells."""
    tests = {"spearman": spearmanr, "pearson": pearsonr}
    rows = [
        (name, *test(graph_time, palantir_time)) for name, test in tests.items()
    ]
    return pd.DataFrame(rows, columns = ["statistic", "value", "p_value"])


def main() -> None:
    """Principal tree fit, cell ordering along it, and comparison of that
    ordering with the accepted Palantir pseudotime."""
    set_seed()
    adata = load_annotated()
    attach_embeddings(adata, (MAIN_EMBEDDING,))
    labels = adata.obs[CELL_TYPE_KEY].astype(str).rename(CELL_TYPE_KEY)
    coords = np.asarray(adata.obsm[UMAP_KEY])

    fit_principal_tree(adata)
    order_cells(adata, select_root(adata, MAIN_ROOT_METHODS, MAIN_EMBEDDING))
    nodes, edges = graph_layout(adata)

    raw = adata.obs["t"].astype(float)
    graph_time = ((raw - raw.min()) / (raw.max() - raw.min())).rename("graph")
    palantir_time = pd.read_csv(
        find_result("pt_palantir_pseudotime.csv"), index_col = 0
    )[RUN].reindex(adata.obs_names).rename("palantir")
    segments = adata.obs["seg"].astype(str).rename("segment")

    agreement = agreement_table(graph_time, palantir_time)
    summary = celltype_table(graph_time, palantir_time, labels)
    save_table(
        pd.DataFrame(
            {
                "pseudotime_raw" : raw,
                "pseudotime"     : graph_time,
                "segment"        : segments,
                "milestone"      : adata.obs["milestones"].astype(str),
                "palantir"       : palantir_time,
                CELL_TYPE_KEY    : labels,
            }
        ),
        "principal_graph_pseudotime",
    )
    save_table(nodes, "principal_graph_nodes")
    save_table(
        segment_table(segments, labels, adata.uns["graph"]),
        "principal_graph_segments",
    )
    save_table(summary, "principal_graph_by_celltype")
    save_table(agreement, "principal_graph_agreement", index = False)

    branch = pd.Series(
        pd.Categorical(
            "branch " + segments,
            categories = [
                f"branch {name}" for name in sorted(segments.unique(), key = int)
            ],
        ),
        index = adata.obs_names,
    )
    figures = (
        (labels, "Cell type", "Principal tree over the UMAP",
         "principal_graph_umap"),
        (graph_time, "pseudotime along the tree, 0 at the root",
         "Principal tree pseudotime", "principal_graph_pseudotime"),
        (branch, "Branch", "Principal tree branches",
         "principal_graph_segments"),
    )
    for colour_by, key, title, name in figures:
        tree_figure(coords, nodes, edges, colour_by, key, title, name)
    plot_celltype_pseudotime(graph_time, labels)
    plot_agreement(graph_time, palantir_time, agreement.at[0, "value"]) # type: ignore

    print(summary.to_string())
    print()
    print(agreement.to_string(index = False))


if __name__ == "__main__":
    main()
    