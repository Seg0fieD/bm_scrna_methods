"""Backbone views of the trajectory: the PAGA spanning tree over the UMAP, the
fate entropy distribution, and the transition-matrix arrow field."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib.axes import Axes

from pt_cellrank import BASE_EMBEDDING, build_kernel
from trajectory_utils import (
    BACKGROUND, CELL_TYPE_KEY, FOREGROUND, MAIN_EMBEDDING, MAIN_ROOT_METHODS,
    UMAP_LEGEND_FONTSIZE, attach_embeddings, find_result, load_annotated,
    run_tag, save_figure, save_table, set_seed, style_bar, style_dark,
)
SCRIPT_NAME     = "trajectory_backbone"
RUN             = run_tag(MAIN_EMBEDDING, MAIN_ROOT_METHODS)
UMAP_KEY        = "X_umap"

POINT_SIZE      = 1.5
POINT_ALPHA     = 0.55
BACKDROP_COLOUR = "#39404f"
EDGE_COLOUR     = "#9aa4b8"
EDGE_WIDTH      = (0.8, 4.5)
NODE_SIZE       = 320
LABEL_SIZE      = 9
CATEGORY_CMAP   = "tab20"
PSEUDOTIME_CMAP = "viridis"
ENTROPY_CMAP    = "magma"

GRID_SIZE       = 42
MIN_BIN_CELLS   = 5
ARROW_LENGTH    = 1.7
ARROW_FLOOR     = 0.08
ARROW_WIDTH     = 0.0020
ARROW_COLOUR    = "#f2c14e"

TREE_FIGSIZE    = (14, 11)
BOX_FIGSIZE     = (12, 9)
CAPTION_SIZE    = 10

ARROW_CAPTION = (
    "Arrows are the CellRank transition matrix projected onto the UMAP. They "
    "show the direction the chosen root cell imposes on the ordering. This is "
    "not RNA velocity and carries no evidence independent of that root."
)


def load_series(name: str, column: str) -> pd.Series:
    """One column of a saved per-cell table, indexed by barcode."""
    frame = pd.read_csv(find_result(f"{name}.csv"), index_col = 0)
    if column not in frame.columns:
        raise KeyError(f"{name}.csv has no column {column}")
    return frame[column]


def load_tree() -> pd.DataFrame:
    """PAGA spanning tree of the primary embedding as a cell-type matrix."""
    path = find_result(f"paga_topology_{MAIN_EMBEDDING}_tree.csv")
    return pd.read_csv(path, index_col = 0)


def centroid_frame(
    adata: AnnData, pseudotime: pd.Series, entropy: pd.Series
) -> pd.DataFrame:
    """Per cell type: the median UMAP position, cell count, median pseudotime
    and mean fate entropy."""
    coords = adata.obsm[UMAP_KEY]
    frame = pd.DataFrame(
        {
            CELL_TYPE_KEY : adata.obs[CELL_TYPE_KEY].astype(str).to_numpy(),
            "umap_1"      : coords[:, 0],
            "umap_2"      : coords[:, 1],
            "pseudotime"  : pseudotime.to_numpy(),
            "entropy"     : entropy.to_numpy(),
        }
    )
    grouped = frame.groupby(CELL_TYPE_KEY, observed = True)
    table = pd.DataFrame(
        {
            "cells"             : grouped.size(),
            "umap_1"            : grouped["umap_1"].median(),
            "umap_2"            : grouped["umap_2"].median(),
            "median_pseudotime" : grouped["pseudotime"].median(),
            "mean_entropy"      : grouped["entropy"].mean(),
        }
    )
    return table.sort_values("median_pseudotime")


def tree_edges(tree: pd.DataFrame, centroids: pd.DataFrame) -> pd.DataFrame:
    """Undirected edges of the spanning tree, oriented from the earlier to the
    later cell type and carrying both endpoint positions."""
    matrix = tree.to_numpy(dtype = float)
    rows, cols = np.nonzero(matrix)
    seen: set[tuple[str, str]] = set()
    records = []
    for row, col in zip(rows, cols):
        left, right = str(tree.index[row]), str(tree.columns[col])
        if left == right:
            continue
        key = (left, right) if left < right else (right, left)
        if key in seen:
            continue
        seen.add(key)
        times = centroids["median_pseudotime"]
        source, target = (left, right)
        if times.get(source, 0.0) > times.get(target, 0.0):
            source, target = target, source
        records.append(
            {
                "source"       : source,
                "target"       : target,
                "connectivity" : float(matrix[row, col]),
                "source_x"     : centroids.at[source, "umap_1"],
                "source_y"     : centroids.at[source, "umap_2"],
                "target_x"     : centroids.at[target, "umap_1"],
                "target_y"     : centroids.at[target, "umap_2"],
            }
        )
    table = pd.DataFrame.from_records(records)
    return table.sort_values(
        "connectivity", ascending = False, ignore_index = True
    )


def edge_widths(values: np.ndarray) -> np.ndarray:
    """Line widths scaled linearly between the fixed width limits."""
    low, high = EDGE_WIDTH
    span = values.max() - values.min()
    if span <= 0:
        return np.full(values.shape, 0.5 * (low + high))
    return low + (high - low) * (values - values.min()) / span


def draw_backbone(axis: Axes, edges: pd.DataFrame, labels: bool = True) -> None:
    """Spanning-tree segments between cell-type centroids, with the cell-type
    names written at the nodes."""
    widths = edge_widths(edges["connectivity"].to_numpy(dtype = float))
    for width, row in zip(widths, edges.itertuples(index = False)):
        axis.plot(
            [row.source_x, row.target_x], # pyright: ignore[reportArgumentType]
            [row.source_y, row.target_y], # pyright: ignore[reportArgumentType]
            color     = EDGE_COLOUR,
            linewidth = width,
            alpha     = 0.9,
            zorder    = 3,
            solid_capstyle = "round",
        )
    if not labels:
        return
    placed: set[str] = set()
    for row in edges.itertuples(index = False):
        for name, x_position, y_position in (
            (row.source, row.source_x, row.source_y),
            (row.target, row.target_x, row.target_y),
        ):
            if name in placed:
                continue
            placed.add(name)
            axis.text(
                x_position,
                y_position,
                name,
                fontsize = LABEL_SIZE,
                color    = BACKGROUND,
                ha       = "center",
                va       = "center",
                zorder   = 6,
                bbox     = {
                    "boxstyle"  : "round,pad=0.25",
                    "facecolor" : FOREGROUND,
                    "edgecolor" : "none",
                    "alpha"     : 0.85,
                },
            )


def blank_axes(axis: Axes) -> None:
    """Axes stripped of ticks and frame, for an embedding panel."""
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)


def plot_tree_cell_type(
    coords: np.ndarray, labels: pd.Series, edges: pd.DataFrame
) -> None:
    """Spanning tree over the UMAP with the cells coloured by annotated cell
    type."""
    figure, axis = plt.subplots(figsize = TREE_FIGSIZE)
    categories = sorted(pd.unique(labels))
    colours = plt.get_cmap(CATEGORY_CMAP)(np.linspace(0.0, 1.0, len(categories)))
    for colour, category in zip(colours, categories):
        mask = (labels == category).to_numpy()
        axis.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s          = POINT_SIZE,
            color      = colour,
            alpha      = POINT_ALPHA,
            linewidths = 0,
            rasterized = True,
            label      = category,
            zorder     = 1,
        )
    draw_backbone(axis, edges)
    axis.set_title("Trajectory backbone over the UMAP", color = FOREGROUND)
    style_dark(figure, axis)
    blank_axes(axis)
    legend = axis.legend(
        fontsize       = UMAP_LEGEND_FONTSIZE,
        markerscale    = 6,
        loc            = "center left",
        bbox_to_anchor = (1.0, 0.5),
        frameon        = False,
        labelspacing   = 0.3,
    )
    for text in legend.get_texts():
        text.set_color(FOREGROUND)
    figure.tight_layout()
    save_figure(figure, "trajectory_backbone_tree_cell_type")


def plot_tree_pseudotime(
    coords: np.ndarray, pseudotime: pd.Series, edges: pd.DataFrame
) -> None:
    """Spanning tree over the UMAP with the cells coloured by pseudotime."""
    figure, axis = plt.subplots(figsize = TREE_FIGSIZE)
    points = axis.scatter(
        coords[:, 0],
        coords[:, 1],
        c          = pseudotime.to_numpy(),
        cmap       = PSEUDOTIME_CMAP,
        s          = POINT_SIZE,
        alpha      = POINT_ALPHA,
        linewidths = 0,
        rasterized = True,
        vmin       = 0.0,
        vmax       = 1.0,
        zorder     = 1,
    )
    draw_backbone(axis, edges)
    axis.set_title("Trajectory backbone against pseudotime", color = FOREGROUND)
    style_dark(figure, axis)
    blank_axes(axis)
    bar = figure.colorbar(points, ax = axis, fraction = 0.035, pad = 0.02)
    bar.set_label("Palantir pseudotime", color = FOREGROUND)
    style_bar(bar)
    figure.tight_layout()
    save_figure(figure, "trajectory_backbone_tree_pseudotime")


def plot_tree_nodes(
    coords: np.ndarray, centroids: pd.DataFrame, edges: pd.DataFrame
) -> None:
    """Spanning tree with the nodes sized by cell count and coloured by median
    pseudotime, over a neutral cell backdrop."""
    figure, axis = plt.subplots(figsize = TREE_FIGSIZE)
    axis.scatter(
        coords[:, 0],
        coords[:, 1],
        s          = POINT_SIZE,
        color      = BACKDROP_COLOUR,
        alpha      = 0.45,
        linewidths = 0,
        rasterized = True,
        zorder     = 1,
    )
    draw_backbone(axis, edges, labels = False)
    sizes = NODE_SIZE * np.sqrt(
        centroids["cells"].to_numpy(dtype = float) / centroids["cells"].max()
    )
    nodes = axis.scatter(
        centroids["umap_1"].to_numpy(),
        centroids["umap_2"].to_numpy(),
        c          = centroids["median_pseudotime"].to_numpy(),
        cmap       = PSEUDOTIME_CMAP,
        s          = sizes,
        vmin       = 0.0,
        vmax       = 1.0,
        edgecolors = FOREGROUND,
        linewidths = 0.6,
        zorder     = 5,
    )
    for name, row in centroids.iterrows():
        axis.annotate(
            str(name),
            (row["umap_1"], row["umap_2"]),
            textcoords = "offset points",
            xytext     = (0, 11),
            fontsize   = LABEL_SIZE,
            color      = FOREGROUND,
            ha         = "center",
            zorder     = 6,
        )
    axis.set_title(
        "Cell-type backbone, nodes by median pseudotime", color = FOREGROUND
    )
    style_dark(figure, axis)
    blank_axes(axis)
    bar = figure.colorbar(nodes, ax = axis, fraction = 0.035, pad = 0.02)
    bar.set_label("median pseudotime", color = FOREGROUND)
    style_bar(bar)
    figure.tight_layout()
    save_figure(figure, "trajectory_backbone_tree_nodes")


def plot_entropy(entropy: pd.Series, labels: pd.Series) -> None:
    """Fate entropy distribution of every cell type, ordered by median."""
    frame = pd.DataFrame(
        {"entropy": entropy.to_numpy(), CELL_TYPE_KEY: labels.to_numpy()}
    )
    order = (
        frame.groupby(CELL_TYPE_KEY, observed = True)["entropy"]
        .median()
        .sort_values()
    )
    groups = [
        frame.loc[frame[CELL_TYPE_KEY] == name, "entropy"].to_numpy()
        for name in order.index
    ]
    figure, axis = plt.subplots(figsize = BOX_FIGSIZE)
    parts = axis.boxplot(
        groups,
        orientation  = "horizontal",
        showfliers   = False,
        patch_artist = True,
    )
    colours = plt.get_cmap(ENTROPY_CMAP)(np.linspace(0.15, 0.95, len(groups)))
    for patch, colour in zip(parts["boxes"], colours):
        patch.set_facecolor(colour)
        patch.set_edgecolor(FOREGROUND)
    for key in ("whiskers", "caps", "medians"):
        for line in parts[key]:
            line.set_color(FOREGROUND)
    axis.set_yticks(range(1, len(groups) + 1))
    axis.set_yticklabels(list(order.index), fontsize = 8)
    axis.set_xlabel("Palantir fate entropy", color = FOREGROUND)
    axis.set_title("Fate entropy by cell type", color = FOREGROUND)
    style_dark(figure, axis)
    for text in axis.get_yticklabels():
        text.set_color(FOREGROUND)
    figure.tight_layout()
    save_figure(figure, "trajectory_backbone_entropy")


def project_transitions(matrix, coords: np.ndarray) -> np.ndarray:
    """Per-cell displacement in embedding space: transition probabilities over
    unit neighbour directions, less the isotropic expectation of that
    neighbourhood."""
    sparse = matrix.tocoo()
    rows, cols, values = sparse.row, sparse.col, sparse.data.astype(float)
    offsets = coords[cols] - coords[rows]
    lengths = np.linalg.norm(offsets, axis = 1)
    keep = lengths > 0
    rows, cols, values = rows[keep], cols[keep], values[keep]
    directions = offsets[keep] / lengths[keep, None]

    degree = np.bincount(rows, minlength = coords.shape[0]).astype(float)
    degree[degree == 0] = 1.0
    weights = values - 1.0 / degree[rows]

    field = np.zeros_like(coords, dtype = float)
    np.add.at(field, rows, weights[:, None] * directions)
    return field


def grid_field(
    coords: np.ndarray, field: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mean displacement per occupied grid cell: the grid positions and the two
    displacement components, bins below the occupancy floor dropped."""
    edges = [
        np.linspace(coords[:, axis].min(), coords[:, axis].max(), GRID_SIZE + 1)
        for axis in (0, 1)
    ]
    index_x = np.clip(np.digitize(coords[:, 0], edges[0]) - 1, 0, GRID_SIZE - 1)
    index_y = np.clip(np.digitize(coords[:, 1], edges[1]) - 1, 0, GRID_SIZE - 1)
    flat = index_y * GRID_SIZE + index_x

    counts = np.bincount(flat, minlength = GRID_SIZE ** 2).astype(float)
    sum_x = np.bincount(flat, weights = field[:, 0], minlength = GRID_SIZE ** 2)
    sum_y = np.bincount(flat, weights = field[:, 1], minlength = GRID_SIZE ** 2)

    occupied = counts >= MIN_BIN_CELLS
    centres_x = 0.5 * (edges[0][:-1] + edges[0][1:])
    centres_y = 0.5 * (edges[1][:-1] + edges[1][1:])
    mesh_x, mesh_y = np.meshgrid(centres_x, centres_y)

    shift_x = sum_x[occupied] / counts[occupied]
    shift_y = sum_y[occupied] / counts[occupied]
    spacing = float(min(np.diff(edges[0]).mean(), np.diff(edges[1]).mean()))
    shift_x, shift_y = scale_arrows(shift_x, shift_y, spacing)

    keep = np.hypot(shift_x, shift_y) >= ARROW_FLOOR * ARROW_LENGTH * spacing
    return (
        mesh_x.ravel()[occupied][keep],
        mesh_y.ravel()[occupied][keep],
        shift_x[keep],
        shift_y[keep],
    )


def scale_arrows(
    shift_x: np.ndarray, shift_y: np.ndarray, spacing: float
) -> tuple[np.ndarray, np.ndarray]:
    """Displacements rescaled so the longest typical arrow spans a fixed
    multiple of the grid spacing, leaving the figure independent of the raw
    displacement magnitude."""
    lengths = np.hypot(shift_x, shift_y)
    reference = float(np.percentile(lengths, 95)) if lengths.size else 0.0
    if reference <= 0.0:
        return shift_x, shift_y
    factor = ARROW_LENGTH * spacing / reference
    return shift_x * factor, shift_y * factor


def plot_arrows(
    coords: np.ndarray, pseudotime: pd.Series, field: np.ndarray
) -> None:
    """Grid-averaged transition field over the UMAP, on a pseudotime
    backdrop."""
    grid_x, grid_y, shift_x, shift_y = grid_field(coords, field)
    figure, axis = plt.subplots(figsize = TREE_FIGSIZE)
    points = axis.scatter(
        coords[:, 0],
        coords[:, 1],
        c          = pseudotime.to_numpy(),
        cmap       = PSEUDOTIME_CMAP,
        s          = POINT_SIZE,
        alpha      = 0.35,
        linewidths = 0,
        rasterized = True,
        vmin       = 0.0,
        vmax       = 1.0,
        zorder     = 1,
    )
    axis.quiver(
        grid_x,
        grid_y,
        shift_x,
        shift_y,
        color       = ARROW_COLOUR,
        angles      = "xy",
        scale_units = "xy",
        scale       = 1.0,
        width       = ARROW_WIDTH,
        zorder      = 4,
    )
    axis.set_title("Transition direction on the UMAP", color = FOREGROUND)
    style_dark(figure, axis)
    blank_axes(axis)
    bar = figure.colorbar(points, ax = axis, fraction = 0.035, pad = 0.02)
    bar.set_label("Palantir pseudotime", color = FOREGROUND)
    style_bar(bar)
    figure.text(
        0.5,
        0.012,
        ARROW_CAPTION,
        ha       = "center",
        fontsize = CAPTION_SIZE,
        color    = FOREGROUND,
        wrap     = True,
    )
    figure.tight_layout(rect = (0.0, 0.05, 1.0, 1.0))
    save_figure(figure, "trajectory_backbone_arrows")


def main() -> None:
    """Backbone figures and the centroid and edge tables of the primary run."""
    set_seed()
    adata = load_annotated()
    attach_embeddings(adata, (BASE_EMBEDDING,))

    labels = adata.obs[CELL_TYPE_KEY].astype(str)
    coords = adata.obsm[UMAP_KEY]
    order = adata.obs_names
    pseudotime = load_series("pt_palantir_pseudotime", RUN).reindex(order)
    entropy = load_series("pt_palantir_entropy", RUN).reindex(order)

    centroids = centroid_frame(adata, pseudotime, entropy)
    edges = tree_edges(load_tree(), centroids)
    save_table(centroids, SCRIPT_NAME, "trajectory_backbone_centroids")
    save_table(edges, SCRIPT_NAME, "trajectory_backbone_tree_edges", index = False)

    plot_tree_cell_type(coords, labels, edges)
    plot_tree_pseudotime(coords, pseudotime, edges)
    plot_tree_nodes(coords, centroids, edges)
    plot_entropy(entropy, labels)

    time_key = f"pseudotime_{RUN}"
    adata.obs[time_key] = pseudotime.to_numpy()
    kernel = build_kernel(adata, time_key)
    field = project_transitions(kernel.transition_matrix, coords)
    plot_arrows(coords, pseudotime, field)

    print(centroids.to_string())


if __name__ == "__main__":
    main()