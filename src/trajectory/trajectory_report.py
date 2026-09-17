"""Summary figures and the Palantir agains CellRank comparison for the trajectory stage."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from trajectory_utils import (
    BACKGROUND, CELL_TYPE_KEY, FOREGROUND, MAIN_EMBEDDING, 
    MAIN_ROOT_METHODS, RESULT_DIR, UMAP_LEGEND_FONTSIZE,
    load_annotated, run_tag, save_figure, save_table, set_seed,
    style_bar, style_dark,
)

FATE_TAG        = "k8"
UMAP_KEY        = "X_umap"
RUN             = run_tag(MAIN_EMBEDDING, MAIN_ROOT_METHODS)

POINT_SIZE      = 1.5
POINT_ALPHA     = 0.6
ANNOTATION_SIZE = 6
CATEGORY_CMAP   = "tab20"
PSEUDOTIME_CMAP = "viridis"
ENTROPY_CMAP    = "magma"
FATE_CMAP       = "magma"
AGREEMENT_CMAP  = "RdBu_r"
AGREEMENT_VLIM  = 1.0
OVERVIEW_SIZE   = (18, 13)
HEATMAP_SIZE    = (12, 9)
BOX_SIZE        = (12, 9)


def load_pseudotime():
    """Return the Palantir pseudotime of the primary run."""
    frame = pd.read_csv(RESULT_DIR / "pt_palantir_pseudotime.csv", index_col = 0)
    return frame[RUN]


def load_matrix(name, index):
    """Return a saved per cell matrix as a frame labelled with its feature names."""
    matrix   = np.load(RESULT_DIR / f"{name}.npy")
    features = pd.read_csv(RESULT_DIR / f"{name}_features.csv")["feature"].astype(str)
    return pd.DataFrame(matrix, index=index, columns=features.tolist())


def terminal_labels(columns, labels):
    """Return display names for terminal columns, resolving barcodes to their cell type."""
    return [labels[column] if column in labels.index else column for column in columns]


def normalised_entropy(frame):
    """Return the fate entropy of every cell scaled to the interval zero to one."""
    values = np.clip(frame.to_numpy(dtype = float), 1e-12, None)
    values = values / values.sum(axis = 1, keepdims = True)
    raw    = -(values * np.log(values)).sum(axis = 1)
    return pd.Series(raw / np.log(values.shape[1]), index = frame.index)


def standardise(values):
    """Return the columns of a matrix centred and scaled to unit variance."""
    centred = values - values.mean(axis = 0)
    scale   = centred.std(axis = 0)
    scale[scale == 0] = 1.0
    return centred / scale


def scatter_continuous(figure, axis, coords, values, cmap, title, label):
    """Draw one embedding panel coloured by a continuous value."""
    points = axis.scatter(
        coords[:, 0], coords[:, 1], c = values, s = POINT_SIZE, cmap = cmap,
        alpha = POINT_ALPHA, linewidths = 0, rasterized = True,
    )
    axis.set_title(title, color = FOREGROUND)
    axis.set_xticks([])
    axis.set_yticks([])
    style_dark(figure, axis)
    bar = figure.colorbar(points, ax = axis, fraction = 0.035, pad = 0.02)
    bar.set_label(label)
    style_bar(bar)


def scatter_categorical(figure, axis, coords, labels, title):
    """Draw one embedding panel coloured by a categorical label."""
    categories = sorted(pd.unique(labels))
    colours    = plt.get_cmap(CATEGORY_CMAP)(np.linspace(0.0, 1.0, len(categories)))
    for colour, category in zip(colours, categories):
        mask = labels == category
        axis.scatter(
            coords[mask, 0], coords[mask, 1], s = POINT_SIZE, color = colour,
            alpha = POINT_ALPHA, linewidths = 0, rasterized = True, label = category,
        )
    axis.set_title(title, color = FOREGROUND)
    axis.set_xticks([])
    axis.set_yticks([])
    style_dark(figure, axis)
    legend = axis.legend(
        fontsize=UMAP_LEGEND_FONTSIZE, markerscale = 6, loc = "center left",
        bbox_to_anchor = (1.0, 0.5), frameon = False, labelspacing = 0.3,
    )
    for text in legend.get_texts():
        text.set_color(FOREGROUND)


def plot_overview(coords, labels, pseudotime, fates, entropy):
    """Draw the embedding panels for cell type, pseudotime, assigned fate and fate certainty."""
    figure, axes = plt.subplots(2, 2, figsize = OVERVIEW_SIZE)
    figure.patch.set_facecolor(BACKGROUND)
    scatter_categorical(figure, axes[0, 0], coords, labels.to_numpy(), "Annotated cell type")
    scatter_continuous(
        figure, axes[0, 1], coords, pseudotime.to_numpy(),
        PSEUDOTIME_CMAP, "Palantir pseudotime", "pseudotime",
    )
    scatter_categorical(
        figure, axes[1, 0], coords, fates.idxmax(axis = 1).to_numpy(),
        "Most probable terminal state",
    )
    scatter_continuous(
        figure, axes[1, 1], coords, entropy.to_numpy(),
        ENTROPY_CMAP, "Fate entropy", "normalised entropy",
    )
    figure.tight_layout()
    return figure


def fate_matrix(fates, labels):
    """Return the mean fate probability of every cell type for every terminal state."""
    frame = fates.copy()
    frame[CELL_TYPE_KEY] = labels
    means = frame.groupby(CELL_TYPE_KEY, observed = True).mean()
    key = pd.DataFrame(
        {
            "terminal" : means.to_numpy().argmax(axis = 1),
            "value"    : -means.to_numpy().max(axis = 1),
        },
        index=means.index,
    )
    return means.loc[key.sort_values(["terminal", "value"]).index]


def agreement(palantir, cellrank):
    """Return the correlation of every Palantir branch probability with every CellRank fate."""
    left  = standardise(palantir.to_numpy(dtype = float))
    right = standardise(cellrank.to_numpy(dtype = float))
    return pd.DataFrame(
        left.T @ right / left.shape[0], index = palantir.columns, columns = cellrank.columns
    )


def plot_heatmap(matrix, title, cmap, label, vmin, vmax):
    """Draw an annotated heatmap of a cell type or method comparison matrix."""
    figure, axis = plt.subplots(figsize = HEATMAP_SIZE)
    image = axis.imshow(matrix.to_numpy(), cmap = cmap, vmin = vmin, vmax = vmax, aspect = "auto")
    axis.set_xticks(range(matrix.shape[1]))
    axis.set_xticklabels(matrix.columns, rotation = 45, ha = "right", fontsize = 8)
    axis.set_yticks(range(matrix.shape[0]))
    axis.set_yticklabels(matrix.index, fontsize = 8)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            axis.text(
                column, row, f"{matrix.to_numpy()[row, column]:.2f}",
                ha = "center", va = "center", fontsize = ANNOTATION_SIZE, color = FOREGROUND,
            )
    axis.set_title(title, color = FOREGROUND)
    style_dark(figure, axis)
    for text in axis.get_xticklabels() + axis.get_yticklabels():
        text.set_color(FOREGROUND)
    bar = figure.colorbar(image, ax = axis, fraction = 0.03, pad = 0.02)
    bar.set_label(label)
    style_bar(bar)
    figure.tight_layout()
    return figure


def plot_pseudotime(pseudotime, labels):
    """Draw the pseudotime distribution of every cell type ordered by median."""
    frame  = pd.DataFrame({"pseudotime": pseudotime.to_numpy(), CELL_TYPE_KEY: labels.to_numpy()})
    order  = frame.groupby(CELL_TYPE_KEY, observed=True)["pseudotime"].median().sort_values()
    groups = [
        frame.loc[frame[CELL_TYPE_KEY] == name, "pseudotime"].to_numpy() for name in order.index
    ]
    figure, axis = plt.subplots(figsize=BOX_SIZE)
    parts = axis.boxplot(groups, orientation = "horizontal", showfliers = False, patch_artist = True)
    colours = plt.get_cmap(PSEUDOTIME_CMAP)(np.linspace(0.0, 1.0, len(groups)))
    for patch, colour in zip(parts["boxes"], colours):
        patch.set_facecolor(colour)
        patch.set_edgecolor(FOREGROUND)
    for key in ("whiskers", "caps", "medians"):
        for line in parts[key]:
            line.set_color(FOREGROUND)
    axis.set_yticks(range(1, len(groups) + 1))
    axis.set_yticklabels(list(order.index), fontsize = 8)
    axis.set_xlabel("Palantir pseudotime", color = FOREGROUND)
    axis.set_title("Pseudotime distribution by cell type", color = FOREGROUND)
    style_dark(figure, axis)
    for text in axis.get_yticklabels():
        text.set_color(FOREGROUND)
    figure.tight_layout()
    return figure


def most_common(values):
    """Return the most frequent value of a series."""
    return values.value_counts().idxmax()


def summary_table(pseudotime, fates, entropy, labels):
    """Return per cell type pseudotime, dominant terminal state and fate certainty."""
    frame = pd.DataFrame(
        {
            CELL_TYPE_KEY: labels.to_numpy(),
            "pseudotime": pseudotime.to_numpy(),
            "entropy": entropy.to_numpy(),
            "top_terminal": fates.idxmax(axis = 1).to_numpy(),
            "top_probability": fates.max(axis = 1).to_numpy(),
        }
    )
    grouped = frame.groupby(CELL_TYPE_KEY, observed = True)
    table = pd.DataFrame(
        {
            "n_cells"              : grouped.size(),
            "median_pseudotime"    : grouped["pseudotime"].median(),
            "mean_entropy"         : grouped["entropy"].mean(),
            "mean_top_probability" : grouped["top_probability"].mean(),
            "dominant_terminal"    : grouped["top_terminal"].agg(most_common),
        }
    )
    return table.sort_values("median_pseudotime")


def main():
    """Build the summary figures and tables for the trajectory stage."""
    set_seed()
    adata  = load_annotated()
    labels = adata.obs[CELL_TYPE_KEY].astype(str)
    coords = adata.obsm[UMAP_KEY]

    pseudotime       = load_pseudotime().reindex(adata.obs_names)
    cellrank         = load_matrix(f"pt_cellrank_fate_{FATE_TAG}", adata.obs_names)
    palantir         = load_matrix(f"pt_palantir_fate_{RUN}", adata.obs_names)
    palantir.columns = terminal_labels(palantir.columns, labels)
    entropy          = normalised_entropy(cellrank)

    save_figure(
        plot_overview(coords, labels, pseudotime, cellrank, entropy),
        "trajectory_report_overview",
    )

    matrix = fate_matrix(cellrank, labels)
    save_table(matrix, "trajectory_report_fate_matrix")
    save_figure(
        plot_heatmap(
            matrix,
            f"Mean CellRank fate probability by cell type ({FATE_TAG})",
            FATE_CMAP, "mean probability", 0.0, 1.0,
        ),
        "trajectory_report_fate_matrix",
    )

    save_figure(plot_pseudotime(pseudotime, labels), "trajectory_report_pseudotime_by_celltype")

    correlation = agreement(palantir, cellrank)
    save_table(correlation, "trajectory_report_method_agreement")
    save_figure(
        plot_heatmap(
            correlation,
            "Palantir branch probability against CellRank fate probability",
            AGREEMENT_CMAP, "Pearson correlation", -AGREEMENT_VLIM, AGREEMENT_VLIM,
        ),
        "trajectory_report_method_agreement",
    )

    summary = summary_table(pseudotime, cellrank, entropy, labels)
    save_table(summary, "trajectory_report_summary")
    print(summary.to_string())
    print()
    print(correlation.round(2).to_string())


if __name__ == "__main__":
    main()