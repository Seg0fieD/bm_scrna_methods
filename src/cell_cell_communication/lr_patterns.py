"""Signalling programs from factorisation of the interaction matrix."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path
from sklearn.decomposition import NMF

import cell_communication_utils as ccu

SCRIPT_NAME      = "lr_patterns"
RESOURCE         = "consensus"
DIRECTIONS       = ("outgoing", "incoming")
DIRECTION_KEYS   = {"outgoing": "source", "incoming": "target"}
RANK_RANGE       = (2, 3, 4, 5, 6)
N_PATTERNS       = 4
TOP_LOADINGS     = 4
MAX_ITER         = 2000
NODE_GAP         = 0.012
LABEL_FONTSIZE   = 10
FIGURE_DIRECTION = "outgoing"


def load_significant(resource: str) -> pd.DataFrame:
    """
        Significant interactions for one resource: sender, receiver and the
        ligand-receptor pair behind each one.
    """
    path = ccu.RESULTS_DIR / "lr_rank" / f"lr_rank_{resource}_significant.csv"
    return pd.read_csv(path)


def interaction_matrix(frame: pd.DataFrame, direction: str) -> pd.DataFrame:
    """
        Cell-type-by-pair matrix of significant interaction counts, taken over
        senders for the outgoing direction and receivers for the incoming one.
    """
    key = DIRECTION_KEYS[direction]
    pair = frame["ligand_complex"] + " -> " + frame["receptor_complex"]
    return pd.crosstab(frame[key], pair)


def scaled_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    """
        Interaction matrix with every pair divided by its own maximum, so no
        single ligand dominates the factorisation by size alone.
    """
    return counts / counts.max(axis = 0)


def normalise_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Matrix with every row rescaled to sum to one."""
    totals = frame.sum(axis = 1).replace(0, 1)
    return frame.div(totals, axis = 0)


def rank_scan(matrix: pd.DataFrame) -> pd.DataFrame:
    """
        Reconstruction error and iteration count at each candidate number of
        patterns.
    """
    records = []
    for rank in RANK_RANGE:
        model = NMF(
            n_components = rank,
            init         = "nndsvda",
            random_state = ccu.SEED,
            max_iter     = MAX_ITER,
        )
        model.fit(matrix.to_numpy())
        records.append(
            {
                "n_patterns"          : rank,
                "reconstruction_error": float(model.reconstruction_err_),
                "n_iter"              : int(model.n_iter_),
            }
        )
    return pd.DataFrame(records)


def factorise(matrix: pd.DataFrame, n_patterns: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
        Cell-type-by-pattern usage and pattern-by-pair loading matrices from one
        non-negative factorisation.
    """
    model = NMF(
        n_components = n_patterns,
        init         = "nndsvda",
        random_state = ccu.SEED,
        max_iter     = MAX_ITER,
    )
    weights = model.fit_transform(matrix.to_numpy())
    names = [f"pattern {position + 1}" for position in range(n_patterns)]

    usage = pd.DataFrame(weights, index = matrix.index, columns = names)
    usage.index.name = "cell_type"

    loadings = pd.DataFrame(model.components_, index = names, columns = matrix.columns)
    loadings.index.name = "pattern"

    return usage, loadings


def usage_table(usage: pd.DataFrame, totals: pd.Series, direction: str) -> pd.DataFrame:
    """
        Long-form usage: one row per cell type and pattern, with that pattern's
        share of the cell type and the cell type's total weight before
        rescaling.
    """
    tidy = usage.reset_index().melt(
        id_vars = "cell_type", var_name = "pattern", value_name = "usage"
    )
    tidy["total_usage"] = tidy["cell_type"].map(totals)
    tidy.insert(0, "direction", direction)
    return tidy


def loading_table(loadings: pd.DataFrame, direction: str) -> pd.DataFrame:
    """
        Long-form loadings: one row per pattern and ligand-receptor pair,
        carrying that pair's share of the pattern.
    """
    tidy = loadings.reset_index().melt(
        id_vars = "pattern", var_name = "interaction", value_name = "loading"
    )
    tidy.insert(0, "direction", direction)
    return tidy


def pattern_colour(position: int):
    """Colour for one pattern, taken from the qualitative palette by position."""
    return plt.get_cmap("tab10")(position % 10)


def column_scale(columns: list[np.ndarray], gap: float) -> float:
    """Shared vertical scale making the tallest stacked column fill the frame."""
    scales = []
    for heights in columns:
        usable = 1.0 - gap * max(len(heights) - 1, 0)
        scales.append(usable / float(heights.sum()))
    return min(scales)


def stacked_positions(
    heights: np.ndarray, gap: float, scale: float
) -> tuple[np.ndarray, np.ndarray]:
    """
        Top coordinate and drawn height of every node in one column, stacked
        with a fixed gap at a shared scale and centred in the frame.
    """
    drawn = heights * scale
    tops = np.zeros(len(heights))
    cursor = (1.0 - (float(drawn.sum()) + gap * max(len(heights) - 1, 0))) / 2
    for position, height in enumerate(drawn):
        tops[position] = cursor
        cursor += height + gap
    return tops, drawn


def ribbon(ax, span: tuple[float, float], tops: tuple[float, float], height: float, colour) -> None:
    """Curved band joining two stacked nodes, drawn as one filled path."""
    left, right = span
    top_left, top_right = tops
    middle = (left + right) / 2
    vertices = [
        (left, top_left),
        (middle, top_left),
        (middle, top_right),
        (right, top_right),
        (right, top_right + height),
        (middle, top_right + height),
        (middle, top_left + height),
        (left, top_left + height),
        (left, top_left),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    ax.add_patch(
        PathPatch(Path(vertices, codes), facecolor = colour, edgecolor = "none", alpha = 0.55)
    )


def leading_interactions(loadings: pd.DataFrame) -> tuple[list[str], dict[str, list[str]]]:
    """
        Highest-loading pairs of every pattern: the pooled ordered list, and the
        per-pattern lists behind it.
    """
    per_pattern: dict[str, list[str]] = {}
    pooled: list[str] = []
    for pattern in loadings.index:
        ordered = loadings.loc[pattern].sort_values(ascending = False)
        chosen = list(ordered.index[:TOP_LOADINGS])
        per_pattern[pattern] = chosen
        for name in chosen:
            if name not in pooled:
                pooled.append(name)
    return pooled, per_pattern


def plot_programs(usage: pd.DataFrame, loadings: pd.DataFrame) -> Figure:
    """
        River plot of one direction: cell types on the left, shared programs in
        the middle, and the interactions those programs use on the right.
    """
    patterns = list(usage.columns)
    left = usage.to_numpy()
    pattern_height = left.sum(axis = 0)

    interactions, per_pattern = leading_interactions(loadings)
    right = np.zeros((len(patterns), len(interactions)))
    for position, pattern in enumerate(patterns):
        share = loadings.loc[pattern, per_pattern[pattern]]
        share = share / float(share.sum())
        for name, value in share.items():
            right[position, interactions.index(name)] = value * pattern_height[position]

    cell_heights = left.sum(axis = 1)
    interaction_heights = right.sum(axis = 0)
    scale = column_scale([cell_heights, pattern_height, interaction_heights], NODE_GAP)

    cell_top, cell_size = stacked_positions(cell_heights, NODE_GAP, scale)
    pattern_top, pattern_size = stacked_positions(pattern_height, NODE_GAP, scale)
    interaction_top, interaction_size = stacked_positions(interaction_heights, NODE_GAP, scale)

    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)

    cell_cursor = cell_top.copy()
    pattern_cursor = pattern_top.copy()
    for row in range(len(cell_heights)):
        for position in range(len(patterns)):
            value = float(left[row, position]) * scale
            if value <= 0:
                continue
            ribbon(
                ax, (0.215, 0.455),
                (cell_cursor[row], pattern_cursor[position]),
                value, pattern_colour(position),
            )
            cell_cursor[row] += value
            pattern_cursor[position] += value

    pattern_cursor = pattern_top.copy()
    interaction_cursor = interaction_top.copy()
    for position in range(len(patterns)):
        for column in range(len(interactions)):
            value = float(right[position, column]) * scale
            if value <= 0:
                continue
            ribbon(
                ax, (0.545, 0.785),
                (pattern_cursor[position], interaction_cursor[column]),
                value, pattern_colour(position),
            )
            pattern_cursor[position] += value
            interaction_cursor[column] += value

    for row, name in enumerate(usage.index):
        ax.add_patch(
            Rectangle((0.200, cell_top[row]), 0.015, cell_size[row], color = ccu.DARK_FOREGROUND)
        )
        ax.text(
            0.192, cell_top[row] + cell_size[row] / 2, name,
            ha = "right", va = "center",
            fontsize = LABEL_FONTSIZE, color = ccu.DARK_FOREGROUND,
        )

    for position, name in enumerate(patterns):
        ax.add_patch(
            Rectangle(
                (0.455, pattern_top[position]), 0.090, pattern_size[position],
                color = pattern_colour(position),
            )
        )
        ax.text(
            0.500, pattern_top[position] + pattern_size[position] / 2, name,
            ha = "center", va = "center",
            fontsize = LABEL_FONTSIZE, color = ccu.DARK_BACKGROUND,
        )

    for column, name in enumerate(interactions):
        ax.add_patch(
            Rectangle(
                (0.785, interaction_top[column]), 0.015, interaction_size[column],
                color = ccu.DARK_FOREGROUND,
            )
        )
        ax.text(
            0.808, interaction_top[column] + interaction_size[column] / 2, name,
            ha = "left", va = "center",
            fontsize = LABEL_FONTSIZE, color = ccu.DARK_FOREGROUND,
        )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(1.02, -0.02)
    ax.axis("off")
    ax.set_title(
        f"Shared {FIGURE_DIRECTION} signalling programs ({RESOURCE})",
        color = ccu.DARK_FOREGROUND,
    )
    ccu.style_dark(fig)
    fig.tight_layout()
    return fig


def main() -> None:
    """
        Factorisation of the outgoing and incoming interaction matrices, the
        usage, loading and rank tables, and the river plot of one direction.
    """
    ccu.set_seed()
    significant = load_significant(RESOURCE)

    usage_parts = []
    loading_parts = []
    rank_parts = []
    figure_usage = pd.DataFrame()
    figure_loadings = pd.DataFrame()

    for direction in DIRECTIONS:
        counts = interaction_matrix(significant, direction)
        matrix = scaled_matrix(counts)

        scan = rank_scan(matrix)
        scan.insert(0, "direction", direction)
        rank_parts.append(scan)

        usage, loadings = factorise(matrix, N_PATTERNS)
        totals = usage.sum(axis = 1)
        usage = normalise_rows(usage)
        loadings = normalise_rows(loadings)

        usage_parts.append(usage_table(usage, totals, direction))
        loading_parts.append(loading_table(loadings, direction))

        print(
            f"{direction}: {counts.shape[0]} cell types x {counts.shape[1]} interactions, "
            f"{N_PATTERNS} patterns"
        )

        if direction == FIGURE_DIRECTION:
            figure_usage = usage
            figure_loadings = loadings

    ccu.save_table(
        pd.concat(usage_parts, ignore_index = True), SCRIPT_NAME, f"{SCRIPT_NAME}_usage"
    )
    ccu.save_table(
        pd.concat(loading_parts, ignore_index = True), SCRIPT_NAME, f"{SCRIPT_NAME}_loadings"
    )
    ccu.save_table(
        pd.concat(rank_parts, ignore_index = True), SCRIPT_NAME, f"{SCRIPT_NAME}_rank"
    )

    ccu.save_figure(plot_programs(figure_usage, figure_loadings), f"{SCRIPT_NAME}_programs")


if __name__ == "__main__":
    main()