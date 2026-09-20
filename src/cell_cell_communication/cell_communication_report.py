"""Cross-cutting figures and summary tables of the ranked interactions."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, PathPatch, Wedge
from matplotlib.path import Path

import cell_communication_utils as ccu

SCRIPT_NAME      = "cell_communication_report"
MAIN_RESOURCE    = "consensus"
TOP_INTERACTIONS = 30
RANK_FLOOR       = 1e-12
EDGE_WIDTH       = (0.4, 6.0)
DOT_SIZE         = (20.0, 220.0)
NODE_GAP         = 0.035
ARC_WIDTH        = 0.06

LABEL_FONTSIZE   = 10
TITLE_FONTSIZE   = 14
LEGEND_FONTSIZE  = 12

CIRCLE_XLIM      = (-2.6, 2.6)
CIRCLE_YLIM      = (-1.6, 2.2)
TITLE_Y          = 2.05

SENT_COLOUR      = "#e2504a"
RECEIVED_COLOUR  = "#e6e6e6"

DEVELOPING = (
    "HSPC",
    "Granulocyte precursor",
    "Erythroid early",
    "Erythroid late",
    "Pre-B large",
    "Pre-B small",
)

FINISHED = (
    "Monocyte classical",
    "Monocyte non-classical",
    "Dendritic cell",
    "Plasmacytoid DC",
    "Platelet",
    "B naive",
    "CD4 T naive",
    "CD4 T memory",
    "T activated",
    "CD8 T memory",
    "CD8 T cytotoxic",
    "NK",
    "T proliferating",
)

COMPARTMENT_COLOURS = {"developing": "#f2b134", "finished": "#5ab1bb"}


def load_all(resource: str) -> pd.DataFrame:
    """
        Full interaction table for one resource: every tested ligand-receptor
        pair with its per-method scores and merged ranks.
    """
    path = ccu.RESULTS_DIR / "lr_rank" / f"lr_rank_{resource}_all.csv"
    return pd.read_csv(path)


def load_significant(resource: str) -> pd.DataFrame:
    """
        Significant interactions for one resource: sender, receiver and the
        ligand-receptor pair behind each one.
    """
    path = ccu.RESULTS_DIR / "lr_rank" / f"lr_rank_{resource}_significant.csv"
    return pd.read_csv(path)


def ordered_labels(frame: pd.DataFrame) -> list[str]:
    """
        Cell type labels in the fixed developing-then-finished order, with any
        label not on either list appended alphabetically.
    """
    observed = set(frame["source"]).union(frame["target"])
    ordered = [name for name in DEVELOPING + FINISHED if name in observed]
    ordered += sorted(observed.difference(ordered))
    return ordered


def compartment(label: str) -> str:
    """Developmental compartment a cell type belongs to."""
    if label in DEVELOPING:
        return "developing"
    return "finished"


def counts_matrix(significant: pd.DataFrame, labels: list[str]) -> pd.DataFrame:
    """
        Sender-by-receiver matrix of significant interaction counts, reindexed
        to a fixed label list so runs share their axes.
    """
    counts = pd.crosstab(significant["source"], significant["target"])
    return counts.reindex(index = labels, columns = labels, fill_value = 0)


def interaction_pairs(frame: pd.DataFrame) -> pd.Series:
    """Ligand-receptor pair of every row as one label."""
    return frame["ligand_complex"] + " -> " + frame["receptor_complex"]


def rank_score(values) -> np.ndarray:
    """
        Aggregate rank as a positive score: negative log10, floored so zero
        stays finite.
    """
    return -np.log10(np.clip(np.asarray(values, dtype = float), RANK_FLOOR, 1.0))


def scale_between(values: np.ndarray, low: float, high: float) -> np.ndarray:
    """Values rescaled onto a fixed range, flat input mapped to the midpoint."""
    span = float(values.max()) - float(values.min())
    if span <= 0:
        return np.full(len(values), (low + high) / 2)
    return low + (values - float(values.min())) * (high - low) / span


def edge_width(value: float, low: float, high: float) -> float:
    """Line width of one edge, scaled between the fixed bounds by count."""
    if high <= low:
        return (EDGE_WIDTH[0] + EDGE_WIDTH[1]) / 2
    return EDGE_WIDTH[0] + (value - low) * (EDGE_WIDTH[1] - EDGE_WIDTH[0]) / (high - low)


def label_colour(position: int):
    """Colour for one cell type, taken from the qualitative palette by position."""
    return plt.get_cmap("tab20")(position % 20)


def outward_text(ax, radius: float, angle: float, name: str) -> None:
    """Label placed outside the circle, reading outward from its node."""
    facing_right = np.cos(angle) >= 0
    ax.text(
        radius * np.cos(angle), radius * np.sin(angle), name,
        ha = "left" if facing_right else "right",
        va = "center",
        rotation = np.degrees(angle) if facing_right else np.degrees(angle) - 180,
        rotation_mode = "anchor",
        fontsize = LABEL_FONTSIZE, color = ccu.DARK_FOREGROUND,
    )


def circle_title(ax, text: str) -> None:
    """Figure title placed in the blank band reserved above the circle."""
    ax.text(
        0.0, TITLE_Y, text,
        ha = "center", va = "center",
        fontsize = TITLE_FONTSIZE, color = ccu.DARK_FOREGROUND,
    )


def summary_table(counts: pd.DataFrame) -> pd.DataFrame:
    """
        Per-cell-type totals: interactions sent and received, self-signalling,
        and the split of what each cell type sends across the two compartments.
    """
    developing = [name for name in counts.columns if compartment(name) == "developing"]
    finished = [name for name in counts.columns if compartment(name) == "finished"]

    records = []
    for name in counts.index:
        records.append(
            {
                "cell_type"         : name,
                "compartment"       : compartment(name),
                "sent"              : int(counts.loc[name].sum()),
                "received"          : int(counts[name].sum()),
                "self"              : int(counts.loc[name, name]), # type: ignore
                "sent_to_developing": int(counts.loc[name, developing].sum()),
                "sent_to_finished"  : int(counts.loc[name, finished].sum()),
            }
        )
    frame = pd.DataFrame(records)
    return frame.sort_values("sent", ascending = False, ignore_index = True)


def overlap_table(significant: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
        Interaction-level resource comparison: one row per sender, receiver and
        ligand-receptor pair significant in either resource, flagged per resource.
    """
    key = ["source", "target", "ligand_complex", "receptor_complex"]
    parts = []
    for resource in ccu.RESOURCES:
        frame = significant[resource][key].drop_duplicates().copy()
        frame[resource] = True
        parts.append(frame.set_index(key))

    merged = pd.concat(parts, axis = 1).fillna(False).reset_index()
    for resource in ccu.RESOURCES:
        merged[resource] = merged[resource].astype(bool)

    memberships = []
    for record in merged.itertuples(index = False):
        flags = [name for name in ccu.RESOURCES if getattr(record, name)]
        memberships.append("both" if len(flags) > 1 else flags[0])
    merged["membership"] = memberships
    return merged.sort_values(["membership", "source", "target"], ignore_index = True)


def node_layout(weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
        Start and end angle of every node on the circle, sized by weight and
        separated by a fixed gap.
    """
    total = float(weights.sum())
    usable = 2 * np.pi - NODE_GAP * len(weights)
    spans = weights * usable / total if total > 0 else np.zeros(len(weights))
    starts = np.zeros(len(weights))
    cursor = 0.0
    for position, span in enumerate(spans):
        starts[position] = cursor
        cursor += span + NODE_GAP
    return starts, starts + spans


def arc_points(radius: float, start: float, end: float, steps: int = 24) -> np.ndarray:
    """Points along a circular arc between two angles."""
    angles = np.linspace(start, end, steps)
    return np.column_stack([radius * np.cos(angles), radius * np.sin(angles)])


def chord_path(radius: float, first: tuple[float, float], second: tuple[float, float]) -> Path:
    """Closed ribbon path joining two arc spans through the centre of the circle."""
    start = arc_points(radius, first[0], first[1])
    end = arc_points(radius, second[0], second[1])

    vertices = [tuple(start[0])]
    codes = [Path.MOVETO]
    for point in start[1:]:
        vertices.append(tuple(point))
        codes.append(Path.LINETO)
    vertices += [(0.0, 0.0), (0.0, 0.0), tuple(end[0])]
    codes += [Path.CURVE4, Path.CURVE4, Path.CURVE4]
    for point in end[1:]:
        vertices.append(tuple(point))
        codes.append(Path.LINETO)
    vertices += [(0.0, 0.0), (0.0, 0.0), tuple(start[0])]
    codes += [Path.CURVE4, Path.CURVE4, Path.CURVE4]
    vertices.append(tuple(start[0]))
    codes.append(Path.CLOSEPOLY)
    return Path(vertices, codes)


def plot_circle(counts: pd.DataFrame) -> Figure:
    """
        Circle plot of cross-signalling: cell types as nodes, edge width by
        interaction count, edge colour by sender, self-signalling omitted.
    """
    labels = list(counts.index)
    angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, len(labels), endpoint = False)
    points = np.column_stack([np.cos(angles), np.sin(angles)])

    cross = counts.to_numpy().astype(float).copy()
    np.fill_diagonal(cross, 0.0)
    present = cross[cross > 0]
    low = float(present.min()) if len(present) > 0 else 0.0
    high = float(present.max()) if len(present) > 0 else 1.0

    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    for row in range(len(labels)):
        for column in range(len(labels)):
            if cross[row, column] <= 0:
                continue
            arrow = FancyArrowPatch(
                tuple(points[row]), tuple(points[column]),
                connectionstyle = "arc3,rad=0.2",
                arrowstyle = "-|>", mutation_scale = 9,
                linewidth = edge_width(cross[row, column], low, high),
                color = label_colour(row), alpha = 0.55,
                shrinkA = 7, shrinkB = 9,
            )
            ax.add_patch(arrow)

    totals = counts.to_numpy().sum(axis = 1) + counts.to_numpy().sum(axis = 0)
    sizes = scale_between(totals.astype(float), 40.0, 320.0)
    for position, name in enumerate(labels):
        ax.scatter(
            points[position, 0], points[position, 1],
            s = float(sizes[position]), color = label_colour(position), zorder = 3,
        )
        outward_text(ax, 1.13, angles[position], name)

    ax.set_xlim(CIRCLE_XLIM)
    ax.set_ylim(CIRCLE_YLIM)
    ax.set_aspect("equal")
    ax.axis("off")
    circle_title(ax, f"Significant interactions between cell types ({MAIN_RESOURCE})")
    ccu.style_dark(fig)
    return fig


def plot_chord(counts: pd.DataFrame) -> Figure:
    """
        Chord diagram of cross-signalling: arc length by total interactions,
        ribbon width by interaction count, ribbon colour by sender.
    """
    labels = list(counts.index)
    cross = counts.to_numpy().astype(float).copy()
    np.fill_diagonal(cross, 0.0)

    weights = cross.sum(axis = 1) + cross.sum(axis = 0)
    keep = [position for position in range(len(labels)) if weights[position] > 0]
    kept_labels = [labels[position] for position in keep]
    kept = cross[np.ix_(keep, keep)]
    kept_weights = kept.sum(axis = 1) + kept.sum(axis = 0)

    starts, ends = node_layout(kept_weights)
    scale = (ends - starts) / np.where(kept_weights > 0, kept_weights, 1.0)

    outgoing: dict[tuple[int, int], tuple[float, float]] = {}
    incoming: dict[tuple[int, int], tuple[float, float]] = {}
    for row in range(len(keep)):
        cursor = starts[row]
        for column in range(len(keep)):
            width = kept[row, column] * scale[row]
            if width > 0:
                outgoing[(row, column)] = (cursor, cursor + width)
            cursor += width
        for column in range(len(keep)):
            width = kept[column, row] * scale[row]
            if width > 0:
                incoming[(column, row)] = (cursor, cursor + width)
            cursor += width

    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    radius = 1.0
    for pair, span in outgoing.items():
        target = incoming.get(pair)
        if target is None:
            continue
        ax.add_patch(
            PathPatch(
                chord_path(radius, span, target),
                facecolor = label_colour(keep[pair[0]]), edgecolor = "none", alpha = 0.5,
            )
        )

    for position, name in enumerate(kept_labels):
        ax.add_patch(
            Wedge(
                (0.0, 0.0), radius + ARC_WIDTH,
                np.degrees(starts[position]), np.degrees(ends[position]),
                width = ARC_WIDTH, facecolor = label_colour(keep[position]),
            )
        )
        outward_text(ax, radius + 0.13, (starts[position] + ends[position]) / 2, name)

    ax.set_xlim(CIRCLE_XLIM)
    ax.set_ylim(CIRCLE_YLIM)
    ax.set_aspect("equal")
    ax.axis("off")
    circle_title(ax, f"Interaction flow between cell types ({MAIN_RESOURCE})")
    ccu.style_dark(fig)
    return fig


def plot_heatmap(counts: pd.DataFrame) -> Figure:
    """
        Sender-by-receiver heatmap of interaction counts with marginal bars for
        the total each cell type sends and receives.
    """
    labels = list(counts.index)
    values = counts.to_numpy()

    fig = plt.figure(figsize = ccu.FIGSIZE)
    grid = fig.add_gridspec(
        2, 3, width_ratios = [6, 1, 0.22], height_ratios = [1, 6],
        wspace = 0.06, hspace = 0.06,
        left = 0.14, right = 0.93, top = 0.90, bottom = 0.22,
    )
    top = fig.add_subplot(grid[0, 0])
    main = fig.add_subplot(grid[1, 0])
    side = fig.add_subplot(grid[1, 1])
    scale_axis = fig.add_subplot(grid[1, 2])

    image = main.imshow(values, cmap = "magma", aspect = "auto", vmin = 0)
    main.set_xticks(range(len(labels)))
    main.set_xticklabels(labels, rotation = 90, fontsize = LABEL_FONTSIZE)
    main.set_yticks(range(len(labels)))
    main.set_yticklabels(labels, fontsize = LABEL_FONTSIZE)
    main.set_xlabel("receiver", fontsize = LABEL_FONTSIZE + 2)
    main.set_ylabel("sender", fontsize = LABEL_FONTSIZE + 2)

    top.bar(range(len(labels)), values.sum(axis = 0), color = ccu.DARK_FOREGROUND)
    top.set_xlim(-0.5, len(labels) - 0.5)
    top.set_xticks([])
    top.tick_params(labelsize = LABEL_FONTSIZE)
    top.set_ylabel("received", fontsize = LABEL_FONTSIZE)

    side.barh(range(len(labels)), values.sum(axis = 1), color = ccu.DARK_FOREGROUND)
    side.set_ylim(len(labels) - 0.5, -0.5)
    side.set_yticks([])
    side.tick_params(labelsize = LABEL_FONTSIZE)
    side.set_xlabel("sent", fontsize = LABEL_FONTSIZE)

    bar = fig.colorbar(image, cax = scale_axis)
    bar.set_label("interactions", color = ccu.DARK_FOREGROUND,
                  fontsize = LABEL_FONTSIZE + 1, labelpad = 12)
    bar.ax.tick_params(colors = ccu.DARK_FOREGROUND, labelsize = LABEL_FONTSIZE)

    fig.suptitle(
        f"Interaction counts per sender and receiver ({MAIN_RESOURCE})",
        color = ccu.DARK_FOREGROUND, fontsize = TITLE_FONTSIZE, y = 0.965,
    )
    ccu.style_dark(fig)
    return fig


def plot_hierarchy(counts: pd.DataFrame) -> Figure:
    """
        Hierarchy plot of signalling into each compartment: every cell type as a
        sender on the left, the receiving compartment on the right,
        self-signalling drawn apart from cross-signalling.
    """
    labels = list(counts.index)
    groups = (
        ("into the developing compartment",
         [name for name in labels if compartment(name) == "developing"]),
        ("into the finished compartment",
         [name for name in labels if compartment(name) == "finished"]),
    )

    present = counts.to_numpy()[counts.to_numpy() > 0]
    low = float(present.min()) if len(present) > 0 else 0.0
    high = float(present.max()) if len(present) > 0 else 1.0

    fig, axes = plt.subplots(1, 2, figsize = ccu.FIGSIZE)
    for axis, (title, targets) in zip(axes, groups):
        sender_y = np.linspace(0.0, 1.0, len(labels))
        target_y = np.linspace(0.0, 1.0, len(targets))

        for row, sender in enumerate(labels):
            for column, receiver in enumerate(targets):
                value = float(counts.loc[sender, receiver]) # type: ignore
                if value <= 0:
                    continue
                colour = "#e2504a" if sender == receiver else label_colour(row)
                arrow = FancyArrowPatch(
                    (0.0, sender_y[row]), (1.0, target_y[column]),
                    connectionstyle = "arc3,rad=0.12",
                    arrowstyle = "-|>", mutation_scale = 8,
                    linewidth = edge_width(value, low, high),
                    color = colour, alpha = 0.55,
                    shrinkA = 3, shrinkB = 5,
                )
                axis.add_patch(arrow)

        for row, sender in enumerate(labels):
            axis.scatter(0.0, sender_y[row], s = 30, color = label_colour(row), zorder = 3)
            axis.text(
                -0.05, sender_y[row], sender, ha = "right", va = "center",
                fontsize = LABEL_FONTSIZE, color = ccu.DARK_FOREGROUND,
            )
        for column, receiver in enumerate(targets):
            axis.scatter(1.0, target_y[column], s = 30, color = ccu.DARK_FOREGROUND, zorder = 3)
            axis.text(
                1.05, target_y[column], receiver, ha = "left", va = "center",
                fontsize = LABEL_FONTSIZE, color = ccu.DARK_FOREGROUND,
            )

        axis.set_xlim(-0.85, 1.85)
        axis.set_ylim(-0.06, 1.06)
        axis.axis("off")
        axis.set_title(title, color = ccu.DARK_FOREGROUND, fontsize = LABEL_FONTSIZE + 3, pad = 14)

    fig.suptitle(
        f"Signalling into the two compartments ({MAIN_RESOURCE}), self-signalling in red",
        color = ccu.DARK_FOREGROUND, fontsize = TITLE_FONTSIZE, y = 0.97,
    )
    fig.subplots_adjust(top = 0.86, bottom = 0.04, left = 0.02, right = 0.98, wspace = 0.30)
    ccu.style_dark(fig)
    return fig


def plot_dotplot(significant: pd.DataFrame) -> Figure:
    """
        Dotplot of the top interactions: cell type pair against ligand-receptor
        pair, sized by magnitude, coloured by specificity.
    """
    top = significant.nsmallest(TOP_INTERACTIONS, "magnitude_rank").copy()
    top["pair"] = interaction_pairs(top)
    top["cell_pair"] = top["source"] + " -> " + top["target"]

    pairs = list(top["pair"].value_counts().index)
    cell_pairs = list(top["cell_pair"].value_counts().index)
    x = [cell_pairs.index(name) for name in top["cell_pair"]]
    y = [pairs.index(name) for name in top["pair"]]

    magnitude = rank_score(top["magnitude_rank"])
    specificity = rank_score(top["specificity_rank"])

    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    dots = ax.scatter(
        x, y,
        s = scale_between(magnitude, DOT_SIZE[0], DOT_SIZE[1]),
        c = specificity, cmap = "plasma",
    )
    ax.set_xticks(range(len(cell_pairs)))
    ax.set_xticklabels(cell_pairs, rotation = 90, fontsize = LABEL_FONTSIZE)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs, fontsize = LABEL_FONTSIZE)
    ax.set_xlabel("sender and receiver", fontsize = LABEL_FONTSIZE + 2)
    ax.set_ylabel("ligand and receptor", fontsize = LABEL_FONTSIZE + 2)
    ax.set_title(
        f"Top {TOP_INTERACTIONS} interactions by magnitude ({MAIN_RESOURCE})",
        color = ccu.DARK_FOREGROUND, fontsize = TITLE_FONTSIZE, pad = 14,
    )

    bar = fig.colorbar(dots, ax = ax, fraction = 0.030, pad = 0.055)
    bar.set_label("specificity, negative log10 rank", color = ccu.DARK_FOREGROUND,
                  fontsize = LABEL_FONTSIZE + 1, labelpad = 18)
    bar.ax.tick_params(colors = ccu.DARK_FOREGROUND, labelsize = LABEL_FONTSIZE)

    ccu.style_dark(fig)
    fig.subplots_adjust(left = 0.13, right = 0.88, top = 0.93, bottom = 0.30)
    return fig


def plot_sources_targets(counts: pd.DataFrame) -> Figure:
    """
        Interactions sent and received per cell type, one row each with the two
        totals joined, ordered by the difference between them.
    """
    labels = list(counts.index)
    sent = counts.to_numpy().sum(axis = 1).astype(float)
    received = counts.to_numpy().sum(axis = 0).astype(float)
    net = received - sent
    order = np.argsort(net)
    rows = np.arange(len(order))

    limit = float(max(sent.max(), received.max())) * 1.22

    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    for row, position in enumerate(order):
        ax.plot(
            [sent[position], received[position]], [row, row],
            color = ccu.DARK_FOREGROUND, alpha = 0.3, linewidth = 1.6, zorder = 1,
        )
    ax.scatter(
        sent[order], rows, s = 120, color = SENT_COLOUR,
        label = "sent", zorder = 3,
    )
    ax.scatter(
        received[order], rows, s = 120, color = RECEIVED_COLOUR,
        label = "received", zorder = 3,
    )

    for row, position in enumerate(order):
        ax.text(
            limit * 0.98, row, f"{net[position]:+.0f}",
            ha = "right", va = "center",
            fontsize = LABEL_FONTSIZE, color = ccu.DARK_FOREGROUND,
        )
    ax.text(
        limit * 0.98, len(order) - 0.45, "net",
        ha = "right", va = "center",
        fontsize = LABEL_FONTSIZE, color = ccu.DARK_FOREGROUND,
    )

    ax.set_yticks(rows)
    ax.set_yticklabels([labels[position] for position in order], fontsize = LABEL_FONTSIZE + 1)
    for tick, position in zip(ax.get_yticklabels(), order):
        tick.set_color(COMPARTMENT_COLOURS[compartment(labels[position])])

    ax.set_xlim(-1.5, limit)
    ax.set_ylim(-0.8, len(order))
    ax.tick_params(axis = "x", labelsize = LABEL_FONTSIZE)
    ax.grid(axis = "x", alpha = 0.12, linewidth = 0.7)
    ax.set_axisbelow(True)
    ax.set_xlabel("significant interactions", fontsize = LABEL_FONTSIZE + 2)
    ax.set_title(
        f"Interactions sent and received per cell type ({MAIN_RESOURCE}), "
        f"names coloured by compartment",
        color = ccu.DARK_FOREGROUND, fontsize = TITLE_FONTSIZE, pad = 14,
    )
    ax.legend(fontsize = LEGEND_FONTSIZE, loc = "lower right", framealpha = 0.15)

    ccu.style_dark(fig)
    fig.subplots_adjust(left = 0.16, right = 0.97, top = 0.92, bottom = 0.08)
    return fig


def plot_resource_overlap(overlap: pd.DataFrame) -> Figure:
    """
        Ligand-receptor pairs significant in either resource, with the count
        each resource contributes drawn side by side.
    """
    frame = overlap.copy()
    frame["pair"] = interaction_pairs(frame)
    table = frame.groupby("pair")[list(ccu.RESOURCES)].sum()
    table = table.sort_values(list(ccu.RESOURCES), ascending = True)

    positions = np.arange(len(table))
    height = 0.38

    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    ax.barh(
        positions + height / 2, table[ccu.RESOURCES[0]], height = height,
        color = "#5ab1bb", label = ccu.RESOURCES[0],
    )
    ax.barh(
        positions - height / 2, table[ccu.RESOURCES[1]], height = height,
        color = "#f2b134", label = ccu.RESOURCES[1],
    )
    ax.set_yticks(positions)
    ax.set_yticklabels(table.index, fontsize = LABEL_FONTSIZE)
    ax.tick_params(axis = "x", labelsize = LABEL_FONTSIZE)
    ax.set_xlabel("significant interactions", fontsize = LABEL_FONTSIZE + 2)
    ax.set_ylabel("ligand and receptor", fontsize = LABEL_FONTSIZE + 2)

    shared = int((overlap["membership"] == "both").sum())
    first_only = int((overlap["membership"] == ccu.RESOURCES[0]).sum())
    second_only = int((overlap["membership"] == ccu.RESOURCES[1]).sum())
    ax.set_title(
        f"Resource comparison: {shared} interactions significant in both, "
        f"{first_only} in {ccu.RESOURCES[0]} only, {second_only} in {ccu.RESOURCES[1]} only",
        color = ccu.DARK_FOREGROUND, fontsize = TITLE_FONTSIZE, pad = 14,
    )
    ax.legend(fontsize = LEGEND_FONTSIZE, loc = "lower right", framealpha = 0.15)

    ccu.style_dark(fig)
    fig.subplots_adjust(left = 0.16, right = 0.97, top = 0.93, bottom = 0.08)
    return fig


def main() -> None:
    """
        Summary and resource comparison tables, and the seven cross-cutting
        figures of the stage.
    """
    ccu.set_seed()

    significant = {}
    for resource in ccu.RESOURCES:
        significant[resource] = load_significant(resource)

    labels = ordered_labels(load_all(MAIN_RESOURCE))
    counts = counts_matrix(significant[MAIN_RESOURCE], labels)

    summary = summary_table(counts)
    ccu.save_table(summary, SCRIPT_NAME, f"{SCRIPT_NAME}_summary")

    overlap = overlap_table(significant)
    ccu.save_table(overlap, SCRIPT_NAME, f"{SCRIPT_NAME}_resource_overlap")

    print(
        f"{MAIN_RESOURCE}: {int(counts.to_numpy().sum())} significant interactions over "
        f"{len(labels)} cell types, {int(np.trace(counts.to_numpy()))} of them self-signalling"
    )
    print(
        f"resources: {int((overlap['membership'] == 'both').sum())} shared, "
        f"{len(overlap)} in either"
    )

    ccu.save_figure(plot_circle(counts), f"{SCRIPT_NAME}_circle")
    ccu.save_figure(plot_chord(counts), f"{SCRIPT_NAME}_chord")
    ccu.save_figure(plot_heatmap(counts), f"{SCRIPT_NAME}_heatmap")
    ccu.save_figure(plot_hierarchy(counts), f"{SCRIPT_NAME}_hierarchy")
    ccu.save_figure(plot_dotplot(significant[MAIN_RESOURCE]), f"{SCRIPT_NAME}_dotplot")
    ccu.save_figure(plot_sources_targets(counts), f"{SCRIPT_NAME}_sources_targets")
    ccu.save_figure(plot_resource_overlap(overlap), f"{SCRIPT_NAME}_resource_overlap")


if __name__ == "__main__":
    main()