"""Agreement of the Scaden predictions on sorted bulk RNA-seq (GSE74246) with
the known group shares: pure sorted samples and known mixtures."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import deconvolution_utils as du
import pseudobulk_score as ps

SCRIPT       = "bulk_score"
SOURCE       = "scaden_torch_bulk"
VARIANTS     = ("counts", "length")
TRUTH        = ("corces_bulk", "corces_bulk_truth.csv")
GROUPS       = ["HSPC", "Erythroid", "Myeloid", "B lineage", "T cells", "NK"]
IN_BETWEEN   = ("GMP", "MEP", "CLP")
TYPE_ORDER   = ("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "Ery",
                "Mono", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell")
KIND_COLORS  = {"random": "#5aa9e6", "marrow": "#f4a259"}
ENSEMBLE     = "ensemble"
METRICS      = ["ccc", "pearson", "bias", "rmse"]
SCATTER_COLS = 3
PANEL_SIZE   = 4.6
ANNOT_SIZE   = 12


def load_groups() -> pd.DataFrame:
    """Predicted group shares of every bulk variant, seed, training length and
    sample."""
    predictions = pd.concat(
        {variant: pd.read_csv(du.RESULTS_DIR / SOURCE /
                              f"{SOURCE}_{variant}_predictions.csv",
                              index_col=[0, 1, 2])
         for variant in VARIANTS},
        names=["bulk"],
    )
    return du.broad_groups(predictions)


def pure_shares(groups: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    """Mean predicted group shares per bulk variant, training length and
    sorted type, with the share given to the type's own group."""
    pure = truth[truth["kind"] == "pure"]
    frame = groups.join(pure["sorted_type"], on="sample", how="inner")
    shares = frame.groupby(["bulk", "steps", "sorted_type"])[
        list(du.GROUP_ORDER)].mean()
    own = pure.groupby("sorted_type")[list(du.GROUP_ORDER)].first()
    types = shares.index.get_level_values("sorted_type")
    shares["own_group_share"] = (shares * own.reindex(types).to_numpy()
                                 ).sum(axis=1)
    shares["in_between"] = types.isin(IN_BETWEEN)
    return shares


def mixture_scores(groups: pd.DataFrame,
                   truth: pd.DataFrame) -> pd.DataFrame:
    """Mixture agreement per bulk variant, training length, mixture kind and
    seed, and for the mean over seeds."""
    mixtures = truth[truth["kind"] != "pure"]
    subsets = {kind: (mixtures["kind"] == kind).to_numpy()
               for kind in KIND_COLORS}
    ensemble = groups.groupby(level=["bulk", "steps", "sample"]).mean()
    runs = [(key, frame.droplevel(["bulk", "seed", "steps"]))
            for key, frame in groups.groupby(level=["bulk", "seed", "steps"])]
    runs += [((bulk, ENSEMBLE, steps), frame.droplevel(["bulk", "steps"]))
             for (bulk, steps), frame in ensemble.groupby(
                 level=["bulk", "steps"])]
    tables = []
    for (bulk, seed, steps), frame in runs:
        table = ps.score_table(mixtures[list(du.GROUP_ORDER)],
                               frame.loc[mixtures.index,
                                         list(du.GROUP_ORDER)], subsets)
        tables.append(table.assign(bulk=bulk, seed=str(seed), steps=steps))
    return pd.concat(tables, ignore_index=True)


def summarise(scores: pd.DataFrame) -> pd.DataFrame:
    """Per bulk variant, training length, mixture kind and group: median,
    minimum and maximum over seeds, and the value for the seed mean."""
    keys = ["bulk", "steps", "subset", "cell_type"]
    kept = scores[scores["cell_type"].isin([*GROUPS, "all"])]
    per_seed = kept[kept["seed"] != ENSEMBLE].groupby(keys)[METRICS].agg(
        ["median", "min", "max"])
    per_seed.columns = ["_".join(column) for column in per_seed.columns]
    ensemble = kept[kept["seed"] == ENSEMBLE].set_index(keys)[METRICS]
    return per_seed.join(ensemble.add_suffix(f"_{ENSEMBLE}"))


def plot_pure(shares: pd.DataFrame, variant: str, steps: int) -> Figure:
    """Annotated heatmap of the mean predicted group shares per sorted type
    on a fixed 0 to 1 colour scale."""
    matrix = shares.loc[(variant, steps)].reindex(TYPE_ORDER)[GROUPS]
    figure, axis = plt.subplots(figsize=(13, 3 + 0.6 * len(matrix)),
                                layout="constrained")
    du.style_dark(figure, [axis])
    image = axis.imshow(matrix.to_numpy(), cmap="magma", vmin=0, vmax=1,
                        aspect="auto")
    for (row, col), value in np.ndenumerate(matrix.to_numpy()):
        axis.text(col, row, f"{value:.2f}", ha="center", va="center",
                  fontsize=ANNOT_SIZE,
                  color="black" if value > 0.6 else du.FOREGROUND)
    axis.set_xticks(range(len(GROUPS)), GROUPS)
    axis.set_yticks(range(len(matrix)),
                    [f"{name} (in-between)" if name in IN_BETWEEN else name
                     for name in matrix.index])
    axis.set_title(f"Pure sorted samples ({variant}): mean predicted group "
                   f"share, {steps:,} training steps",
                   fontsize=du.TITLE_FONTSIZE, pad=du.TITLE_PAD)
    bar = figure.colorbar(image, ax=axis, shrink=0.8, pad=0.02)
    bar.set_label("Predicted share", color=du.FOREGROUND,
                  fontsize=du.LABEL_FONTSIZE)
    bar.ax.tick_params(colors=du.FOREGROUND, labelsize=du.TICK_FONTSIZE)
    return figure


def plot_mixtures(groups: pd.DataFrame, truth: pd.DataFrame,
                  scores: pd.DataFrame, variant: str, steps: int) -> Figure:
    """Predicted (mean over seeds) against true group share in the mixtures,
    one square panel per group with the identity line."""
    mixtures = truth[truth["kind"] != "pure"]
    predicted = groups.loc[variant].xs(steps, level="steps").groupby(
        level="sample").mean().loc[mixtures.index]
    ccc = scores[(scores["bulk"] == variant) & (scores["steps"] == steps)
                 & (scores["seed"] == ENSEMBLE)].set_index(
        ["subset", "cell_type"])["ccc"]
    n_rows = int(np.ceil(len(GROUPS) / SCATTER_COLS))
    figure, axes = plt.subplots(
        n_rows, SCATTER_COLS, layout="constrained",
        figsize=(SCATTER_COLS * PANEL_SIZE + 2, n_rows * PANEL_SIZE + 1),
    )
    du.style_dark(figure, list(axes.flat))
    for axis, group in zip(axes.flat, GROUPS, strict=True):
        limit = 1.05 * max(mixtures[group].max(), predicted[group].max())
        axis.plot([0, limit], [0, limit], color=du.FOREGROUND,
                  linestyle="--", linewidth=1)
        for kind, color in KIND_COLORS.items():
            rows = mixtures["kind"] == kind
            axis.scatter(mixtures.loc[rows, group],
                         predicted.loc[rows, group], s=12, alpha=0.6,
                         color=color, linewidths=0)
        axis.set(xlim=(0, limit), ylim=(0, limit), aspect="equal")
        axis.set_title(group, fontsize=du.LABEL_FONTSIZE, pad=8)
        axis.text(0.04, 0.95, "\n".join(
            f"CCC {kind} {ccc[kind, group]:.2f}" for kind in KIND_COLORS),
            transform=axis.transAxes, va="top", fontsize=du.TICK_FONTSIZE,
            color=du.FOREGROUND,
            bbox={"facecolor": du.BACKGROUND, "edgecolor": "none",
                  "alpha": 0.85, "pad": 2})
    figure.supxlabel("True group share", fontsize=du.LABEL_FONTSIZE,
                     color=du.FOREGROUND)
    figure.supylabel("Predicted group share (mean over seeds)",
                     fontsize=du.LABEL_FONTSIZE, color=du.FOREGROUND)
    figure.suptitle(f"Mixtures of sorted bulk profiles ({variant}), "
                    f"{steps:,} training steps", fontsize=du.TITLE_FONTSIZE,
                    color=du.FOREGROUND)
    legend = figure.legend(
        handles=[Line2D([], [], marker="o", linestyle="", markersize=9,
                        color=color, label=kind)
                 for kind, color in KIND_COLORS.items()],
        title="Mixture kind", loc="outside right upper", frameon=False,
        fontsize=du.LEGEND_FONTSIZE, title_fontsize=du.LEGEND_FONTSIZE,
        labelcolor=du.FOREGROUND,
    )
    legend.get_title().set_color(du.FOREGROUND)
    return figure


def main() -> None:
    """Pure-sample and mixture tables for both bulk variants, and figures for
    the longest training."""
    groups = load_groups()
    truth = pd.read_csv(du.RESULTS_DIR.joinpath(*TRUTH), index_col=0)
    shares = pure_shares(groups, truth)
    scores = mixture_scores(groups, truth)
    longest = int(groups.index.get_level_values("steps").max())
    du.save_table(shares, SCRIPT, f"{SCRIPT}_pure")
    du.save_table(scores, SCRIPT, f"{SCRIPT}_mixtures", index=False)
    du.save_table(summarise(scores), SCRIPT, f"{SCRIPT}_summary")
    for variant in VARIANTS:
        du.save_figure(plot_pure(shares, variant, longest), SCRIPT,
                       f"{SCRIPT}_pure_{variant}")
        du.save_figure(plot_mixtures(groups, truth, scores, variant,
                                     longest),
                       SCRIPT, f"{SCRIPT}_mixtures_{variant}")


if __name__ == "__main__":
    main()