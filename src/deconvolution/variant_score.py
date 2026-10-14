"""Agreement of the Scaden training variants with the true test fractions:
shuffle order and training length, summarised over seeds."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import deconvolution_utils as du
import pseudobulk_score as ps

SCRIPT        = "variant_score"
SOURCE        = ("scaden_torch_variants",
                 "scaden_torch_variants_predictions.csv")
RUN_KEYS      = ["order", "seed", "steps"]
ORDER_STYLE   = {
    "buffer": ("#5aa9e6", "1,000-row buffer shuffle (Scaden)"),
    "full":   ("#f4a259", "full shuffle"),
}
TREND_SUBSETS = {
    "donor":           "Donor set, both donors",
    "donor MantonBM1": "Donor set, MantonBM1",
    "donor MantonBM8": "Donor set, MantonBM8",
}
COMPARED      = {
    ("buffer", 5000): ("#5aa9e6", "buffer shuffle, 5,000 steps (Scaden)"),
    ("full", 20000):  ("#f4a259", "full shuffle, 20,000 steps"),
}
OFFSET        = 0.08


def score_runs(predictions: pd.DataFrame, truth: pd.DataFrame,
               subsets: dict[str, np.ndarray]) -> pd.DataFrame:
    """Agreement table of every run (shuffle order, seed and training length)
    at label and group level."""
    tables = []
    for run, frame in predictions.groupby(level=RUN_KEYS, sort=False):
        predicted = frame.droplevel(RUN_KEYS).loc[truth.index, truth.columns]
        for level, transform in ps.LEVELS.items():
            table = ps.score_table(transform(truth), transform(predicted),
                                   subsets)
            tables.append(table.assign(**dict(zip(RUN_KEYS, run,
                                                  strict=True)),
                                       level=level))
    return pd.concat(tables, ignore_index=True)


def summarise(scores: pd.DataFrame) -> pd.DataFrame:
    """Median per-label CCC and pooled RMSE per level, subset, shuffle order
    and training length, as median, minimum and maximum over seeds."""
    keys = ["level", "subset", *RUN_KEYS]
    per_run = pd.concat([
        scores[scores["cell_type"] != "all"].groupby(keys)["ccc"].median()
        .rename("median_ccc"),
        scores[scores["cell_type"] == "all"].set_index(keys)["rmse"]
        .rename("pooled_rmse"),
    ], axis=1)
    summary = per_run.groupby(level=["level", "subset", "order",
                                     "steps"]).agg(["median", "min", "max"])
    summary.columns = ["_".join(column) for column in summary.columns]
    return summary


def plot_trend(summary: pd.DataFrame, level: str) -> Figure:
    """Median per-label CCC against training length per donor subset, one
    line per shuffle order with the range over seeds."""
    figure, axes = plt.subplots(1, len(TREND_SUBSETS), figsize=(15, 6.5),
                                sharey=True, layout="constrained")
    du.style_dark(figure, list(axes))
    for axis, (subset, title) in zip(axes, TREND_SUBSETS.items(),
                                     strict=True):
        for rank, (order, (color, _)) in enumerate(ORDER_STYLE.items()):
            rows = summary.loc[(level, subset, order)]
            x = np.arange(len(rows)) + (rank - 0.5) * OFFSET
            middle = rows["median_ccc_median"].to_numpy()
            axis.errorbar(
                x, middle, color=color, marker="o", markersize=7, lw=2,
                capsize=5,
                yerr=[middle - rows["median_ccc_min"].to_numpy(),
                      rows["median_ccc_max"].to_numpy() - middle],
            )
        axis.set_xticks(range(len(rows)), [f"{steps:,}"
                                           for steps in rows.index])
        axis.set(xlim=(-0.4, len(rows) - 0.6), ylim=(0, 1))
        axis.set_title(title, fontsize=du.TITLE_FONTSIZE, pad=du.TITLE_PAD)
        axis.set_xlabel("Training steps", fontsize=du.LABEL_FONTSIZE)
        axis.grid(color=du.FOREGROUND, alpha=0.15)
    axes[0].set_ylabel(f"Median CCC over {level}s",
                       fontsize=du.LABEL_FONTSIZE)
    figure.suptitle(f"Training variants per {level}: median over seeds, "
                    "bars show the range over seeds",
                    fontsize=du.TITLE_FONTSIZE, color=du.FOREGROUND)
    legend = figure.legend(
        handles=[Line2D([], [], color=color, marker="o", lw=2, label=label)
                 for color, label in ORDER_STYLE.values()],
        title="Shuffle", loc="outside lower center", ncols=2, frameon=False,
        fontsize=du.LEGEND_FONTSIZE, title_fontsize=du.LEGEND_FONTSIZE,
        labelcolor=du.FOREGROUND,
    )
    legend.get_title().set_color(du.FOREGROUND)
    return figure


def plot_labels(scores: pd.DataFrame) -> Figure:
    """Per-label CCC on the donor-composition test set for the Scaden default
    and the longest full-shuffle training, median over seeds."""
    donor = scores[(scores["level"] == "label")
                   & (scores["subset"] == "donor")
                   & (scores["cell_type"] != "all")]
    medians = donor.groupby(["cell_type", "order", "steps"])["ccc"].median()
    labels = list(dict.fromkeys(donor["cell_type"]))
    figure, axis = plt.subplots(figsize=(15, 3 + 0.45 * len(labels)),
                                layout="constrained")
    du.style_dark(figure, [axis])
    positions = np.arange(len(labels))
    values = {run: medians.xs(run, level=["order", "steps"])
              .reindex(labels).to_numpy() for run in COMPARED}
    axis.hlines(positions, *np.sort(np.vstack(list(values.values())),
                                    axis=0), color=du.FOREGROUND, lw=1,
                alpha=0.5)
    for run, (color, _) in COMPARED.items():
        axis.scatter(values[run], positions, s=70, color=color, zorder=3)
    axis.axvline(0, color=du.FOREGROUND, lw=1)
    axis.set_yticks(positions, labels)
    axis.set(xlim=(-0.2, 1), ylim=(len(labels) - 0.5, -0.5))
    axis.set_xlabel("Concordance correlation (CCC), median over seeds",
                    fontsize=du.LABEL_FONTSIZE)
    axis.set_title("Per-label agreement on the donor-composition test set",
                   fontsize=du.TITLE_FONTSIZE, pad=du.TITLE_PAD)
    axis.grid(axis="x", color=du.FOREGROUND, alpha=0.15)
    legend = axis.legend(
        handles=[Line2D([], [], color=color, marker="o", linestyle="",
                        markersize=9, label=label)
                 for color, label in COMPARED.values()],
        title="Training", loc="upper left", bbox_to_anchor=(1.01, 1.0),
        frameon=False, fontsize=du.LEGEND_FONTSIZE,
        title_fontsize=du.LEGEND_FONTSIZE, labelcolor=du.FOREGROUND,
    )
    legend.get_title().set_color(du.FOREGROUND)
    return figure


def main() -> None:
    """Scores of every variant run, summary over seeds and figures."""
    truth, _, subsets = ps.load_truth()
    predictions = pd.read_csv(du.RESULTS_DIR.joinpath(*SOURCE),
                              index_col=[0, 1, 2, 3])
    scores = score_runs(predictions, truth, subsets)
    summary = summarise(scores)
    du.save_table(scores, SCRIPT, SCRIPT, index=False)
    du.save_table(summary, SCRIPT, f"{SCRIPT}_summary")
    for level in ps.LEVELS:
        du.save_figure(plot_trend(summary, level), SCRIPT,
                       f"{SCRIPT}_{level}_trend")
    du.save_figure(plot_labels(scores), SCRIPT, f"{SCRIPT}_labels")


if __name__ == "__main__":
    main()