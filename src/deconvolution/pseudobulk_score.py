"""Agreement between predicted and true cell-type fractions on the test
mixtures, per test set and donor, at label and broad-group level."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import deconvolution_utils as du

SCRIPT        = "pseudobulk_score"
TRUTH         = ("pseudobulk_test", "pseudobulk_test_fractions.csv")
MODELS        = {
    "scaden": ("Scaden",
               ("scaden_fit", "scaden_fit_predictions.csv")),
    "scaden_torch": ("Scaden in PyTorch",
                     ("scaden_torch_fit", "scaden_torch_fit_predictions.csv")),
}

LEVELS        = {"label": pd.DataFrame.copy, "group": du.broad_groups}
DONOR_COLORS  = dict(zip(du.TEST_DONORS, ("#5aa9e6", "#f4a259"),
                         strict=True))
SCATTER_COLS  = 4
PANEL_SIZE    = 4.0
ANNOT_SIZE    = 11


def agreement(true: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    """Pearson correlation, Lin's concordance correlation, RMSE, mean absolute
    error and mean bias of predicted against true fractions."""
    covariance = np.mean((true - true.mean()) * (predicted - predicted.mean()))
    spread = true.var() + predicted.var() + (true.mean() - predicted.mean())**2
    varies = true.std() > 0 and predicted.std() > 0
    return {
        "pearson": float(np.corrcoef(true, predicted)[0, 1]) if varies
        else np.nan,
        "ccc": float(2 * covariance / spread) if spread > 0 else np.nan,
        "rmse": float(np.sqrt(np.mean((predicted - true) ** 2))),
        "mae": float(np.mean(np.abs(predicted - true))),
        "bias": float(np.mean(predicted - true)),
    }


def score_table(truth: pd.DataFrame, predicted: pd.DataFrame,
                subsets: dict[str, np.ndarray]) -> pd.DataFrame:
    """Per subset and cell type agreement with true and predicted means and
    relative error; 'all' pools every cell type of the subset."""
    rows = []
    for subset, mask in subsets.items():
        true, pred = truth[mask], predicted[mask]
        for cell_type in [*truth.columns, "all"]:
            columns = truth.columns if cell_type == "all" else [cell_type]
            values = true[columns].to_numpy().ravel()
            rows.append({
                "subset": subset, "cell_type": cell_type,
                "n": int(mask.sum()), "true_mean": values.mean(),
                "predicted_mean": pred[columns].to_numpy().mean(),
                **agreement(values, pred[columns].to_numpy().ravel()),
            })
    table = pd.DataFrame(rows)
    table["relative_mae"] = (table["mae"] / table["true_mean"]
                             .where(table["true_mean"] > 0))
    return table


def plot_heatmap(table: pd.DataFrame, title: str) -> Figure:
    """Annotated heatmap of concordance correlation: cell type against test
    subset on a fixed 0 to 1 colour scale."""
    matrix = table[table["cell_type"] != "all"].pivot(
        index="cell_type", columns="subset", values="ccc",
    ).reindex(index=table["cell_type"].unique()[:-1],
              columns=table["subset"].unique())
    figure, axis = plt.subplots(
        figsize=(13, 2.5 + 0.55 * len(matrix)), layout="constrained",
    )
    du.style_dark(figure, [axis])
    image = axis.imshow(matrix.to_numpy(), cmap="magma", vmin=0, vmax=1,
                        aspect="auto")
    for (row, col), value in np.ndenumerate(matrix.to_numpy()):
        axis.text(col, row, "n/a" if np.isnan(value) else f"{value:.2f}",
                  ha="center", va="center", fontsize=ANNOT_SIZE,
                  color="black" if value > 0.6 else du.FOREGROUND)
    axis.set_xticks(range(matrix.shape[1]),
                    [name.replace(" ", "\n") for name in matrix.columns])
    axis.set_yticks(range(matrix.shape[0]), matrix.index)
    axis.set_title(title, fontsize=du.TITLE_FONTSIZE, pad=du.TITLE_PAD)
    bar = figure.colorbar(image, ax=axis, shrink=0.8, pad=0.02)
    bar.set_label("Concordance correlation (CCC)", color=du.FOREGROUND,
                  fontsize=du.LABEL_FONTSIZE)
    bar.ax.tick_params(colors=du.FOREGROUND, labelsize=du.TICK_FONTSIZE)
    return figure


def plot_scatter(truth: pd.DataFrame, predicted: pd.DataFrame,
                 donors: pd.Series, table: pd.DataFrame,
                 title: str) -> Figure:
    """Predicted against true fraction per cell type, one square panel each
    with its own range and the identity line, points coloured by donor."""
    n_rows = math.ceil(truth.shape[1] / SCATTER_COLS)
    figure, axes = plt.subplots(
        n_rows, SCATTER_COLS, layout="constrained",
        figsize=(SCATTER_COLS * PANEL_SIZE, n_rows * PANEL_SIZE + 0.8),
    )
    du.style_dark(figure, list(axes.flat))
    ccc = table.set_index("cell_type")["ccc"]
    for axis, cell_type in zip(axes.flat, truth.columns, strict=False):
        limit = 1.05 * max(truth[cell_type].max(),
                           predicted[cell_type].max(), 0.01)
        axis.plot([0, limit], [0, limit], color=du.FOREGROUND,
                  linestyle="--", linewidth=1)
        for donor, color in DONOR_COLORS.items():
            rows = donors == donor
            axis.scatter(truth.loc[rows, cell_type],
                         predicted.loc[rows, cell_type], s=10, alpha=0.6,
                         color=color, linewidths=0)
        axis.set(xlim=(0, limit), ylim=(0, limit), aspect="equal")
        axis.set_title(cell_type, fontsize=du.LABEL_FONTSIZE, pad=8)
        axis.text(0.04, 0.95, f"CCC {ccc[cell_type]:.2f}",
                  transform=axis.transAxes, va="top",
                  fontsize=du.TICK_FONTSIZE, color=du.FOREGROUND,
                  bbox={"facecolor": du.BACKGROUND, "edgecolor": "none",
                        "alpha": 0.85, "pad": 2})
    for axis in axes.flat[truth.shape[1]:]:
        axis.set_visible(False)
    figure.supxlabel("True fraction", fontsize=du.LABEL_FONTSIZE,
                     color=du.FOREGROUND)
    figure.supylabel("Predicted fraction", fontsize=du.LABEL_FONTSIZE,
                     color=du.FOREGROUND)
    figure.suptitle(title, fontsize=du.TITLE_FONTSIZE, color=du.FOREGROUND)
    legend = figure.legend(
        handles=[Line2D([], [], marker="o", linestyle="", markersize=9,
                        color=color, label=donor)
                 for donor, color in DONOR_COLORS.items()],
        title="Donor", loc="outside right upper", frameon=False,
        fontsize=du.LEGEND_FONTSIZE, title_fontsize=du.LEGEND_FONTSIZE,
        labelcolor=du.FOREGROUND,
    )
    legend.get_title().set_color(du.FOREGROUND)
    return figure


def plot_bias(table: pd.DataFrame, title: str) -> Figure:
    """Mean predicted minus mean true fraction per cell type, one bar per
    donor, on a symmetric axis around zero."""
    bias = table[table["cell_type"] != "all"].pivot(
        index="cell_type", columns="subset", values="bias",
    ).reindex(index=table["cell_type"].unique()[:-1])
    bias = bias[[f"donor {donor}" for donor in DONOR_COLORS]]
    figure, axis = plt.subplots(figsize=(14, 2.5 + 0.5 * len(bias)),
                                layout="constrained")
    du.style_dark(figure, [axis])
    positions = np.arange(len(bias))
    for offset, (donor, color) in zip((-0.2, 0.2), DONOR_COLORS.items(),
                                      strict=True):
        axis.barh(positions + offset, bias[f"donor {donor}"], height=0.38,
                  color=color, label=donor)
    limit = 1.1 * np.nanmax(np.abs(bias.to_numpy()))
    axis.axvline(0, color=du.FOREGROUND, linewidth=1)
    axis.set(xlim=(-limit, limit), ylim=(len(bias) - 0.5, -0.5))
    axis.set_yticks(positions, bias.index)
    axis.set_xlabel("Mean predicted minus mean true fraction",
                    fontsize=du.LABEL_FONTSIZE)
    axis.set_title(title, fontsize=du.TITLE_FONTSIZE, pad=du.TITLE_PAD)
    legend = axis.legend(title="Donor", loc="upper left",
                         bbox_to_anchor=(1.01, 1.0), frameon=False,
                         fontsize=du.LEGEND_FONTSIZE,
                         title_fontsize=du.LEGEND_FONTSIZE,
                         labelcolor=du.FOREGROUND)
    legend.get_title().set_color(du.FOREGROUND)
    return figure


def load_truth() -> tuple[pd.DataFrame, pd.DataFrame,
                          dict[str, np.ndarray]]:
    """True test fractions, test-set and donor annotation, and a boolean mask
    per test subset."""
    frame = pd.read_csv(du.RESULTS_DIR.joinpath(*TRUTH), index_col=0)
    meta, truth = frame[["set", "donor"]], frame.drop(columns=["set", "donor"])
    subsets = {name: (meta["set"] == name).to_numpy()
               for name in meta["set"].unique()}
    subsets |= {f"donor {donor}": (meta["donor"] == donor).to_numpy()
                for donor in du.TEST_DONORS}
    return truth, meta, subsets


def main() -> None:
    """Scoring of every model at both levels, with tables and figures."""
    truth, meta, subsets = load_truth()
    donor_rows = subsets["donor"]

    tables = []
    for model, (label, source) in MODELS.items():
        predicted = pd.read_csv(du.RESULTS_DIR.joinpath(*source),
                                index_col=0).loc[truth.index, truth.columns]
        for level, transform in LEVELS.items():
            true_level, pred_level = transform(truth), transform(predicted)
            table = score_table(true_level, pred_level, subsets)
            tables.append(table.assign(model=model, level=level))
            name = f"{SCRIPT}_{model}_{level}"
            du.save_figure(plot_heatmap(
                table, f"{label}: agreement with the true "
                f"fractions per {level} and test subset"), SCRIPT,
                f"{name}_ccc")
            du.save_figure(plot_scatter(
                true_level[donor_rows], pred_level[donor_rows],
                meta.loc[donor_rows, "donor"],
                table[table["subset"] == "donor"],
                f"{label}: predicted against true fraction "
                f"per {level}, donor-composition test set"), SCRIPT,
                f"{name}_scatter")
            du.save_figure(plot_bias(
                table, f"{label}: mean bias per {level}, "
                "donor-composition test set"), SCRIPT, f"{name}_bias")

    scores = pd.concat(tables)[["model", "level", *tables[0].columns[:-2]]]
    du.save_table(scores, SCRIPT, SCRIPT, index=False)


if __name__ == "__main__":
    main()