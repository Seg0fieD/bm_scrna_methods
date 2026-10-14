"""Donor-level checks: haemoglobin transcript share outside erythroid cells
and immediate-early gene expression in T cells."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.image import AxesImage

import deconvolution_utils as du

SCRIPT           = "donor_diagnostics"
DONORS           = tuple(sorted(du.TRAIN_DONORS + du.TEST_DONORS))
HB_PATTERN       = r"^HB[ABDEGMQZ][0-9]?$"
ERYTHROID_LABELS = ("Erythroid early", "Erythroid late")
T_LABELS         = ("CD4 T naive", "CD4 T memory", "CD8 T memory",
                    "CD8 T cytotoxic", "T activated")
ACTIVATION_GENES = ("CD69", "BTG1", "FOS", "JUN", "JUNB", "DUSP1")
MIN_CELLS        = 20
SCALE_FACTOR     = 1e4
ANNOT_SIZE       = 10


def group_summary(values: pd.DataFrame, obs: pd.DataFrame,
                  keep: pd.Series) -> pd.DataFrame:
    """Per donor and cell type of the kept cells: cell count, mean and share
    of non-zero values per column; groups below the minimum size blanked."""
    kept, keys = values[keep], [obs.loc[keep, du.DONOR_KEY],
                                obs.loc[keep, du.CELL_TYPE_KEY]]
    grouped = kept.groupby(keys, observed=True)
    summary = pd.concat(
        [grouped.mean().add_prefix("mean_"),
         (kept > 0).groupby(keys, observed=True).mean().add_prefix("share_")],
        axis=1,
    )
    summary.insert(0, "cells", grouped.size())
    summary.loc[summary["cells"] < MIN_CELLS, summary.columns[1:]] = np.nan
    return summary


def heatmap(axis: Axes, matrix: pd.DataFrame, vmax: float,
            fmt: str) -> AxesImage:
    """Annotated donor-by-feature heatmap on a fixed 0 to vmax scale."""
    image = axis.imshow(matrix.to_numpy(dtype=float), cmap="magma", vmin=0,
                        vmax=vmax, aspect="auto")
    for (row, col), value in np.ndenumerate(matrix.to_numpy(dtype=float)):
        axis.text(col, row, "n/a" if np.isnan(value) else fmt.format(value),
                  ha="center", va="center", fontsize=ANNOT_SIZE,
                  color="black" if value > 0.6 * vmax else du.FOREGROUND)
    axis.set_yticks(range(matrix.shape[0]), matrix.index)
    axis.set_xticks(range(matrix.shape[1]), matrix.columns, rotation=35,
                    ha="right")
    return image


def plot_haemoglobin(summary: pd.DataFrame) -> Figure:
    """Mean haemoglobin share of counts per donor and non-erythroid cell type,
    as percent, on one fixed colour scale."""
    matrix = 100 * summary["mean_haemoglobin"].unstack()
    figure, axis = plt.subplots(figsize=(16, 7), layout="constrained")
    du.style_dark(figure, [axis])
    image = heatmap(axis, matrix, float(np.nanmax(matrix.to_numpy())),
                    "{:.2f}")
    axis.set_title("Haemoglobin transcripts as percent of counts, "
                   "non-erythroid cells", fontsize=du.TITLE_FONTSIZE,
                   pad=du.TITLE_PAD)
    bar = figure.colorbar(image, ax=axis, shrink=0.85, pad=0.01)
    bar.set_label("Percent of counts", color=du.FOREGROUND,
                  fontsize=du.LABEL_FONTSIZE)
    bar.ax.tick_params(colors=du.FOREGROUND, labelsize=du.TICK_FONTSIZE)
    return figure


def plot_activation(summary: pd.DataFrame) -> Figure:
    """Mean log-normalised immediate-early gene expression per donor, one
    panel per T-cell label stacked vertically on one colour scale."""
    means = summary.filter(like="mean_")
    means.columns = means.columns.str.removeprefix("mean_")
    vmax = float(np.nanmax(means.to_numpy()))
    figure, axes = plt.subplots(len(T_LABELS), 1, layout="constrained",
                                figsize=(12, 4.2 * len(T_LABELS)))
    du.style_dark(figure, list(axes))
    for axis, label in zip(axes, T_LABELS, strict=True):
        matrix = (means.xs(label, level=du.CELL_TYPE_KEY)
                  .reindex(DONORS))
        image = heatmap(axis, matrix, vmax, "{:.2f}")
        axis.set_title(label, fontsize=du.LABEL_FONTSIZE, pad=10)
    figure.suptitle("Immediate-early gene expression in T cells per donor",
                    fontsize=du.TITLE_FONTSIZE, color=du.FOREGROUND)
    bar = figure.colorbar(image, ax=list(axes), shrink=0.5, pad=0.02)
    bar.set_label("Mean log1p counts per 10,000", color=du.FOREGROUND,
                  fontsize=du.LABEL_FONTSIZE)
    bar.ax.tick_params(colors=du.FOREGROUND, labelsize=du.TICK_FONTSIZE)
    return figure


def main() -> None:
    """Both donor checks with their tables and figures."""
    adata = du.load_counts(DONORS)
    obs = adata.obs
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    haemoglobin = adata.var_names.str.match(HB_PATTERN)
    print(f"haemoglobin genes: {', '.join(adata.var_names[haemoglobin])}")
    hb_share = pd.DataFrame(
        {"haemoglobin": np.asarray(adata.X[:, haemoglobin].sum(axis=1))
         .ravel() / totals}, index=obs.index,
    )
    genes = [gene for gene in ACTIVATION_GENES if gene in adata.var_names]
    expression = pd.DataFrame(
        np.log1p(adata[:, genes].X.toarray() / totals[:, None]
                 * SCALE_FACTOR), index=obs.index, columns=genes,
    )

    labels = obs[du.CELL_TYPE_KEY]
    hb_summary = group_summary(hb_share, obs,
                               ~labels.isin(ERYTHROID_LABELS))
    activation = group_summary(expression, obs, labels.isin(T_LABELS))
    du.save_table(hb_summary, SCRIPT, "donor_diagnostics_haemoglobin")
    du.save_table(activation, SCRIPT, "donor_diagnostics_activation")
    du.save_figure(plot_haemoglobin(hb_summary), SCRIPT,
                   "donor_diagnostics_haemoglobin")
    du.save_figure(plot_activation(activation), SCRIPT,
                   "donor_diagnostics_activation")


if __name__ == "__main__":
    main()