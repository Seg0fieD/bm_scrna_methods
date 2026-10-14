"""Training-loss curves of the PyTorch Scaden ensemble: mean loss per pass
over the training data and per-step loss over the final passes."""

from __future__ import annotations

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import StrMethodFormatter

import deconvolution_utils as du

SCRIPT       = "scaden_torch_loss"
FIT_SCRIPT   = "scaden_torch_fit"
TRAIN_FILE   = du.INTERIM_DIR / "scaden_fit" / "training.h5ad"
BATCH_SIZE   = 128
ZOOM_PASSES  = 4
MODEL_COLORS = dict(zip(du.SCADEN_ARCHITECTURES,
                        ("#5aa9e6", "#f4a259", "#7bc96f"), strict=True))


def format_axis(axis: Axes, ylabel: str) -> None:
    """Logarithmic loss axis with plain-number tick labels, grid and label."""
    axis.set_yscale("log")
    axis.yaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
    axis.yaxis.set_minor_formatter(StrMethodFormatter("{x:g}"))
    axis.tick_params(which="both", colors=du.FOREGROUND,
                     labelsize=du.TICK_FONTSIZE)
    axis.set_ylabel(ylabel, fontsize=du.LABEL_FONTSIZE)
    axis.grid(color=du.FOREGROUND, alpha=0.15, which="major")


def add_legend(axis: Axes, colors: dict[str, str], title: str) -> None:
    """Line legend in the reserved column right of an axis."""
    legend = axis.legend(
        handles=[Line2D([], [], color=color, lw=2.5, label=label)
                 for label, color in colors.items()],
        title=title, loc="upper left", bbox_to_anchor=(1.01, 1.0),
        frameon=False, fontsize=du.LEGEND_FONTSIZE,
        title_fontsize=du.LEGEND_FONTSIZE, labelcolor=du.FOREGROUND,
    )
    legend.get_title().set_color(du.FOREGROUND)


def plot_trend(passes: pd.Series) -> Figure:
    """Per-pass mean training loss of every sub-model: mean over seeds as a
    line, range over seeds as a band."""
    spread = passes.groupby(["model", "pass_index"]).agg(["mean", "min",
                                                          "max"])
    figure, axis = plt.subplots(figsize=(15, 8), layout="constrained")
    du.style_dark(figure, [axis])
    for name, color in MODEL_COLORS.items():
        curve = spread.loc[name]
        x = curve.index + 1
        axis.fill_between(x, curve["min"], curve["max"], color=color,
                          alpha=0.3, lw=0)
        axis.plot(x, curve["mean"], lw=2, color=color)
    format_axis(axis, "Mean squared error (mean per pass)")
    axis.set_title("Mean training loss per pass over the training data "
                   "(line: mean over seeds, band: range over seeds)",
                   fontsize=du.TITLE_FONTSIZE, pad=du.TITLE_PAD)
    axis.set_xlabel("Pass over the training data", fontsize=du.LABEL_FONTSIZE)
    add_legend(axis, MODEL_COLORS, "Sub-model")
    return figure


def plot_last_passes(losses: pd.DataFrame, seed: int,
                     steps_per_pass: float) -> Figure:
    """Per-step training loss of every sub-model for one seed over the final
    passes, with pass boundaries marked."""
    last = losses["step"].max()
    start = last + 1 - int(ZOOM_PASSES * steps_per_pass)
    window = losses[(losses["seed"] == seed) & (losses["step"] >= start)]
    figure, axis = plt.subplots(figsize=(15, 8), layout="constrained")
    du.style_dark(figure, [axis])
    for name, color in MODEL_COLORS.items():
        curve = window[window["model"] == name]
        axis.plot(curve["step"], curve["loss"], lw=1.4, color=color)
    for boundary in np.arange(np.ceil(start / steps_per_pass),
                              last / steps_per_pass) * steps_per_pass:
        axis.axvline(boundary, color=du.FOREGROUND, ls=":", lw=1.2)
    format_axis(axis, "Mean squared error (per step)")
    axis.set_title(f"Seed {seed}: training loss per step over the last "
                   f"{ZOOM_PASSES} passes (dotted lines: pass boundaries)",
                   fontsize=du.TITLE_FONTSIZE, pad=du.TITLE_PAD)
    axis.set_xlabel("Training step", fontsize=du.LABEL_FONTSIZE)
    add_legend(axis, MODEL_COLORS, "Sub-model")
    return figure


def main() -> None:
    """Loss figures from the saved per-step training loss: per-pass trend and
    per-step loss over the final passes."""
    losses = pd.read_csv(du.RESULTS_DIR / FIT_SCRIPT
                         / "scaden_torch_fit_loss.csv")
    steps_per_pass = ad.read_h5ad(TRAIN_FILE, backed="r").n_obs / BATCH_SIZE
    print(f"{steps_per_pass:.1f} steps per pass, "
          f"{losses['step'].max() + 1:,} steps")
    pass_index = (losses["step"] // steps_per_pass).rename("pass_index")
    passes = losses.groupby(["model", "seed", pass_index])["loss"].mean()
    du.save_figure(plot_trend(passes), SCRIPT, "scaden_torch_loss_passes")
    du.save_figure(plot_last_passes(losses, losses["seed"].min(),
                                    steps_per_pass),
                   SCRIPT, "scaden_torch_loss_last_passes")


if __name__ == "__main__":
    main()