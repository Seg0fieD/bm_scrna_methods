"""
    Per-lane cell barcode lists for velocyto: the barcodes kept after quality control, 
    split by sequencing lane, in the form the aligned reads carry. 
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from matplotlib.ticker import PercentFormatter

from trajectory_utils import (
    ANNOTATED_H5AD, BACKGROUND, CELL_TYPE_KEY, FOREGROUND, RESULT_DIR,
    save_figure, save_table, style_dark,
)

BARCODE_PATTERN = r"^(?P<barcode>[ACGTN]{16}-\d+)-(?P<lane>[^-]+)$"

FIGURE_SIZE   = (15, 8)
PLOT_RIGHT    = 0.74
BAR_HEIGHT    = 0.72
CATEGORY_CMAP = "tab20"
TITLE_SIZE    = 17
LABEL_SIZE    = 12
TICK_SIZE     = 11
LEGEND_SIZE   = 10


def kept_cells() -> pd.DataFrame:
    """Kept cells as bare barcode, sequencing lane and cell type, parsed from
    ``<barcode>-<gem group>-<lane>`` observation names; any name outside that
    form, or a barcode repeated within a lane, is rejected."""
    adata = sc.read_h5ad(ANNOTATED_H5AD, backed = "r")
    obs = adata.obs[[CELL_TYPE_KEY]].astype(str)
    adata.file.close()
    cells = obs.index.to_series().str.extract(BARCODE_PATTERN)
    if cells.isna().any(axis = None):
        bad = cells.index[cells.isna().any(axis = 1)][:3].tolist()
        raise ValueError(f"observation names outside the expected form: {bad}")
    if cells.duplicated().any():
        raise ValueError("a barcode occurs twice within one lane")
    return cells.join(obs)


def write_barcode_lists(cells: pd.DataFrame) -> pd.DataFrame:
    """One headerless file per lane, one barcode per line, as read by
    ``velocyto run -b``; the per-lane summary with its dominant cell type."""
    rows = []
    for lane, group in cells.groupby("lane"):
        path = RESULT_DIR / f"velocity_barcodes_{lane}.tsv"
        group["barcode"].to_csv(path, index = False, header = False)
        shares = group[CELL_TYPE_KEY].value_counts(normalize = True)
        rows.append(
            {
                "lane"           : lane,
                "cells"          : len(group),
                "dominant_type"  : shares.index[0],
                "dominant_share" : round(shares.iloc[0], 3),
                "file"           : path.name,
            }
        )
        print(f"wrote {path} ({len(group)} barcodes)")
    return pd.DataFrame(rows)


def plot_composition(cells: pd.DataFrame) -> None:
    """Cell type composition of every lane as stacked shares, cell types in
    order of overall abundance, each bar labelled with its kept-cell total."""
    counts = pd.crosstab(cells["lane"], cells[CELL_TYPE_KEY])
    counts = counts[counts.sum().sort_values(ascending = False).index]
    shares = counts.div(counts.sum(axis = 1), axis = 0)

    figure, axis = plt.subplots(figsize = FIGURE_SIZE)
    figure.subplots_adjust(left = 0.09, right = PLOT_RIGHT)
    shares.plot.barh(
        stacked   = True,
        ax        = axis,
        width     = BAR_HEIGHT,
        color     = plt.get_cmap(CATEGORY_CMAP).colors[: shares.shape[1]], # type: ignore
        edgecolor = BACKGROUND,
        linewidth = 0.4,
        legend    = False,
    )
    for row, total in enumerate(counts.sum(axis = 1)):
        axis.text(
            1.01, row, f"n = {total:,}", va = "center", fontsize = TICK_SIZE,
            color = FOREGROUND,
        )
    axis.set_xlim(0.0, 1.0)
    axis.invert_yaxis()
    axis.xaxis.set_major_formatter(PercentFormatter(1.0))
    axis.tick_params(labelsize = TICK_SIZE)
    axis.set_xlabel("share of the lane's kept cells", fontsize = LABEL_SIZE)
    axis.set_ylabel("lane", fontsize = LABEL_SIZE)
    axis.set_title(
        "Kept cells per lane, by cell type", fontsize = TITLE_SIZE, pad = 14
    )
    style_dark(figure, axis)
    axis.xaxis.label.set_color(FOREGROUND)
    axis.yaxis.label.set_color(FOREGROUND)

    legend = figure.legend(
        *axis.get_legend_handles_labels(),
        title          = "Cell type",
        loc            = "upper left",
        bbox_to_anchor = (PLOT_RIGHT + 0.08, 0.92),
        frameon        = False,
        fontsize       = LEGEND_SIZE,
        title_fontsize = LEGEND_SIZE + 2,
        labelspacing   = 0.45,
    )
    legend.get_title().set_color(FOREGROUND)
    for text in legend.get_texts():
        text.set_color(FOREGROUND)
    save_figure(figure, "velocity_barcodes_composition")


def main() -> None:
    """Per-lane barcode files, their summary table and the lane composition
    figure."""
    cells = kept_cells()
    summary = write_barcode_lists(cells)
    save_table(summary, "velocity_barcodes_summary", index = False)
    plot_composition(cells)
    print(summary.to_string(index = False))


if __name__ == "__main__":
    main()