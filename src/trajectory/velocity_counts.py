"""
    Spliced, unspliced and ambiguous counts from per-lane velocyto looms,
    merged and restricted to the annotated cells, with per-lane checks.
"""

from __future__ import annotations
import re
from operator import itemgetter
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from trajectory_utils import (ANNOTATED_H5AD, BACKGROUND, FOREGROUND,
                              CELL_TYPE_KEY, FIGSIZE, PROJECT_ROOT,
                              save_figure, save_table, style_dark)


LOOM_DIR        = PROJECT_ROOT / "data" / "interim" / "velocity_loom"
OUTPUT_H5AD     = PROJECT_ROOT / "data" / "interim" / "bm_velocity.h5ad"
SCRIPT_NAME     = Path(__file__).stem

LANE_KEY        = "lane"
LANES           = tuple(f"MantonBM{index}" for index in range(1, 9))
LAYERS          = ("spliced", "unspliced", "ambiguous")
EMBEDDINGS      = ("X_pca", "X_umap")
LOOM_CELL       = re.compile(r"^(?P<lane>[^:]+):(?P<barcode>[ACGT]{16})x$")
FLAG_DEVIATION  = 0.25
TOP_CELL_TYPES  = 3
LABEL_FONTSIZE  = 12
TITLE_FONTSIZE  = 14
MEDIAN_COLOR    = "#f2a541"

LOG_FIELDS      = {
    "reads_skipped_no_barcode": (
        re.compile(r"(\d+) reads were skipped because no apropiate"),
        itemgetter(-1)),
    "reads_counted_in_cells": (
        re.compile(r"containing \d+ cells and (\d+) reads"), sum),
    "reads_repeat_masked": (
        re.compile(r"(\d+) reads not considered because fully enclosed"),
        sum),
    "reads_scanned_millions_floor": (
        re.compile(r"Read first (\d+) million reads"), max),
}


def read_loom(lane: str) -> AnnData:
    """
        Cells-by-genes AnnData of one lane: CSR count layers, cells named in
        the annotated form, genes indexed by Ensembl accession.
    """
    path = LOOM_DIR / f"{lane}.loom"
    with h5py.File(path, "r") as handle:
        layers = {name: sparse.csr_matrix(handle["layers"][name][:, :].T)
                  for name in LAYERS}
        cell_ids = pd.Series(handle["col_attrs"]["CellID"].asstr()[:])
        symbols = handle["row_attrs"]["Gene"].asstr()[:]
        accessions = handle["row_attrs"]["Accession"].asstr()[:]
    parsed = cell_ids.str.extract(LOOM_CELL)
    if parsed.isna().any().any() or not parsed["lane"].eq(lane).all():
        raise ValueError(f"{path.name}: cell names outside the expected form")
    obs_names = pd.Index(parsed["barcode"] + "-1-" + lane)
    if obs_names.has_duplicates:
        raise ValueError(f"{path.name}: repeated cell barcodes")
    print(f"{lane}: {len(obs_names):,} cells x {len(accessions):,} genes")
    return AnnData(X = layers["spliced"].copy(), layers = layers,
                   obs = pd.DataFrame(index = obs_names),
                   var = pd.DataFrame({"gene_symbol": symbols},
                                    index = accessions))


def parse_log(lane: str) -> dict[str, str | float]:
    """
        Read accounting of one velocyto log: barcode-skipped, counted,
        repeat-masked and scanned reads.
    """
    text = (LOOM_DIR / f"{lane}.log").read_text()
    row: dict[str, str | float] = {LANE_KEY: lane}
    for field, (pattern, reduce_values) in LOG_FIELDS.items():
        values = [int(match) for match in pattern.findall(text)]
        row[field] = reduce_values(values) if values else np.nan
    return row


def layer_totals(adata: AnnData) -> pd.DataFrame:
    """
        Per-cell molecule totals per count layer and the unspliced fraction
        of all counted molecules.
    """
    totals = pd.DataFrame(
        {f"n_{name}": np.asarray(adata.layers[name].sum(axis = 1)).ravel()
         for name in LAYERS}, index = adata.obs_names)
    total = totals.sum(axis = 1)
    totals["unspliced_fraction"] = totals["n_unspliced"] / total.where(
        total > 0)
    return totals


def recovery(expected: pd.DataFrame, recovered: pd.DataFrame,
             keys: str | list[str]) -> pd.DataFrame:
    """    
        Expected against recovered cell counts per group, with the recovered
        fraction.
    """
    table = pd.DataFrame({
        "expected": expected.groupby(keys, observed = True).size(),
        "recovered": recovered.groupby(keys, observed = True).size(),
    }).fillna(0).astype(int)
    table["recovered_fraction"] = table["recovered"] / table["expected"]
    return table.reset_index()


def lane_summary(obs: pd.DataFrame,
                 logs: list[dict[str, str | float]]) -> pd.DataFrame:
    """
        Per-lane median counts, unspliced fraction and its deviation from
        the across-lane median, dominant cell types, read accounting.
    """
    grouped = obs.groupby(LANE_KEY, observed = True)
    sums = grouped[[f"n_{name}" for name in LAYERS]].sum()
    summary = grouped[[f"n_{name}" for name in LAYERS]
                      + ["unspliced_fraction"]].median().add_prefix("median_")
    summary["pooled_unspliced_fraction"] = (
        sums["n_unspliced"] / sums.sum(axis = 1))
    summary["deviation_from_lane_median"] = (
        summary["pooled_unspliced_fraction"]
        / summary["pooled_unspliced_fraction"].median() - 1)
    summary["flagged"] = (
        summary["deviation_from_lane_median"].abs() > FLAG_DEVIATION)
    shares = pd.crosstab(obs[LANE_KEY], obs[CELL_TYPE_KEY], normalize = "index")
    summary["dominant_cell_types"] = [
        "; ".join(f"{label} {share:.2f}" for label, share
                  in row.nlargest(TOP_CELL_TYPES).items())
        for _, row in shares.loc[summary.index].iterrows()]
    summary = summary.join(pd.DataFrame(logs).set_index(LANE_KEY))
    summary["reads_masked_fraction_of_counted"] = (
        summary["reads_repeat_masked"] / summary["reads_counted_in_cells"])
    return summary.reset_index()


def lane_figure(obs: pd.DataFrame, summary: pd.DataFrame) -> None:
    """
        Box plot of per-cell unspliced fraction per lane on a fixed axis,
        with the median of the per-lane medians.
    """
    lanes = summary[LANE_KEY].tolist()
    values = [obs.loc[obs[LANE_KEY] == lane, "unspliced_fraction"].dropna()
              for lane in lanes]
    line = {"color": FOREGROUND}
    figure, axis = plt.subplots(figsize = FIGSIZE)
    axis.boxplot(values, showfliers = False, boxprops = line, whiskerprops = line,
                 capprops = line, medianprops={"color": MEDIAN_COLOR})
    axis.axhline(summary["median_unspliced_fraction"].median(),
                 color = FOREGROUND, linestyle = "--", linewidth = 1,
                 label = "median of the\nlane medians")
    axis.set_xticks(range(1, len(lanes) + 1),
                    [f"{lane}\nn = {len(value):,}"
                     for lane, value in zip(lanes, values)])
    axis.set_ylim(0, np.ceil(obs["unspliced_fraction"].quantile(0.995) * 20)
                  / 20)
    axis.set_ylabel("unspliced / all counted molecules per cell",
                    color = FOREGROUND, fontsize = LABEL_FONTSIZE)
    axis.set_title("Unspliced fraction per cell, by sequencing lane",
                   fontsize = TITLE_FONTSIZE)
    axis.tick_params(labelsize = LABEL_FONTSIZE)
    style_dark(figure, axis)
    axis.legend(loc = "upper left", bbox_to_anchor = (1.01, 1), frameon = False,
                fontsize = LABEL_FONTSIZE, labelcolor = FOREGROUND,
                facecolor = BACKGROUND)
    save_figure(figure, f"{SCRIPT_NAME}_lane_unspliced")


def main() -> None:
    """Merged velocity object, recovery tables, lane summary and figure."""
    annotated = ad.read_h5ad(ANNOTATED_H5AD, backed = "r")
    expected = annotated.obs[[CELL_TYPE_KEY]].copy()
    embeddings = {key: np.asarray(annotated.obsm[key]) for key in EMBEDDINGS}
    annotated.file.close()
    expected[LANE_KEY] = expected.index.str.rsplit("-", n = 1).str[1]

    merged = ad.concat([read_loom(lane) for lane in LANES],
                       join = "inner", merge = "same")
    present = expected.index.intersection(merged.obs_names, sort=False)
    velocity = merged[present].copy()
    del merged
    velocity.obs = expected.loc[present].join(layer_totals(velocity))
    positions = expected.index.get_indexer(present)
    for key, matrix in embeddings.items():
        velocity.obsm[key] = matrix[positions]

    by_cell_type = recovery(expected, velocity.obs, CELL_TYPE_KEY)
    by_cell_type["median_unspliced_fraction"] = by_cell_type[
        CELL_TYPE_KEY].map(velocity.obs.groupby(
            CELL_TYPE_KEY, observed = True)["unspliced_fraction"].median())
    by_lane_cell_type = recovery(expected, velocity.obs,
                                 [LANE_KEY, CELL_TYPE_KEY]).join(
        velocity.obs.groupby([LANE_KEY, CELL_TYPE_KEY], observed = True)
        ["unspliced_fraction"].median()
        .rename("median_unspliced_fraction"),
        on = [LANE_KEY, CELL_TYPE_KEY])
    summary = lane_summary(velocity.obs, [parse_log(lane) for lane in LANES])
    tables = {
        "recovery_lane": recovery(expected, velocity.obs, LANE_KEY),
        "recovery_cell_type": by_cell_type,
        "recovery_lane_cell_type": by_lane_cell_type,
        "lane_summary": summary,
    }
    for name, table in tables.items():
        save_table(table, SCRIPT_NAME, f"{SCRIPT_NAME}_{name}", index = False)
    lane_figure(velocity.obs, summary)

    pd.set_option("display.width", 250)
    print(f"\nrecovered {velocity.n_obs:,} of {len(expected):,} cells, "
          f"{velocity.n_vars:,} genes")
    for name in ("recovery_lane", "recovery_cell_type"):
        print(tables[name].round(3).to_string(index = False))
    print(summary.round(3).to_string(index = False))
    activated = velocity.obs.loc[
        velocity.obs[CELL_TYPE_KEY] == "T activated", LANE_KEY]
    print("\nT activated share by lane:\n"
          + activated.value_counts(normalize = True).round(3).to_string())

    velocity.write_h5ad(OUTPUT_H5AD, compression = "gzip")
    print(f"wrote {OUTPUT_H5AD}")


if __name__ == "__main__":
    main()
