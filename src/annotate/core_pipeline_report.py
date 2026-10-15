"""Summary tables of the core pipeline: cells per step and donor, quality
control rules, doublet cutoffs, variable genes and cells per label."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT    = Path(__file__).resolve().parents[2]
INTERIM = ROOT / "data" / "interim"
SCRIPT  = Path(__file__).stem
TAB     = ROOT / "results" / "annotate" / SCRIPT

STEPS = {
    "merged": INTERIM / "bm_merged.h5ad",
    "qc_filtered": INTERIM / "bm_qcfiltered.h5ad",
    "scrublet_removed": INTERIM / "bm_scored.h5ad",
    "top4_removed": INTERIM / "bm_clean.h5ad",
    "annotated": ROOT / "data" / "processed" / "bm_annotated.h5ad",
}
QC_METRICS = INTERIM / "bm_qcmetrics.h5ad"
HVG        = INTERIM / "bm_normalized_hvg.h5ad"

QC_RULES = pd.DataFrame(
    [
        ("low counts", "log1p_total_counts", 5, "low", -np.inf, True),
        ("low genes", "log1p_n_genes_by_counts", 5, "low", -np.inf, True),
        ("high top-20 share", "pct_counts_in_top_20_genes", 5, "high",
         -np.inf, True),
        ("high mitochondrial", "pct_counts_mt", 3, "high", 8.0, True),
        ("high complexity", "log1p_n_genes_by_counts", 5, "high", -np.inf,
         False),
    ],
    columns=["rule", "metric", "n_mads", "side", "floor", "removes"],
)
DOUBLET_RATE   = 0.04
MEDIAN_COLUMNS = ["n_genes_by_counts", "total_counts", "pct_counts_mt",
                  "pct_counts_ribo", "pct_counts_hb",
                  "pct_counts_in_top_20_genes"]
HVG_COLUMNS    = ["highly_variable_rank", "highly_variable_nbatches", "means",
                  "variances_norm", "mt", "ribo", "hb"]


def metadata(path: Path) -> ad.AnnData:
    """Cell and gene annotation of a saved data file, without the matrix."""
    backed = ad.read_h5ad(path, backed="r")
    meta = ad.AnnData(obs=backed.obs.copy(), var=backed.var.copy(),
                      uns=dict(backed.uns))
    backed.file.close()
    return meta


def qc_rule_table(obs: pd.DataFrame, kept: pd.Index) -> pd.DataFrame:
    """Median-absolute-deviation bound and flagged cell count per quality
    control rule, checked against the cells the filter kept."""
    removed = np.zeros(len(obs), dtype=bool)
    bounds, flagged_counts = [], []
    for rule in QC_RULES.itertuples():
        values = obs[rule.metric].to_numpy()
        median = np.median(values)
        spread = rule.n_mads * np.median(np.abs(values - median))
        if rule.side == "low":
            bound = median - spread
            flagged = values < bound
        else:
            bound = median + spread
            flagged = values > bound
        flagged &= values > rule.floor
        bounds.append(bound)
        flagged_counts.append(int(flagged.sum()))
        if rule.removes:
            removed |= flagged
    if not obs.index[~removed].equals(kept):
        raise RuntimeError("recomputed QC rules do not match the kept cells")
    return QC_RULES.assign(bound=bounds, cells_flagged=flagged_counts)


def doublet_table(meta: dict[str, ad.AnnData]) -> pd.DataFrame:
    """Scrublet threshold and calls, and the top-score cutoff and removals,
    per donor; checked against the cells each step kept."""
    scored = meta["scrublet_removed"]
    counts = {step: meta[step].obs["donor"].value_counts()
              for step in ("qc_filtered", "scrublet_removed", "top4_removed")}
    scores = scored.obs["doublet_score"]
    cutoff = scores.groupby(scored.obs["donor"], observed=True).quantile(
        1 - DOUBLET_RATE)
    above = scores > scored.obs["donor"].map(cutoff).astype(float)
    table = pd.DataFrame({
        "cells_scored": counts["qc_filtered"],
        "scrublet_threshold": pd.Series({
            donor: float(batch["threshold"])
            for donor, batch in scored.uns["scrublet"]["batches"].items()}),
        "scrublet_calls": counts["qc_filtered"] - counts["scrublet_removed"],
        "top4_cutoff": cutoff,
        "top4_removed": (counts["scrublet_removed"]
                         - counts["top4_removed"]),
    }).rename_axis("donor").sort_index()
    if not above.groupby(scored.obs["donor"], observed=True).sum().equals(
            table["top4_removed"]):
        raise RuntimeError("recomputed doublet cutoff does not match")
    return table


def save_table(table: pd.DataFrame, name: str) -> None:
    """Results table as CSV in the script's results directory; logs the
    output path."""
    TAB.mkdir(parents=True, exist_ok=True)
    path = TAB / f"{SCRIPT}_{name}.csv"
    table.to_csv(path)
    print("wrote", path)


def main() -> None:
    """Core pipeline summary tables, read from the saved data files."""
    meta = {step: metadata(path) for step, path in STEPS.items()}
    qc = metadata(QC_METRICS)
    hvg = metadata(HVG).var

    steps = pd.DataFrame(
        {"cells": [m.n_obs for m in meta.values()],
         "genes": [m.n_vars for m in meta.values()]},
        index=pd.Index(list(meta), name="step"),
    )
    donors = pd.DataFrame(
        {step: m.obs["donor"].value_counts() for step, m in meta.items()}
    ).rename_axis("donor").sort_index()
    labels = meta["annotated"].obs
    cell_types = pd.crosstab(labels["cell_type"], labels["donor"])

    save_table(steps, "steps")
    save_table(donors, "donors")
    save_table(qc.obs.groupby("donor", observed=True)[MEDIAN_COLUMNS]
               .median(), "qc_medians")
    save_table(qc_rule_table(qc.obs, meta["qc_filtered"].obs_names)
               .set_index("rule"), "qc_rules")
    save_table(doublet_table(meta), "doublets")
    save_table(hvg.loc[hvg["highly_variable"], HVG_COLUMNS]
               .sort_values("highly_variable_rank"), "hvg")
    save_table(cell_types.assign(total=cell_types.sum(axis=1)), "cell_types")


if __name__ == "__main__":
    main()