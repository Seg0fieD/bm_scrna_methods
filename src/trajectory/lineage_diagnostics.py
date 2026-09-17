"""Diagnostics on the B naive pseudotime tail and on raw count support for selected driver genes."""

import numpy as np
import pandas as pd
from scipy import sparse

from trajectory_utils import RESULT_DIR, load_annotated, save_table


FATE_TAG          = "k8"
PSEUDOTIME_COLUMN = "diffmap_pca_marker"
LINEAGE           = "B naive"
UNEXPECTED_LABEL  = "Pre-B small"
FATE_THRESHOLD    = 0.60
N_BINS            = 40
TAIL_BINS         = (38, 39)
CHECK_GENES       = ("EPCAM", "PTCRA", "HIST1H2BJ", "HIST1H3H", "HIST1H2AC")
QC_CANDIDATES     = (
                  "n_genes_by_counts",
                  "n_genes",
                  "total_counts",
                  "n_counts",
                  "pct_counts_mt",
                  "doublet_score",
                  "scrublet_score",
                  "predicted_doublet",
                  )


def load_pseudotime():
    """Return the Palantir pseudotime that defines the gene trend axis."""
    frame = pd.read_csv(RESULT_DIR / "pt_palantir_pseudotime.csv", index_col=0)
    return frame[PSEUDOTIME_COLUMN]


def load_fates(index, tag):
    """Return the CellRank fate probabilities as a frame of cells by terminal state."""
    matrix = np.load(RESULT_DIR / f"pt_cellrank_fate_{tag}.npy")
    features = pd.read_csv(RESULT_DIR / f"pt_cellrank_fate_{tag}_features.csv")
    return pd.DataFrame(matrix, index=index, columns=features.iloc[:, -1].astype(str))


def committed_bins(pseudotime, probability):
    """Bin the cells committed to the lineage into pseudotime bins of equal cell count."""
    ordered = pseudotime.loc[probability.index[probability >= FATE_THRESHOLD]]
    return pd.Series(pd.qcut(ordered, N_BINS, labels = False), index=ordered.index, name = "bin")


def available_qc(adata):
    """Return the quality control columns present in the object."""
    return [column for column in QC_CANDIDATES if column in adata.obs.columns]


def tail_cells(bins, labels, label):
    """Return the cells carrying the given label within the final pseudotime bins."""
    tail = bins.index[bins.isin(TAIL_BINS)]
    return tail[labels.loc[tail] == label]


def tail_table(adata, bins, pseudotime, fates, labels):
    """Assemble per cell diagnostics for the unexpected label in the final pseudotime bins."""
    cells = tail_cells(bins, labels, UNEXPECTED_LABEL)
    frame = pd.DataFrame(index=cells)
    frame.index.name = "cell"
    frame["bin"] = bins.loc[cells]
    frame["pseudotime"] = pseudotime.loc[cells]
    frame["lineage_fate"] = fates.loc[cells, LINEAGE]
    frame["top_fate"] = fates.loc[cells].idxmax(axis = 1)
    frame["top_fate_probability"] = fates.loc[cells].max(axis = 1)
    for column in available_qc(adata):
        frame[column] = adata.obs.loc[cells, column].values
    return frame.sort_values("pseudotime")


def comparison(adata, bins, pseudotime, fates, labels):
    """Compare the tail cells against their own label and against the whole lineage."""
    groups = {
        f"{UNEXPECTED_LABEL} in tail bins": tail_cells(bins, labels, UNEXPECTED_LABEL),
        f"{LINEAGE} in tail bins": tail_cells(bins, labels, LINEAGE),
        f"all {UNEXPECTED_LABEL}": labels.index[labels == UNEXPECTED_LABEL],
        f"all committed {LINEAGE}": bins.index,
    }
    columns = available_qc(adata)
    records = []
    for name, cells in groups.items():
        record = {
            "group": name,
            "n_cells": int(len(cells)),
            "median_pseudotime": float(pseudotime.loc[cells].median()),
            "median_lineage_fate": float(fates.loc[cells, LINEAGE].median()),
        }
        for column in columns:
            values = pd.to_numeric(adata.obs.loc[cells, column], errors = "coerce")
            record[f"median_{column}"] = float(values.median())
        records.append(record)
    return pd.DataFrame.from_records(records)


def raw_counts(adata):
    """Return the raw count matrix, its gene names, and the object slot it came from."""
    if "counts" in adata.layers:
        return adata.layers["counts"], pd.Index(adata.var_names), "layers['counts']"
    if adata.raw is not None:
        return adata.raw.X, pd.Index(adata.raw.var_names), "raw.X"
    return adata.X, pd.Index(adata.var_names), "X (not raw counts)"


def n_expressing(values):
    """Return the number of cells with a non-zero count."""
    return int((values > 0).sum())


def mean_expressing(values):
    """Return the mean count across the cells that express the gene."""
    positive = values[values > 0]
    return float(positive.mean()) if positive.size else 0.0


def gene_summary(adata, genes):
    """Summarise detection rate and count level of the requested genes in each cell type."""
    matrix, names, source = raw_counts(adata)
    print(f"count source: {source}")
    present = [gene for gene in genes if gene in names]
    for gene in set(genes) - set(present):
        print(f"absent from the matrix: {gene}")

    block = matrix[:, [names.get_loc(gene) for gene in present]]
    dense = block.toarray() if sparse.issparse(block) else np.asarray(block)
    frame = pd.DataFrame(dense, index = adata.obs_names, columns = present)
    frame["cell_type"] = adata.obs["cell_type"].astype(str).values

    long = frame.melt(id_vars = "cell_type", var_name = "gene", value_name = "count")
    grouped = long.groupby(["gene", "cell_type"], observed = True)["count"]
    summary = grouped.agg(
        n_cells         = "size",
        n_expressing    = n_expressing,
        mean_all        = "mean",
        mean_expressing = mean_expressing,
        max_count       = "max",
    ).reset_index()
    summary["pct_expressing"] = summary["n_expressing"] / summary["n_cells"]
    return summary.sort_values(["gene", "pct_expressing"], ascending = [True, False])


def main():
    """Run both diagnostics and write their tables."""
    adata = load_annotated()
    labels = adata.obs["cell_type"].astype(str)
    pseudotime = load_pseudotime().reindex(adata.obs_names)
    fates = load_fates(adata.obs_names, FATE_TAG)
    bins = committed_bins(pseudotime, fates[LINEAGE])

    tail = tail_table(adata, bins, pseudotime, fates, labels)
    summary = comparison(adata, bins, pseudotime, fates, labels)
    save_table(tail.reset_index(), "lineage_diagnostics_tail")
    save_table(summary, "lineage_diagnostics_tail_summary")

    print(f"{UNEXPECTED_LABEL} cells in bins {TAIL_BINS}: {len(tail)}")
    print(tail.to_string())
    print()
    print(summary.to_string(index = False))
    print()

    genes = gene_summary(adata, CHECK_GENES)
    save_table(genes, "lineage_diagnostics_genes")
    print(genes[genes["pct_expressing"] > 0.01].to_string(index=False))


if __name__ == "__main__":
    main()