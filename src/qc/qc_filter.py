"""Filter cells and genes post QC plots visualization"""

from pathlib import Path

import numpy as np
import scanpy as sc

IN = Path("data/interim/bm_qcmetrics.h5ad")
OUT = Path("data/interim/bm_qcfiltered.h5ad")

MT_HARD_LIMIT = 8.0
MIN_CELLS_PER_GENE = 3


def mad_flags(a, metric, n_mads):
    """Return low and high outlier masks for one metric."""
    x = a.obs[metric].to_numpy()
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return x < med - n_mads * mad, x > med + n_mads * mad


def main():
    a = sc.read_h5ad(IN)
    n_start = a.n_obs

    low_counts, _ = mad_flags(a, "log1p_total_counts", 5)
    low_genes, high_genes = mad_flags(a, "log1p_n_genes_by_counts", 5)
    _, high_top20 = mad_flags(a, "pct_counts_in_top_20_genes", 5)
    _, high_mt = mad_flags(a, "pct_counts_mt", 3)

    mt_bad = high_mt & (a.obs["pct_counts_mt"].to_numpy() > MT_HARD_LIMIT)
    drop = low_counts | low_genes | high_top20 | mt_bad

    a.obs["high_complexity"] = high_genes

    for name, mask in [("low counts", low_counts), ("low genes", low_genes),
                       ("high top20", high_top20), ("high mito", mt_bad)]:
        print(f"{name:12s} {int(mask.sum()):6d}")
    print(f"{'total drop':12s} {int(drop.sum()):6d}")
    print("flagged high_complexity, kept:", int((high_genes & ~drop).sum()))

    a = a[~drop].copy()
    sc.pp.filter_genes(a, min_cells=MIN_CELLS_PER_GENE)

    print(f"\ncells {n_start} -> {a.n_obs}")
    print("genes left:", a.n_vars)
    print("\nper donor")
    print(a.obs["donor"].value_counts().sort_index())

    a.write_h5ad(OUT, compression="gzip")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()