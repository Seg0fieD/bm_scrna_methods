"""Compute QC metrics on the merged data and plot them per donor."""

from pathlib import Path

import numpy as np
import pandas as pd 
import scanpy as sc
import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings("ignore", message = "Variable names are not unique")

IN = Path("data/interim/bm_merged.h5ad")
OUT = Path("data/interim/bm_qcmetrics.h5ad")
FIG = Path("figures/qc")

HB_PREFIXES = ("HBA", "HBB", "HBD", "HBG", "HBM", "HBQ", "HBZ")

def flag_gene_groups(a):
    """Mark mitochondrial, ribosomal and haemoglobin genes in var."""

    a.var["mt"]   = a.var_names.str.startswith("MT-")
    a.var["ribo"] = a.var_names.str.startswith(("RPS", "RPL"))
    a.var["hb"]   = a.var_names.str.startswith(HB_PREFIXES)

def mad_bounds(values, n_mads = 3):
    """Return lower and upper cutoffs at n_mads median absolute deviations."""

    med = np.median(values)
    mad = np.median(np.abs(values - med))
    return med - n_mads * mad, med + n_mads * mad

def report(a):
    """Print per-donor medians and suggested MAD cutoffs."""

    cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt",
            "pct_counts_ribo", "pct_counts_hb", "pct_counts_in_top_20_genes"]

    print("\nper-donor medians")
    print(a.obs.groupby("donor", observed = True)[cols].median().round(2))

    print("\nsuggested cutoff (3 MAD, whole dataset)")
    rows = []
    for col in ["log1p_n_genes_by_counts", "log1p_total_counts",
                "pct_counts_in_top_20_genes", "pct_counts_mt"]:
        low, high = mad_bounds(a.obs[col].to_numpy())
        n_out = int(((a.obs[col] < low) | (a.obs[col] > high)).sum())
        rows.append({"metric": col, "low": round(low, 2),
                     "high": round(high, 2), "cells_outside": n_out})
    print(pd.DataFrame(rows).to_string(index = False))

def plot(a):
    """Save violin and scatter plots of the QC metrics."""

    axes = sc.pl.violin(
                a , 
                ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb"],
                groupby = "donor", rotation = 90, stripplot = False, show = False,
            )
    
    axes[0].figure.savefig(
        FIG / "violin_by_donor.png", dpi = 150, bbox_inches = "tight"
    )

    ax = sc.pl.scatter(
            a, x = "total_counts", y = "n_genes_by_counts",
            color = "pct_counts_mt",
            show = False
        )
    ax.figure.savefig(
        FIG / "scatter_counts_vs_genes.png", dpi = 300, bbox_inches="tight"
    )

    plt.close("all")




def main():
    a = sc.read_h5ad(IN)
    flag_gene_groups(a)
    print("flagged genes:", a.var[["mt", "ribo", "hb"]].sum().to_dict())

    sc.pp.calculate_qc_metrics(
        a, 
        qc_vars = ["mt", "ribo", "hb"], percent_top = [20],
        log1p = True, inplace = True,
    )

    report(a)

    FIG.mkdir(parents = True, exist_ok = True)
    # sc.settings.figdir = FIG
    # sc.settings.autoshow = False
    plot(a)

    a.write_h5ad(OUT, compression = "gzip")
    print("\nwrote ", OUT)


if __name__ == "__main__":
    main()
