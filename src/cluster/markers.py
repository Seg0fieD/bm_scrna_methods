"""Find marker genes per cluster and check which clusters look alike."""

import warnings
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
from pandas.errors import PerformanceWarning

warnings.filterwarnings("ignore", category=PerformanceWarning)

IN = Path("data/interim/bm_clustered.h5ad")
OUT = Path("data/interim/bm_markers.h5ad")
FIG = Path("figures/annotate")
TAB = Path("results/tables")

GROUP = "leiden_1.0"
N_TOP = 25
N_COMPARE = 20
SHARED_LIMIT = 10
N_DOT = 3


def cluster_summary(a):
    """Size and mean QC values per cluster."""
    cols = ["doublet_score", "pct_counts_mt", "n_genes_by_counts",
            "pct_counts_hb", "pct_counts_ribo"]
    out = a.obs.groupby(GROUP, observed=True)[cols].mean().round(2)
    out.insert(0, "cells", a.obs[GROUP].value_counts().sort_index())
    return out


def clean_markers(a):
    """Marker table without ribosomal, mitochondrial and haemoglobin genes."""
    table = sc.get.rank_genes_groups_df(a, group=None)
    drop = a.var_names[a.var["ribo"] | a.var["mt"] | a.var["hb"]]
    return table[~table["names"].isin(drop)]


def top_genes(table, n):
    """Dict of cluster to its top n genes, in order."""
    out = {}
    for cluster in sorted(table["group"].unique(), key=int):
        rows = table[table["group"] == cluster]
        out[cluster] = list(rows["names"][:n])
    return out


def find_similar(tops):
    """Print cluster pairs that share many of their top genes."""
    print(f"\ncluster pairs sharing >= {SHARED_LIMIT} of top {N_COMPARE}")
    found = False
    for one, two in combinations(tops, 2):
        shared = set(tops[one]) & set(tops[two])
        if len(shared) >= SHARED_LIMIT:
            print(f"  {one} & {two}: {len(shared)} -> {sorted(shared)}")
            found = True
    if not found:
        print("  none")


def main():
    a = sc.read_h5ad(IN)
    print("clusters:", a.obs[GROUP].nunique())

    summary = cluster_summary(a)
    print("\nper cluster")
    print(summary.to_string())

    sc.tl.rank_genes_groups(a, groupby=GROUP, method="wilcoxon")
    sc.tl.dendrogram(a, groupby=GROUP)

    TAB.mkdir(parents=True, exist_ok=True)
    sc.get.rank_genes_groups_df(a, group=None).to_csv(
        TAB / f"markers_{GROUP}_all.csv", index=False
    )

    table = clean_markers(a)
    table.to_csv(TAB / f"markers_{GROUP}_clean.csv", index=False)
    summary.to_csv(TAB / f"cluster_summary_{GROUP}.csv")

    tops = top_genes(table, N_TOP)
    print("\ntop 10 per cluster, housekeeping genes removed")
    for cluster in tops:
        print(f"  {cluster}: {', '.join(tops[cluster][:10])}")

    find_similar(top_genes(table, N_COMPARE))

    FIG.mkdir(parents=True, exist_ok=True)
    dot = {}
    for cluster in tops:
        dot[cluster] = tops[cluster][:N_DOT]

    sc.pl.dotplot(a, dot, groupby=GROUP, dendrogram=True,
                  standard_scale="var", show=False)
    plt.gcf().savefig(FIG / "dotplot_top3.png", dpi = 250, bbox_inches="tight")
    plt.close("all")

    a.write_h5ad(OUT, compression="gzip")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()