"""Cross-check the manual labels with CellTypist."""

from pathlib import Path

import celltypist
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from celltypist import models

IN = Path("data/processed/bm_annotated.h5ad")
OUT = Path("data/processed/bm_annotated.h5ad")
FIG = Path("figures/annotate")
TAB = Path("results/tables")

MODEL = "Immune_All_Low.pkl"
GROUP = "leiden_1.0"


def main():
    a = sc.read_h5ad(IN)

    models.download_models(model=[MODEL])
    result = celltypist.annotate(
        a, model=MODEL, majority_voting=True, over_clustering=GROUP
    )

    a.obs["celltypist"] = result.predicted_labels["majority_voting"].values

    table = pd.crosstab(a.obs["cell_type"], a.obs["celltypist"])
    TAB.mkdir(parents=True, exist_ok=True)
    table.to_csv(TAB / "celltypist_vs_manual.csv")

    print("\nbest CellTypist match per manual label")
    for label in a.obs["cell_type"].cat.categories:
        row = table.loc[label]
        best = row.idxmax()
        share = 100 * row[best] / row.sum()
        print(f"  {label:24s} -> {best}  ({share:.0f}%)")

    FIG.mkdir(parents = True, exist_ok = True)
    sc.pl.umap(a, color = "celltypist", legend_loc = "on data",
               legend_fontsize = 6, show = False)
    plt.gcf().savefig(FIG / "umap_celltypist.png", dpi=300,
                      bbox_inches = "tight")
    plt.close("all")

    a.write_h5ad(OUT, compression = "gzip")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()