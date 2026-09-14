"""Score doublets with Scrublet per donor and remove the predicted ones."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

IN = Path("data/interim/bm_qcfiltered.h5ad")
OUT = Path("data/interim/bm_clean.h5ad")
FIG = Path("figures/qc")


def main():
    a = sc.read_h5ad(IN)
    n_start = a.n_obs

    sc.pp.scrublet(a, batch_key="donor", random_state=0)

    print("\npredicted doublets per donor")
    print(pd.crosstab(a.obs["donor"], a.obs["predicted_doublet"]))

    hc = a.obs["high_complexity"].to_numpy()
    dbl = a.obs["predicted_doublet"].to_numpy()
    print("\nhigh_complexity cells:", int(hc.sum()))
    print("of those called doublet:", int((hc & dbl).sum()))

    FIG.mkdir(parents=True, exist_ok=True)
    ax = sc.pl.violin(
        a, "doublet_score", groupby="donor", rotation=90,
        stripplot=False, show=False,
    )
    ax.figure.savefig(
        FIG / "violin_doublet_score.png", dpi=150, bbox_inches="tight"
    )
    plt.close("all")

    a = a[~a.obs["predicted_doublet"]].copy() # keeps the cells scrublet did not call doublets
    sc.pp.filter_genes(a, min_cells=3)

    print(f"\ncells {n_start} -> {a.n_obs}")
    print("genes left:", a.n_vars)

    a.write_h5ad(OUT, compression="gzip")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()