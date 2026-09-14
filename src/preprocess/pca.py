"""Run PCA on the HVGs and check how strong the donor batch effect is."""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc 

IN = Path("data/interim/bm_normalized_hvg.h5ad")
OUT = Path("data/interim/bm_pca.h5ad")
FIG = Path("figures/pca")


N_COMPS = 50

def main():
    a = sc.read_h5ad(IN)

    sc.pp.pca(a, n_comps = N_COMPS, mask_var = "highly_variable",
              svd_solver = "arpack")

    FIG.mkdir(parents = True, exist_ok = True)

    vr = a.uns["pca"]["variance_ratio"]
    fig, ax = plt.subplots(figsize = (10, 6))
    ax.plot(range(1, len(vr) + 1), vr, marker = "o", markersize = 3)
    ax.set_yscale("log")
    ax.set_xlabel("PC")
    ax.set_ylabel("Variance ratio")
    fig.savefig(FIG / "variance_ratio.png", dpi = 150, bbox_inches = "tight")

    sc.pl.pca(
        a, color=["donor", "pct_counts_mt", "doublet_score", "HBB"],
        ncols=2, show=False,
    )
    plt.gcf().savefig(FIG / "pca_overview.png", dpi=250, bbox_inches="tight")

    plt.close("all")


    a.write_h5ad(OUT, compression = "gzip")
    print("wrote ", OUT)


if __name__ == "__main__":
    main()

