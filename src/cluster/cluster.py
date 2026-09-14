"""Build the neighbour graph, run UMAP and cluster at three resolutions."""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

IN = Path("data/interim/bm_pca.h5ad")
OUT = Path("data/interim/bm_clustered.h5ad")
FIG = Path("figures/cluster")

N_PCS = 30
N_NEIGHBORS = 15 
RESOLUTIONS = [0.4, 0.6, 0.8, 1.0, 1.2, 1.5]
QC_KEYS = ["donor", "doublet_score", "pct_counts_mt",
           "n_genes_by_counts", "HBB", "MPO"]

def save_umap(a, key, name, legend_loc=None):
    """Create UMAP coloured by key and save it."""

    sc.pl.umap(a, color = key, legend_loc = legend_loc, show = False)
    plt.gcf().savefig(FIG / f"{name}.png", dpi = 250, bbox_inches = "tight")
    plt.close("all")

def main():
    a = sc.read_h5ad(IN)

    sc.pp.neighbors(a, n_neighbors = N_NEIGHBORS, n_pcs = N_PCS)
    sc.tl.umap(a)

    keys = []
    for res in RESOLUTIONS:
        key = f"leiden_{res}"
        sc.tl.leiden(
            a, resolution = res, key_added = key,
            flavor = "igraph", n_iterations = 2, directed = False,
        )
        print(key, "clusters:", a.obs[key].nunique())
        keys.append(key)

    FIG.mkdir(parents = True, exist_ok = True)

    for key in keys:
        save_umap(a, key, f"umap_{key}", legend_loc = "on data")

    for key in QC_KEYS:
        save_umap(a, key, f"umap_{key}")

    sc.pl.umap(a, color = keys, ncols = 3, legend_loc = "on data", show = False)
    plt.gcf().savefig(FIG / "umap_leiden_grid.png", dpi = 250,
                      bbox_inches = "tight")
    plt.close("all")

    a.write_h5ad(OUT, compression = "gzip")
    print("wrote", OUT)


if __name__ == "__main__":
    main()