"""Call doublets by per-donor score rank instead of Scrublet's auto threshold"""

from pathlib import Path 

import scanpy as sc 

IN = Path("data/interim/bm_scored.h5ad")
OUT = Path("data/interim/bm_clean.h5ad")


DOUBLET_RATE = 0.04

def main():
    a       = sc.read_h5ad(IN)
    n_start = a.n_obs

    cutoff = a.obs.groupby("donor", observed = True)["doublet_score"].quantile(
        1 - DOUBLET_RATE
    )

    print("per-donor cutoff")
    print(cutoff.round(3))

    limit = a.obs["donor"].map(cutoff).to_numpy()
    drop = a.obs["doublet_score"].to_numpy() > limit

    print("\ndropped per donor")
    print(a.obs.loc[drop, "donor"].value_counts().sort_index())

    a = a[~drop].copy()
    sc.pp.filter_genes(a, min_cells = 3)

    print(f"\ncells {n_start} -> {a.n_obs}")
    print("genes left: ", a.n_vars)

    a.write_h5ad(OUT, compression = "gzip")
    print("\nwrote", OUT)

if __name__ == "__main__":
    main()
