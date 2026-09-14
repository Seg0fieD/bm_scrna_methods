"""Normalize, log-transform and pick highly variable genes."""

from pathlib import Path

import scanpy as sc

IN  = Path("data/interim/bm_clean.h5ad") 
OUT = Path("data/interim/bm_normalized_hvg.h5ad")

N_HVG = 2000

def main():
    a = sc.read_h5ad(IN)
    a.layers["counts"] = a.X.copy()

    sc.pp.normalize_total(a, target_sum = 1e4)
    sc.pp.log1p(a)

    sc.pp.highly_variable_genes(
        a, n_top_genes = N_HVG, flavor ="seurat_v3",
        layer = "counts", batch_key = "donor",
    )

    hvg = a.var["highly_variable"]
    print(f"HVGs:  {int(hvg.sum())}")

    for group in ["mt", "ribo", "hb"]:
        print(f" {group}: {int((hvg & a.var[group]).sum())}")

    print("\ntop 20 by rank")
    print(list(a.var.loc[hvg].sort_values("highly_variable_rank").index[:20]))

    a.write_h5ad(OUT, compression = "gzip")
    print("\nwrote", OUT)



if __name__ == "__main__":
    main()
