"""Pearson residual embedding: analytic Pearson residual on variable genes, 30 Components"""

import scanpy as sc 

from bench_common import N_COMPS, SEED, load_annotated_data, save_embedding

def main():
    adata = load_annotated_data()    
    adata = adata[:, adata.var ["highly_variable"]].copy()
    adata.X = adata.layers["counts"].copy()
    print(f"{adata.n_obs} cells x {adata.n_vars} variable genes")

    sc.experimental.pp.normalize_pearson_residuals(adata)
    sc.pp.pca(adata, n_comps = N_COMPS, svd_solver = "arpack", 
              random_state = SEED)
    save_embedding("pearson", adata.obsm["X_pca"])


if __name__ == "__main__":
    main()
    