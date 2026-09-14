"""PCA baseline embeddings: log-normalized variable genes, no scaling, 30 Components"""

import scanpy as sc 

from bench_common import N_COMPS, SEED, load_annotated_data, save_embedding

def main():
    adata = load_annotated_data()
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"{adata.n_obs} cells x {adata.n_vars} variable genes")

    sc.pp.pca(adata, n_comps = N_COMPS, svd_solver = "arpack", random_state = SEED)
    save_embedding("pca", adata.obsm["X_pca"])



if __name__ == "__main__":
    main()
