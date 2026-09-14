"""Diffusion map embeddings build on top of the PCA and scVI embeddings"""

import anndata as ad
import numpy as np
import scanpy as sc 

from bench_common import N_COMPS, N_NEIGHBORS, SEED, load_embedding, save_embedding

BASES = ("pca", "scvi")

def main():
    for base in BASES:
        X = load_embedding(base)
        adata = ad.AnnData(np.zeros((X.shape[0], 1), dtype = np.float32))
        adata.obsm["X_base"] = X

        sc.pp.neighbors(adata, n_neighbors = N_NEIGHBORS, 
                        use_rep = "X_base", random_state = SEED)
        sc.tl.diffmap(adata, n_comps = N_COMPS + 1)

        save_embedding(f"diffmap_{base}", adata.obsm["X_diffmap"][:, 1:])

if __name__ == "__main__":
    main()
    
