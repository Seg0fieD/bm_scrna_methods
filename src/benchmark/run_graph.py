"""Build a neighbour graph, UMAP, and clustering for every saved embedding, using identical settings"""

import anndata as ad 
import numpy as np
import scanpy as sc

from bench_common import ( BATCH_KEY, LABEL_KEY, N_NEIGHBORS, RESOLUTION, SEED, 
                          TAB_DIR, list_embeddings, load_annotated_data, 
                          load_embedding )

OUT = TAB_DIR / "bm_benchmark_1.h5ad"

def main():
    source = load_annotated_data()
    adata = ad.AnnData(
        np.zeros((source.n_obs, 1), dtype = np.float32),
        obs = source.obs[[LABEL_KEY, BATCH_KEY]].copy(),
    )

    adata.obs_names = source.obs_names
    del source

    for name in list_embeddings():
        rep = f"X_{name}"
        adata.obsm[rep] = load_embedding(name)

        sc.pp.neighbors(
            adata, n_neighbors = N_NEIGHBORS, use_rep = rep, 
            key_added = name, random_state = SEED
        )
        sc.tl.umap(adata, neighbors_key = name, random_state = SEED)
        adata.obsm[f"X_umap_{name}"] = adata.obsm["X_umap"]

        sc.tl.leiden(
            adata,
            resolution    = RESOLUTION,
            key_added     = f"leiden_{name}",
            neighbors_key = name,
            flavor        = "igraph",
            n_iterations  = 2, 
            directed      = False,
            random_state  = SEED,
        )
        print(f"{name} : {adata.obs[f'leiden_{name}'].nunique()} clusters")

    del adata.obsm["X_umap"]
    adata.write_h5ad(OUT)
    print(f"wrote {OUT}")

if __name__ == "__main__":
    main()       