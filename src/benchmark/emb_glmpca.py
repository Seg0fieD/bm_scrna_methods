"""GLM-PCA embedding: Poisson factor model on raw counts, 30 components"""
import numpy as np
from glmpca import glmpca
from scipy.sparse import issparse

from benchmark_utils import N_COMPS,  load_annotated_data, save_embedding # SEED, EMB_DIR,

# N_CELLS = 2000

def main() -> None:
    """GLM-PCA embedding saved for the benchmark."""
    adata = load_annotated_data()
    adata = adata[:, adata.var["highly_variable"]].copy()

    # rng = np.random.default_rng(SEED)
    # keep = np.sort(rng.choice(adata.n_obs, size = N_CELLS, replace = False))
    # adata = adata[keep].copy()

    raw = adata.layers["counts"]
    counts = raw.toarray() if issparse(raw) else np.asarray(raw)

    Y = counts.T
    Y = Y[Y.sum(axis = 1) > 0]
    print(f"{Y.shape[1]} cells x {Y.shape[0]} genes, fitting {N_COMPS} factors")


    res = glmpca.glmpca(Y, N_COMPS, fam="poi", penalty=10, verbose=True)
    save_embedding("glmpca", res["factors"])
    # np.save(EMB_DIR / "glmpca_cells.npy", adata.obs_names.to_numpy())


if __name__ == "__main__":
    main()




