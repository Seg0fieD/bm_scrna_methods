"""scVI embedding: Variational autoencoder fitted to raw counts with donor as batch, 30 latent dimensions"""

import scvi 

from bench_common import BATCH_KEY, N_COMPS, SEED, load_annotated_data, save_embedding


def main():
    scvi.settings.seed = SEED

    adata = load_annotated_data()
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"{adata.n_obs} cells x {adata.n_vars} variable genes")


    scvi.model.SCVI.setup_anndata(adata, layer = "counts", batch_key = BATCH_KEY)
    model = scvi.model.SCVI(adata, n_latent = N_COMPS)
    model.train(batch_size=128, accelerator="gpu", devices=1)

    save_embedding("scvi", model.get_latent_representation())


if __name__ == "__main__":
    main()