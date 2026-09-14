"""Consensus NMF embedding: gene programs fitted to raw counts, cell usage matrix as embeddings"""

from cnmf import cNMF

from bench_common import N_COMPS, SEED, TAB_DIR, load_annotated_data, save_embedding

N_ITER = 20 
DENSITY_THRESHOLD = 2.00


def main():
    adata = load_annotated_data()
    adata = adata[:, adata.var["highly_variable"]].copy()
    adata.X = adata.layers["counts"].copy()
    del adata.layers["counts"]
    print(f"{adata.n_obs} cells x {adata.n_vars} variable genes")

    cnmf_dir = TAB_DIR / "cnmf"
    cnmf_dir.mkdir(parents = True, exist_ok = True)
    counts_fn = cnmf_dir / "cnmf_counts.h5ad"
    adata.write_h5ad(counts_fn)

    run = cNMF(output_dir= str(cnmf_dir), name = "bm")
    run.prepare(
        counts_fn         = str(counts_fn),
        components        = [N_COMPS],
        n_iter            = N_ITER,
        seed              = SEED, 
        num_highvar_genes = adata.n_vars,
    )
    run.factorize(worker_i = 0, total_workers = 1)
    run.combine()
    run.consensus(k = N_COMPS, density_threshold = DENSITY_THRESHOLD, 
                  show_clustering = True)

    usage, spectra_scores, spectra_tpm, top_genes = run.load_results(
        K = N_COMPS, density_threshold = DENSITY_THRESHOLD)

    usage = usage.loc[adata.obs_names]

    spectra_scores.to_csv(TAB_DIR / "cnmf_spectra_scores.csv")
    spectra_tpm.to_csv(TAB_DIR / "cnmf_spectra_tpm.csv")

    save_embedding("cnmf", usage.to_numpy())
    top_genes.to_csv(TAB_DIR / "cnmf_top_genes.csv")


if __name__ == "__main__":
    main()
    
    