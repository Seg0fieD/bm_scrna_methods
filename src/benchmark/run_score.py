"""Score every embedding with scib-metrics and store in results table."""

import scanpy as sc 
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation

from bench_common import BATCH_KEY, FIG_DIR, LABEL_KEY, TAB_DIR, list_embeddings, load_annotated_data


IN     = TAB_DIR / "bm_benchmark_1.h5ad"
OUT    = TAB_DIR / "scib_scores.csv"
N_JOBS = 4


def main():

    graph = sc.read_h5ad(IN)
    source = load_annotated_data()

    adata = source[graph.obs_names, source.var["highly_variable"]].copy()
    del source
    if "counts" in adata.layers:
        del adata.layers["counts"]

    adata.obs = graph.obs.copy()
    for key, val in graph.obsm.items():
        adata.obsm[key] = val

    keys = [f"X_{name}" for name in list_embeddings()]
    print(f"Scoring {len(keys)} embeddings on {adata.n_obs} cells")

    bench = Benchmarker(
        adata,
        batch_key                = BATCH_KEY,
        label_key                = LABEL_KEY,
        embedding_obsm_keys      = keys, 
        bio_conservation_metrics = BioConservation(),
        batch_correction_metrics = BatchCorrection(pcr_comparison = False),
        n_jobs                   = N_JOBS,
    )
    bench.benchmark()

    results = bench.get_results(min_max_scale = False )
    results.to_csv(OUT)
    bench.get_results(min_max_scale = True).to_csv(TAB_DIR / "scib_scores_scaled.csv")
    print(results)

    bench.plot_results_table(min_max_scale = False, save_dir = str(FIG_DIR))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()


