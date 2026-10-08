"""Build the comparison figure and ranked score table for the dimensionality-reduction benchmark"""

import matplotlib.pyplot as plt 
import pandas as pd
import scanpy as sc
from anndata import AnnData

from benchmark_utils import ( BATCH_KEY, FIG_DIR, GRAPH_H5AD, LABEL_KEY, SCORE_DIR, LEGEND_FONTSIZE,
                          list_embeddings, save_figure )


IN     = GRAPH_H5AD
SCORES = SCORE_DIR / "scib_scores.csv"
NCOLS  = 4

def umap_grid(adata: AnnData, names: list[str], color: str, filename: str) -> None:
    """One UMAP panel per embedding, all coloured by the same column."""
    nrows = -(-len(names) // NCOLS)
    fig, axes = plt.subplots(nrows = nrows, ncols = NCOLS, figsize = (22, 5 * nrows))
    axes = axes.ravel()


    for idx, name in enumerate(names):
        last = idx == len(names) - 1
        adata.obsm["X_umap"] = adata.obsm[f"X_umap_{name}"]
        sc.pl.umap(
            adata, color = color,
            ax = axes[idx], show = False, 
            title = name,  frameon = False,
            size = 3, legend_loc = "right margin" if last else None,
            legend_fontsize = LEGEND_FONTSIZE,
        ) 
    for ax in axes[len(names):]:
        ax.axis("off")

    del adata.obsm["X_umap"]
    save_figure(FIG_DIR / filename, fig)

def ranked_scores() -> pd.DataFrame:
    """Score table ranked by total score, the metric-type row dropped."""
    scores = pd.read_csv(SCORES, index_col = 0)
    scores = scores.drop(index = "Metric Type", errors = "ignore").astype(float)
    scores = scores.sort_values("Total", ascending = False )
    scores.round(3).to_csv(SCORE_DIR / "scib_scores_ranked.csv")
    print(scores[["Bio conservation", "Batch correction", "Total"]].round(3))
    return scores

def score_scatter(scores: pd.DataFrame) -> None:
    """Bio conservation against batch correction, one point per embedding."""
    fig, ax = plt.subplots(figsize = (9, 7))
    ax.scatter(scores["Batch correction"], scores["Bio conservation"], s = 70)

    for name, row in scores.iterrows():
        ax.annotate(
            str(name).replace("X_", ""),
                              (row["Batch correction"], row["Bio conservation"]),
                              xytext = (6, 6),
                              textcoords = "offset points",
                              fontsize = 9,
                              )
    ax.set_xlabel("Batch correction")
    ax.set_ylabel("Bio conservation")
    ax.set_title("Bio conservation vs Batch correction")
    ax.grid(alpha = 0.4)
    save_figure(FIG_DIR / "score_scatter.png", fig)

def main() -> None:
    """Benchmark comparison figures and the ranked score table written out."""
    adata = sc.read_h5ad(IN)
    names = list_embeddings()


    umap_grid(adata, names, LABEL_KEY, "umap_grid_cell_type.png")
    umap_grid(adata, names, BATCH_KEY, "umap_grid_donor.png")

    score_scatter(ranked_scores())


if __name__ == "__main__":
    main()










