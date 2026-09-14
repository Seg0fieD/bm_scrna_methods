"""Pre-ranked gene set enrichment analysis per cell type against GO, KEGG, and Reactome."""

import gseapy 
import pandas as pd 

import pathway_utils as pu

PERMUTATIONS  = 1000
MIN_SET_SIZE  = 15
MAX_SET_SIZE  = 500
THREADS       = 4

def load_rankings():
    """Return the complete gene ranking per cell type, ordered by descending test statistic"""
    frame = pd.read_csv(
        pu.RESULTS / "de_celltype_full.csv",
        usecols = [pu.LABEL_KEY, "gene", "scores"]
    )
    ranking = {}
    for label, group in frame.groupby(pu.LABEL_KEY, observed = True):
        ordered = group.sort_values("scores", ascending = False)
        ranking[label] = ordered[["gene", "scores"]].reset_index(drop = True)
    return ranking

def run_gsea(ranking, gene_sets):
    """Run pre-ranked enrichment on one ranking and return the results table"""
    result = gseapy.prerank(
        rnk             = ranking,
        gene_sets       = gene_sets, 
        permutation_num = PERMUTATIONS,
        min_size        = MIN_SET_SIZE,
        max_size        = MAX_SET_SIZE,
        threads         = THREADS,
        seed            = pu.SEED,
        outdir          = None, 
        no_plot         = True,
        verbose         = False,
    )

    return result.res2d.drop(columns = ["Name"], errors = "ignore")

def results_for_library(key, rankings):
    """Return the enrichment table for one gene set collection, reusing a completed run if present"""
    path = pu.RESULTS / f"gsea_ranked_{key}.csv"
    if path.exists():
        print(f"{key} : reusing {path.name}")
        return pd.read_csv(path)

    gene_sets = pu.load_gene_sets(key)
    collected = []
    with pu.step(f"gene set enrichment against {key}"):
        for label, ranking in rankings.items():
            frame = run_gsea(ranking, gene_sets)
            frame.insert(0, pu.LABEL_KEY, label)
            frame.insert(1, "library", key)
            collected.append(frame)
    results = pd.concat(collected, ignore_index = True)
    pu.save_table(results, f"gsea_ranked_{key}")
    return results


def main():
    pu.setup_plots()

    rankings = load_rankings()
    print(f"{len(rankings)} cell types ranked")

    collected = []
    for key in pu.GENE_SET_LIBRARIES:
        collected.append(results_for_library(key, rankings))

    results = pd.concat(collected, ignore_index = True)
    pu.save_table(results, "gsea_ranked_all")

    significant = results[results["FDR q-val"] < pu.FDR_CUTOFF]
    ordered     = significant.sort_values([pu.LABEL_KEY, "library", "FDR q-val"])
    pu.save_table(ordered, "gsea_ranked_significant")

    counts = significant.groupby([pu.LABEL_KEY, "library"], observed = True).size().unstack(fill_value = 0)
    print(counts.to_string())

if __name__ == "__main__":
    main()

