"""One-versus-rest differential expression per cell type, ranked across all genes."""

import scanpy as sc

import pathway_utils as pu

def rank_genes(adata):
    """
        Perform a one-versus-rest Wilcoxon rank-sum test per cell type and 
        return the results as a single long-format table.
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby = pu.LABEL_KEY,
        method  = "wilcoxon",
        use_raw = False,
        pts     = True,
    )

    frame = sc.get.rank_genes_groups_df(adata, group = None)
    return frame.rename(columns = {"group" : pu.LABEL_KEY, "names" : "gene"})


def top_genes(frame):
    """
        Select the highest-scoring significantly up-regulated genes per cell type for 
        over-representation analysis. 
    """
    raised = frame[(frame["pvals_adj"] < pu.FDR_CUTOFF) & (frame["logfoldchanges"] > 0)]
    ordered = raised.sort_values([pu.LABEL_KEY, "scores"], ascending = [True, False])
    return ordered.groupby(pu.LABEL_KEY, observed = True).head(pu.TOP_N_MARKER)


def main():
    pu.setup_plots()
    adata = pu.load_annotated()

    del adata.layers[pu.COUNTS_LAYER]

    kept = pu.drop_confounding_genes(adata.var_names.to_list())
    adata = adata[:, kept].copy()

    with pu.step("wilcoxon test per cell type"):
        full = rank_genes(adata)

    pu.save_table(full, "de_celltype_full")

    top = top_genes(full)
    pu.save_table(top, "de_celltype_top")

    print(top.groupby(pu.LABEL_KEY, observed = True).size().sort_values(ascending = False).to_string())

if __name__ == "__main__":
    main()
