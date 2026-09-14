"""Over-representation analysis of cell type marker genes against GO, KEGG and Reactome. """

import gseapy
import pandas as pd 

import pathway_utils as pu

def load_marker_lists():
    """ 
        Return the top marker genes per cell type as a mapping of 
        cell type label to gene list.
    """
    frame = pd.read_csv(pu.RESULTS / "de_celltype_top.csv")
    lists = {}
    for label, group in frame.groupby(pu.LABEL_KEY, observed = True):
        lists[label] = group["gene"].to_list()
    return lists


def load_background():
    """ 
        Return every gene tested in the differential expression step, 
        used as the statistical universe.
    """
    frame = pd.read_csv(pu.RESULTS / "de_celltype_full.csv", usecols = ["gene"])
    return frame["gene"].unique().tolist()

def run_ora(gene_list, gene_sets, background):
    """
        Test one gene list against one gene set collection and return the enrichment table.
    """
    result = gseapy.enrich(
        gene_list = gene_list,
        gene_sets = gene_sets,
        background = background,
        outdir = None,
    )
    return result.results.drop(columns = ["Gene_set"], errors = "ignore")

def main():
    pu.setup_plots()

    markers    = load_marker_lists()
    background = load_background()
    print(f"{len(markers)} cell types, background of {len(background):,} genes")

    collected = []

    for key in pu.GENE_SET_LIBRARIES:
        gene_sets = pu.load_gene_sets(key)
        with pu.step(f"over-representation against {key}"):
            for label, genes in markers.items():
                frame = run_ora(genes, gene_sets, background)
                frame.insert(0, pu.LABEL_KEY, label)
                frame.insert(1, "library", key)
                collected.append(frame)

    results = pd.concat(collected, ignore_index = True)
    pu.save_table(results, "ora_markers_all")

    significant = results[results["Adjusted P-value"] < pu.FDR_CUTOFF]
    ordered     = significant.sort_values([pu.LABEL_KEY, "library", "Adjusted P-value"])
    pu.save_table(ordered, "ora_markers_significant")

    counts = significant.groupby([pu.LABEL_KEY, "library"], observed = True).size().unstack(fill_value = 0)
    print(counts.to_string())

if __name__ == "__main__":
    main()

