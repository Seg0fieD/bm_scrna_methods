# **********need to write 
"""Over-representation analysis of cNMF gene expression programs against GO, KEGG and Reactome."""

import gseapy
import pandas as pd

import pathway_utils as pu

CNMF_OVERDISPERSED = pu.ROOT / "results" / "benchmark" / "cnmf" / "bm" / "bm.overdispersed_genes.txt"


def load_programs():
    """Return the top genes of each cNMF program as a mapping of program label to gene list."""
    frame = pd.read_csv(pu.CNMF_TOP_GENES, index_col = 0)
    programs = {}
    for column in frame.columns:
        programs[f"program_{int(column):02d}"] = frame[column].dropna().to_list()
    return programs


def load_background():
    """Return the overdispersed genes the programs were fitted on, used as the statistical universe."""
    genes = pd.read_csv(CNMF_OVERDISPERSED, header = None)[0].to_list()
    return genes


def run_ora(gene_list, gene_sets, background):
    """Test one gene list against one gene set collection and return the enrichment table."""
    result = gseapy.enrich(
        gene_list  = gene_list,
        gene_sets  = gene_sets,
        background = background,
        outdir     = None,
    )
    return result.results.drop(columns = ["Gene_set"], errors = "ignore")


def main():
    pu.setup_plots()

    programs   = load_programs()
    background = load_background()
    print(f"{len(programs)} programs, background of {len(background):,} genes")

    collected = []
    for key in pu.GENE_SET_LIBRARIES:
        gene_sets = pu.load_gene_sets(key)
        with pu.step(f"over-representation against {key}"):
            for label, genes in programs.items():
                frame = run_ora(genes, gene_sets, background)
                frame.insert(0, "program", label)
                frame.insert(1, "library", key)
                collected.append(frame)

    results = pd.concat(collected, ignore_index = True)
    pu.save_table(results, "ora_cnmf_all")

    significant = results[results["Adjusted P-value"] < pu.FDR_CUTOFF]
    ordered     = significant.sort_values(["program", "library", "Adjusted P-value"])
    pu.save_table(ordered, "ora_cnmf_significant")

    best = ordered[ordered["library"] == "Reactome"].groupby("program").head(1)
    print(best[["program", "Term", "Overlap", "Adjusted P-value"]].to_string(index = False))


if __name__ == "__main__":
    main()