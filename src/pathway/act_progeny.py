"""Per-cell signalling pathway activity inferred from PROGENy footprint regulon"""

import decoupler as dc

import pathway_utils as pu

TOP_TARGETS = 500

def load_progeny():
    """Return the PROGENy regulon table, keeping the strongest footprint genes per pathway"""
    return dc.op.progeny(organism = "human", top = TOP_TARGETS)

def score_cells(adata, net):
    """Estimate per-cell pathways activity with a multivariate linear model and retun the score matrix"""
    dc.mt.mlm(data = adata, net = net, verbose = True)
    return adata.obsm["score_mlm"]

def mean_by_cell_type(scores, labels):
    """Return the mean activity of every pathway within each cell type"""
    frame = scores.copy()
    frame[pu.LABEL_KEY] = labels.values
    summary = frame.groupby(pu.LABEL_KEY, observed = True).mean()
    return summary.reset_index()

def main():
    pu.setup_plots()
    adata = pu.load_annotated()

    net = load_progeny()
    print(f"PROGENy: {net["source"].nunique()} pathways, {len(net):,} gene weights")

    with pu.step("pathway activity per cell"):
        scores = score_cells(adata, net)

    pu.save_cell_scores(scores, "act_progeny")

    summary = mean_by_cell_type(scores, adata.obs[pu.LABEL_KEY])
    pu.save_table(summary, "act_progeny_by_celltype")

    print(summary.set_index(pu.LABEL_KEY).round(2).to_string())


if __name__ == "__main__":
    main()

