"""Per-cell transcription factor activity inferred from DoRothEA regulons"""

import decoupler as dc 
import pandas as pd

import pathway_utils as pu

CONFIDENCE_LEVELS = ["A", "B", "C"]
TOP_FACTORS_SHOWN  = 5


def load_dorothea():
    """Return the DoRothEA regulon table restricted to the high-confidence levels."""
    return dc.op.dorothea(organism = "human", levels = CONFIDENCE_LEVELS)

def score_cells(adata, net):
    """
        Estimate per-cell transcription factor activity with a univariate linear model
        and return the score matrix
    """
    dc.mt.ulm(data = adata, net = net, verbose = True)
    return adata.obsm["score_ulm"]

def mean_by_cell_type(scores, labels):
    """Return the mean activity of every transcription factor within each cell type"""
    frame = scores.copy()
    frame[pu.LABEL_KEY] = labels.values
    summary = frame.groupby(pu.LABEL_KEY, observed = True).mean()
    return summary.reset_index()

def top_factor(summary):
    """Return the highest-scoring transcription factors per cell type as a readable table"""
    rows = []
    for label, row in summary.set_index(pu.LABEL_KEY).iterrows():
        ordered = row.sort_values(ascending = False).head(TOP_FACTORS_SHOWN)
        rows.append({pu.LABEL_KEY: label, "top_factors": ", ".join(ordered.index)})
    return pd.DataFrame(rows)

def main():
    pu.setup_plots()
    adata = pu.load_annotated()

    net = load_dorothea()
    print(f"DoRothEA: {net['source'].nunique():,} factors, {len(net):,} interactions")

    with pu.step("transcription factor activity per cell"):
        scores = score_cells(adata, net)

    summary = mean_by_cell_type(scores, adata.obs[pu.LABEL_KEY])
    pu.save_cell_scores(scores, "act_dorothea")
    pu.save_table(summary, "act_dorothea_by_celltype")

    print(top_factor(summary).to_string(index = False))

if __name__ == "__main__":
    main()

