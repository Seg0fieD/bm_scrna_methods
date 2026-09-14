"""Figures and cross-method comparison tables for the pathway analysis."""

# need to he write it ..

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

import pathway_utils as pu

MIN_LEADING_EDGE     = 5
TOP_FACTORS_PER_TYPE = 3
MIN_RAW_RANGE        = 1.5


def scale_across_cell_types(frame):
    """Standardise each feature across cell types."""
    return (frame - frame.mean()) / frame.std()


def drop_flat_features(frame, min_range = MIN_RAW_RANGE):
    """Drop features whose raw range across cell types falls below `min_range`.

    Standardisation would otherwise amplify negligible variation to full scale.
    """
    spread  = frame.max() - frame.min()
    keep    = spread[spread >= min_range].index
    dropped = sorted(set(frame.columns) - set(keep))
    if dropped:
        print(f"dropped {len(dropped)} flat features: {', '.join(dropped)}")
    return frame[keep]


def leading_edge_size(frame):
    """Add a column giving the number of genes in each leading edge."""
    frame = frame.copy()
    frame["leading_edge_size"] = frame["Lead_genes"].fillna("").str.split(";").str.len()
    return frame


def select_top_factors(scaled):
    """Return the union of each cell type's highest-scoring transcription factors."""
    selected = []
    for label in scaled.index:
        top = scaled.loc[label].sort_values(ascending = False).head(TOP_FACTORS_PER_TYPE)
        selected.extend(top.index)
    return list(dict.fromkeys(selected))


def plot_progeny_umap(adata, scores):
    """Draw one UMAP panel per signalling pathway, coloured by per-cell activity.

    Colour limits are clipped to the 1st and 99th percentile so that extreme cells do not
    compress the visible range.
    """
    columns = []
    for pathway in scores.columns:
        column = f"progeny_{pathway}"
        adata.obs[column] = scores[pathway].values
        columns.append(column)

    sc.pl.umap(
        adata,
        color           = columns,
        ncols           = 4,
        cmap            = "RdBu_r",
        vcenter         = 0,
        vmin            = "p1",
        vmax            = "p99",
        show            = False,
        legend_fontsize = pu.LEGEND_FONTSIZE,
    )

    figure = plt.gcf()
    for axes in figure.axes:
        title = axes.get_title()
        if title.startswith("progeny_"):
            axes.set_title(title.replace("progeny_", ""))

    figure.set_size_inches(20, 16)
    pu.save_figure(figure, "umap_progeny")


def plot_heatmap(frame, name, title):
    """Draw a cell type by feature heatmap, sized to the table, and write it to the figures directory."""
    height = max(6, 0.40 * len(frame))
    width  = max(8, 0.45 * frame.shape[1])

    figure, axes = plt.subplots(figsize = (width, height))
    sns.heatmap(frame, cmap = "RdBu_r", center = 0, ax = axes,
                cbar_kws = {"label": "scaled activity"})
    axes.set_title(title)
    axes.set_xlabel("")
    axes.set_ylabel("")
    pu.save_figure(figure, name)


def overlap_table(ora, gsea):
    """Count terms shared by both gene set methods and unique to each, per cell type and library."""
    rows = []
    for (label, library), ora_group in ora.groupby([pu.LABEL_KEY, "library"], observed = True):
        gsea_group = gsea[(gsea[pu.LABEL_KEY] == label) & (gsea["library"] == library)]
        ora_terms  = set(ora_group["Term"])
        gsea_terms = set(gsea_group["Term"])
        rows.append({
            pu.LABEL_KEY : label,
            "library"    : library,
            "both"       : len(ora_terms & gsea_terms),
            "ora_only"   : len(ora_terms - gsea_terms),
            "gsea_only"  : len(gsea_terms - ora_terms),
        })
    return pd.DataFrame(rows)


def plot_overlap(overlap):
    """Draw the counts from `overlap_table` as stacked bars per cell type, summed over libraries."""
    totals = overlap.groupby(pu.LABEL_KEY, observed = True)[["both", "ora_only", "gsea_only"]].sum()

    figure, axes = plt.subplots(figsize = pu.FIGSIZE)
    totals.plot(kind = "barh", stacked = True, ax = axes)
    axes.set_title("Terms found by both gene set methods, and by each alone")
    axes.set_xlabel("significant terms")
    axes.set_ylabel("")
    axes.legend(title = "")
    pu.save_figure(figure, "ora_vs_gsea_overlap")


def main():
    pu.setup_plots()
    adata = pu.load_annotated()

    progeny = pu.load_cell_scores("act_progeny", adata.obs_names)
    with pu.step("pathway activity maps"):
        plot_progeny_umap(adata, progeny)

    progeny_summary = pd.read_csv(pu.RESULTS / "act_progeny_by_celltype.csv", index_col = 0)
    plot_heatmap(
        scale_across_cell_types(drop_flat_features(progeny_summary)),
        "heatmap_progeny",
        "Signalling pathway activity, scaled within each pathway",
    )

    dorothea_summary = pd.read_csv(pu.RESULTS / "act_dorothea_by_celltype.csv", index_col = 0)
    scaled_factors   = scale_across_cell_types(dorothea_summary)
    plot_heatmap(
        scaled_factors[select_top_factors(scaled_factors)],
        "heatmap_dorothea",
        "Transcription factor activity, scaled within each factor",
    )

    ora  = pd.read_csv(pu.RESULTS / "ora_markers_significant.csv")
    gsea = leading_edge_size(pd.read_csv(pu.RESULTS / "gsea_ranked_significant.csv"))
    pu.save_table(gsea, "gsea_ranked_with_leading_edge")

    strong  = gsea[(gsea["NES"] > 0) & (gsea["leading_edge_size"] >= MIN_LEADING_EDGE)]
    overlap = overlap_table(ora, strong)
    pu.save_table(overlap, "ora_vs_gsea_overlap")
    plot_overlap(overlap)

    print(overlap.groupby(pu.LABEL_KEY, observed = True)[["both", "ora_only", "gsea_only"]].sum().to_string())


if __name__ == "__main__":
    main()