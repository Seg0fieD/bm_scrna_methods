"""Label clusters from marker genes and merge the ones that match."""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

IN = Path("data/interim/bm_markers.h5ad")
OUT = Path("data/processed/bm_annotated.h5ad")
FIG = Path("figures/annotate")

GROUP = "leiden_1.0"

LABELS = {
    "0": "CD4 T naive",
    "2": "CD4 T naive",
    "1": "CD4 T memory",
    "20": "T activated",
    "6": "CD8 T memory",
    "7": "CD8 T cytotoxic",
    "11": "NK",
    "3": "Monocyte classical",
    "4": "Monocyte non-classical",
    "14": "Dendritic cell",
    "17": "Plasmacytoid DC",
    "5": "Granulocyte precursor",
    "9": "HSPC",
    "8": "Pre-B large",
    "10": "Pre-B small",
    "12": "B naive",
    "13": "B naive",
    "15": "Erythroid late",
    "16": "Erythroid early",
    "18": "Platelet",
    "19": "T proliferating",
}

ORDER = [
    "HSPC", "Granulocyte precursor", "Monocyte classical",
    "Monocyte non-classical", "Dendritic cell", "Plasmacytoid DC",
    "Erythroid early", "Erythroid late", "Platelet",
    "Pre-B large", "Pre-B small", "B naive",
    "CD4 T naive", "CD4 T memory", "T activated",
    "CD8 T memory", "CD8 T cytotoxic", "NK", "T proliferating",
]

CANONICAL = {
    "HSPC": ["SOX4", "PRSS57"],
    "Granulocyte": ["MPO", "AZU1", "ELANE"],
    "Monocyte": ["LYZ", "S100A8", "LST1", "AIF1"],
    "DC": ["FCER1A", "CST3"],
    "pDC": ["JCHAIN", "IRF7", "MZB1"],
    "Erythroid": ["HBB", "AHSP", "KLF1", "CA1"],
    "Platelet": ["PF4", "PPBP", "GP9"],
    "B lineage": ["VPREB1", "IGLL1", "CD79A", "MS4A1", "IGHD"],
    "T": ["CD3D", "CD3E", "TRAC", "IL7R"],
    "CD8 / NK": ["CD8A", "GZMK", "GZMH", "GNLY", "KLRF1"],
    "Cycling": ["MKI67", "PCLAF", "STMN1"],
}


def main():
    a = sc.read_h5ad(IN)

    a.obs["cell_type"] = (
        a.obs[GROUP].map(LABELS).astype("category").cat.reorder_categories(ORDER)
    )

    print("cells per label")
    print(a.obs["cell_type"].value_counts().reindex(ORDER).to_string())
    print("\nunlabelled:", int(a.obs["cell_type"].isna().sum()))

    FIG.mkdir(parents = True, exist_ok=True)

    sc.pl.umap(a, color = "cell_type", legend_loc = "on data",
               legend_fontsize = 6, show = False)
    plt.gcf().savefig(FIG / "umap_cell_type.png", dpi = 350,
                      bbox_inches = "tight")
    plt.close("all")

    sc.pl.dotplot(a, CANONICAL, groupby = "cell_type",
                  standard_scale = "var", show = False)
    plt.gcf().savefig(FIG / "dotplot_canonical.png", dpi = 350,
                      bbox_inches = "tight")
    plt.close("all")

    OUT.parent.mkdir(parents = True, exist_ok = True)
    a.write_h5ad(OUT, compression = "gzip")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()