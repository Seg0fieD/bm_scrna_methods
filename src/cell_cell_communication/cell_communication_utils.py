"""Shared paths, settings and helpers for the ligand-receptor analysis scripts."""

from pathlib import Path
import random

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import liana as li
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ANNOTATED_H5AD     = PROJECT_ROOT / "data" / "processed" / "bm_annotated.h5ad"
GENE_SETS_DIR      = PROJECT_ROOT / "data" / "external" / "gene_sets"
RESOURCE_CACHE_DIR = PROJECT_ROOT / "data" / "external" / "lr_resources"
RESULTS_DIR        = PROJECT_ROOT / "results" / "cell_cell_communication"
FIGURES_DIR        = PROJECT_ROOT / "figures" / "cell_cell_communication"

SEED                = 7
GROUP_KEY           = "cell_type"
DONOR_KEY           = "donor"

RESOURCES           = ("consensus", "cellphonedb")
GENE_SET_LIBRARIES  = ("GO_BP", "KEGG", "Reactome")

EXPR_PROP           = 0.1
MIN_CELLS           = 10
N_PERMS             = 1000
MAGNITUDE_CUTOFF    = 0.05
SPECIFICITY_CUTOFF  = 0.05

DPI                 = 300
FIGSIZE             = (17, 10)
LEGEND_FONTSIZE     = 6
DARK_BACKGROUND     = "#12141a"
DARK_FOREGROUND     = "#e6e6e6"

SCORE_COLUMNS       = {
    "lr_means"          : ("CellPhoneDB", "magnitude", True),
    "cellphone_pvals"   : ("CellPhoneDB", "specificity", False),
    "expr_prod"         : ("NATMI","magnitude", True),
    "spec_weight"       : ("NATMI", "specificity", True),
    "scaled_weight"     : ("Connectome", "specificity",True),
    "lr_logfc"          : ("log2FC", "specificity", True),
    "lrscore"           : ("SingleCellSignalR", "magnitude",True),
    "lr_probs"          : ("CellChat", "magnitude", True),
    "cellchat_pvals"    : ("CellChat", "specificity", False),
    "lr_gmeans"         : ("Geometric mean", "magnitude", True),
    "gmean_pvals"       : ("Geometric mean", "specificity", False),
}


def set_seed(seed: int = SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)

def load_annotated() -> ad.AnnData:
    """Load the annotated dataset"""
    adata = sc.read_h5ad(ANNOTATED_H5AD)
    print(f"loaded {ANNOTATED_H5AD} : {adata.n_obs} cells x {adata.n_vars} genes")
    return adata

def cell_type_labels(adata: ad.AnnData) -> list[str]:
    """Return the cell type labels in category order form."""
    return list(adata.obs[GROUP_KEY].cat.categories)

def load_resource(resource_name: str) -> pd.DataFrame:
    """Return a ligand-receptor resources as dataframe, download once and caches for later"""
    RESOURCE_CACHE_DIR.mkdir(parents = True, exist_ok = True)
    cached = RESOURCE_CACHE_DIR / f"{resource_name}.csv"
    if cached.exists():
        resource = pd.read_csv(cached)
    else:
        resource = li.rs.select_resource(resource_name)
        resource.to_csv(cached, index = False)
        print(f"cached resource -> {cached}")
    print(f"{resource_name}: {len(resource)} interactions")
    return resource

def load_gene_sets(library: str) -> dict[str, list[str]]:
    """Read a cached GMT file and return each gene set name mapped to its gene list."""
    path = GENE_SETS_DIR / f"{library}.gmt"
    gene_sets = {}
    with path.open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2:
                gene_sets[fields[0]] = fields[2:]
    print(f"{library}: {len(gene_sets)} gene sets")
    return gene_sets

def rank_interactions(adata: ad.AnnData, resource_name: str) -> pd.DataFrame:
    """Run the LIANA consensus ranking between cell types on cached resources and return interaction table."""
    resource = load_resource(resource_name)
    li.mt.rank_aggregate(
        adata         = adata,
        groupby       = GROUP_KEY,
        resource      = resource,
        resource_name = resource_name,
        expr_prop     = EXPR_PROP,
        min_cells     = MIN_CELLS,
        n_perms       = N_PERMS,
        seed          = SEED,
        use_raw       = False,
        verbose       = True,
        inplace       = True,
    )
    return adata.uns["liana_res"].copy()


def filter_significant(frame: pd.DataFrame) -> pd.DataFrame:
    """Return the interactions passing both the magnitude and specificity rank cutoffs"""
    keep = (frame["magnitude_rank"] <= MAGNITUDE_CUTOFF) & (
        frame["specificity_rank"] <= SPECIFICITY_CUTOFF)
    return frame.loc[keep].sort_values("magnitude_rank").reset_index(drop = True)

def interaction_counts(frame: pd.DataFrame, labels: list[str] | None = None) -> pd.DataFrame:    
    """    
        Return a sender-by-receiver matrix of interaction counts per cell type pair
        reindexed to `labels` so runs share theri axes.
    """
    counts = pd.crosstab(frame["source"], frame["target"])
    if labels is None:
        labels = sorted(set(counts.index) | set(counts.columns))
    return counts.reindex(index = labels, columns = labels, fill_value = 0)

def present_score_columns(frame: pd.DataFrame) -> dict[str, tuple[str, str, bool]]:
    """
        Return the per-method score columns from ranking table, mapped to method name,
        the score type, better interaction rank based on higher value ranks.
    """
    return {
        column: SCORE_COLUMNS[column]
        for column in SCORE_COLUMNS
        if column in frame.columns
    }

def save_table(frame: pd.DataFrame, script_name: str, name: str, index: bool = False) -> Path:
    """write a results table as CSV and output path"""
    path = RESULTS_DIR / script_name / f"{name}.csv"
    frame.to_csv(path, index=index)
    print(f"wrote {path} ({len(frame)} rows)")
    return path


def style_dark(fig : Figure) -> None:
    """Apply cosmetic dark background and foreground colors to a figure."""
    fig.patch.set_facecolor(DARK_BACKGROUND)
    for ax in fig.get_axes():
        ax.set_facecolor(DARK_BACKGROUND)
        ax.title.set_color(DARK_FOREGROUND)
        ax.xaxis.label.set_color(DARK_FOREGROUND)
        ax.yaxis.label.set_color(DARK_FOREGROUND)
        ax.tick_params(colors = DARK_FOREGROUND)
        for spine in ax.spines.values():
            spine.set_color(DARK_FOREGROUND)
        legend = ax.get_legend()
        if legend is not None:
            legend.get_frame().set_facecolor(DARK_BACKGROUND)
            legend.get_frame().set_edgecolor(DARK_FOREGROUND)
            for text in legend.get_texts():
                text.set_color(DARK_FOREGROUND)


def save_figure(fig: Figure, name: str) -> Path:
    """Write a figure as PNG to figure directory at the fixed resolution."""
    FIGURES_DIR.mkdir(parents=True, exist_ok = True)
    path = FIGURES_DIR / f"{name}.png"
    fig.savefig(path, dpi = DPI, bbox_inches = "tight", facecolor = fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {path}")
    return path

