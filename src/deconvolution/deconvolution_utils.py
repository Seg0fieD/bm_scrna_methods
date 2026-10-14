"""Shared paths, settings and helpers for cell-type deconvolution with
Scaden."""

from __future__ import annotations

import os
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from scipy import sparse

PROJECT_ROOT   = Path(__file__).resolve().parents[2]
ANNOTATED_H5AD = PROJECT_ROOT / "data" / "processed" / "bm_annotated.h5ad"
INTERIM_DIR    = PROJECT_ROOT / "data" / "interim" / "deconvolution"
EXTERNAL_DIR   = PROJECT_ROOT / "data" / "external" / "deconvolution"
RESULTS_DIR    = PROJECT_ROOT / "results" / "deconvolution"
FIGURES_DIR    = PROJECT_ROOT / "figures" / "deconvolution"

SEED          = 7
CELL_TYPE_KEY = "cell_type"
DONOR_KEY     = "donor"
COUNTS_LAYER  = "counts"

TRAIN_DONORS = ("MantonBM2", "MantonBM3", "MantonBM4",
                "MantonBM5", "MantonBM6", "MantonBM7")
TEST_DONORS  = ("MantonBM1", "MantonBM8")

BROAD_GROUPS = {
    "HSPC":                   "HSPC",
    "Erythroid early":        "Erythroid",
    "Erythroid late":         "Erythroid",
    "Granulocyte precursor":  "Myeloid",
    "Monocyte classical":     "Myeloid",
    "Monocyte non-classical": "Myeloid",
    "Dendritic cell":         "Dendritic",
    "Plasmacytoid DC":        "Dendritic",
    "Platelet":               "Platelet",
    "Pre-B large":            "B lineage",
    "Pre-B small":            "B lineage",
    "B naive":                "B lineage",
    "CD4 T naive":            "T cells",
    "CD4 T memory":           "T cells",
    "CD8 T memory":           "T cells",
    "CD8 T cytotoxic":        "T cells",
    "T activated":            "T cells",
    "T proliferating":        "T cells",
    "NK":                     "NK",
}
GROUP_ORDER = ("HSPC", "Erythroid", "Myeloid", "Dendritic", "Platelet",
               "B lineage", "T cells", "NK")

SCADEN_ARCHITECTURES = {
    "m256":  ((256, 128, 64, 32),    (0.0, 0.0, 0.0, 0.0)),
    "m512":  ((512, 256, 128, 64),   (0.0, 0.3, 0.2, 0.1)),
    "m1024": ((1024, 512, 256, 128), (0.0, 0.6, 0.3, 0.1)),
}

DPI             = 300
BACKGROUND      = "#12141a"
FOREGROUND      = "#e6e6e6"
TITLE_FONTSIZE  = 16
LABEL_FONTSIZE  = 13
TICK_FONTSIZE   = 11
LEGEND_FONTSIZE = 12
TITLE_PAD       = 14


def set_seed() -> None:
    """Seeding of the Python and NumPy global random number generators."""
    random.seed(SEED)
    np.random.seed(SEED)


def require_legacy_keras() -> None:
    """Guard for the Keras 2 backend, without which Scaden fails to save a
    model and silently predicts from untrained weights."""
    if os.environ.get("TF_USE_LEGACY_KERAS") != "1":
        raise RuntimeError(
            "TF_USE_LEGACY_KERAS=1 must be set in the environment"
        )


def load_counts(donors: tuple[str, ...]) -> ad.AnnData:
    """Raw-count AnnData of the given donors: CSR float32 counts, cell type
    and donor annotation, in annotated row order."""
    source = ad.read_h5ad(ANNOTATED_H5AD)
    if COUNTS_LAYER not in source.layers:
        raise KeyError(f"layer '{COUNTS_LAYER}' not found in {ANNOTATED_H5AD}")
    unknown = set(donors) - set(source.obs[DONOR_KEY])
    if unknown:
        raise ValueError(f"donors absent from the dataset: {sorted(unknown)}")
    keep = source.obs[DONOR_KEY].isin(donors).to_numpy()
    counts = sparse.csr_matrix(source.layers[COUNTS_LAYER][keep],
                               dtype=np.float32)
    print(f"loaded {counts.shape[0]:,} cells x {counts.shape[1]:,} genes "
          f"from {len(donors)} donors")
    return ad.AnnData(
        X=counts,
        obs=source.obs.loc[keep, [CELL_TYPE_KEY, DONOR_KEY]].copy(),
        var=pd.DataFrame(index=source.var_names),
    )


def random_fractions(n_samples: int, n_types: int, rng: np.random.Generator,
                     partial: bool = False) -> np.ndarray:
    """Cell-type fractions as normalised uniform weights; when partial, each
    sample retains a random subset of one to n_types - 1 cell types."""
    weights = rng.random((n_samples, n_types))
    if partial:
        for row in weights:
            row[rng.permutation(n_types)[rng.integers(1, n_types):]] = 0.0
    return weights / weights.sum(axis=1, keepdims=True)


def mix_cells(counts: sparse.csr_matrix, labels: np.ndarray,
              cell_types: list[str], fractions: np.ndarray, n_cells: int,
              rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Pseudobulk profiles from multinomial cell-type counts and sampling with
    replacement within type: summed raw counts and realised fractions."""
    members = [np.flatnonzero(labels == cell_type) for cell_type in cell_types]
    empty = [cell_types[i] for i, idx in enumerate(members) if idx.size == 0]
    if empty and fractions[:, [cell_types.index(c) for c in empty]].any():
        raise ValueError(f"non-zero fraction for absent cell types: {empty}")
    drawn = rng.multinomial(n_cells, fractions)
    rows, cols = [], []
    for sample, sample_counts in enumerate(drawn):
        for idx, n in zip(members, sample_counts, strict=True):
            if n:
                cols.append(rng.choice(idx, size=n))
                rows.append(np.full(n, sample))
    selection = sparse.csr_matrix(
        (np.ones(sum(map(len, cols)), dtype=np.float32),
         (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(fractions), counts.shape[0]),
    )
    return (selection @ counts).toarray(), drawn / n_cells

def pseudobulk(adata: ad.AnnData, cell_types: list[str],
               fractions: np.ndarray, n_cells: int,
               rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pseudobulk expression (samples by genes) and realised cell-type
    fractions for the given target fractions."""
    expression, realised = mix_cells(
        adata.X, adata.obs[CELL_TYPE_KEY].to_numpy(), cell_types,
        fractions, n_cells, rng,
    )
    index = [f"sample_{i}" for i in range(len(fractions))]
    return (pd.DataFrame(expression, index=index, columns=adata.var_names),
            pd.DataFrame(realised, index=index, columns=cell_types))


def write_training(expression: pd.DataFrame, fractions: pd.DataFrame,
                   dataset: str, path: Path) -> None:
    """Scaden training file: summed counts per sample, cell-type fractions and
    dataset name in obs, the ordered cell-type list in uns."""
    data = ad.AnnData(
        X=expression.to_numpy(dtype=np.float32),
        obs=fractions.assign(ds=dataset),
        var=pd.DataFrame(index=expression.columns),
    )
    data.uns["cell_types"] = list(fractions.columns)
    data.write_h5ad(path)


def signature_genes(bulk: pd.DataFrame, var_cutoff: float) -> list[str]:
    """Genes with sample variance above the cutoff in a bulk matrix (genes by
    samples), sorted for a reproducible feature order."""
    return sorted(bulk.index[bulk.var(axis=1) > var_cutoff])


def rmse(predicted: np.ndarray, truth: np.ndarray) -> float:
    """Root-mean-square error over all samples and cell types."""
    return float(np.sqrt(((predicted - truth) ** 2).mean()))


def broad_groups(fractions: pd.DataFrame) -> pd.DataFrame:
    """Cell-type fractions summed into the broad lineage groups, in fixed
    group order."""
    unknown = set(fractions.columns) - set(BROAD_GROUPS)
    if unknown:
        raise KeyError(f"cell types without a broad group: {sorted(unknown)}")
    grouped = fractions.T.groupby(fractions.columns.map(BROAD_GROUPS)).sum().T
    return grouped.reindex(columns=list(GROUP_ORDER), fill_value=0.0)


def script_dir(root: Path, script: str) -> Path:
    """Per-script output directory under the given root, created on
    demand."""
    path = root / script
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_table(frame: pd.DataFrame, script: str, name: str,
               index: bool = True) -> Path:
    """Results table as CSV in the script's results directory; logs the
    output path and row count."""
    path = script_dir(RESULTS_DIR, script) / f"{name}.csv"
    frame.to_csv(path, index=index)
    print(f"wrote {path.relative_to(PROJECT_ROOT)} ({len(frame):,} rows)")
    return path


def save_figure(figure: Figure, script: str, name: str) -> Path:
    """Figure as PNG at the project resolution in the script's figure
    directory; closes the figure and logs the output path."""
    path = script_dir(FIGURES_DIR, script) / f"{name}.png"
    figure.savefig(path, dpi=DPI, bbox_inches="tight",
                   facecolor=figure.get_facecolor())
    plt.close(figure)
    print(f"wrote {path.relative_to(PROJECT_ROOT)}")
    return path


def style_dark(figure: Figure, axes: list[Axes]) -> None:
    """Dark canvas with light titles, labels, ticks and spines on every
    axis."""
    figure.patch.set_facecolor(BACKGROUND)
    for axis in axes:
        axis.set_facecolor(BACKGROUND)
        axis.title.set_color(FOREGROUND)
        axis.xaxis.label.set_color(FOREGROUND)
        axis.yaxis.label.set_color(FOREGROUND)
        axis.tick_params(colors=FOREGROUND, labelsize=TICK_FONTSIZE)
        for spine in axis.spines.values():
            spine.set_color(FOREGROUND)

def fraction_summary(fractions: pd.DataFrame,
                     groups: pd.Series) -> pd.DataFrame:
    """Per cell type and mixture group: share of mixtures containing the type
    and mean realised fraction."""
    present = (fractions > 0).groupby(groups, observed=True).mean()
    mean = fractions.groupby(groups, observed=True).mean()
    return pd.concat([present.T.add_prefix("present_"),
                      mean.T.add_prefix("mean_fraction_")], axis=1)


def plot_fraction_boxes(fractions: pd.DataFrame, groups: pd.Series,
                        colors: dict[str, str], title: str) -> Figure:
    """Horizontal box plots of fractions per cell type, one box per mixture
    group, on a fixed 0 to 1 axis with the legend in a reserved column."""
    n_types, n_groups = fractions.shape[1], len(colors)
    width = 0.8 / n_groups
    figure, axis = plt.subplots(figsize=(15, 3 + 0.21 * n_types * n_groups))
    style_dark(figure, [axis])
    for rank, (group, color) in enumerate(colors.items()):
        axis.boxplot(
            fractions[groups == group].to_numpy(),
            positions=np.arange(n_types) + (rank - (n_groups - 1) / 2) * width,
            widths=width * 0.9, orientation="horizontal", patch_artist=True,
            boxprops={"facecolor": color, "edgecolor": FOREGROUND},
            medianprops={"color": FOREGROUND},
            whiskerprops={"color": FOREGROUND},
            capprops={"color": FOREGROUND},
            flierprops={"marker": ".", "markersize": 2,
                        "markeredgecolor": color},
        )
    axis.set_yticks(np.arange(n_types), fractions.columns)
    axis.set_ylim(n_types - 0.5, -0.5)
    axis.set_xlim(0, 1)
    axis.set_xlabel("Realised cell-type fraction", fontsize=LABEL_FONTSIZE)
    axis.set_title(title, fontsize=TITLE_FONTSIZE, pad=TITLE_PAD)
    legend = axis.legend(
        handles=[Patch(facecolor=color, edgecolor=FOREGROUND, label=group)
                 for group, color in colors.items()],
        title="Mixture kind", fontsize=LEGEND_FONTSIZE,
        title_fontsize=LEGEND_FONTSIZE, loc="upper left",
        bbox_to_anchor=(1.02, 1.0), frameon=False, labelcolor=FOREGROUND,
    )
    legend.get_title().set_color(FOREGROUND)
    return figure