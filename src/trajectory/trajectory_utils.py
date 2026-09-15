"""Shared paths, settings and helper functions for bone marrow Trajectory analysis."""

from __future__ import annotations
import random
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.colorbar import Colorbar

import numpy as np
import pandas as pd
import scanpy as sc 
from anndata import AnnData

PROJECT_ROOT   = Path(__file__).resolve().parents[2]
ANNOTATED_H5AD = PROJECT_ROOT / "data" / "processed" / "bm_annotated.h5ad"
EMBEDDING_DIR  = PROJECT_ROOT / "results" / "benchmark" / "embeddings"
RESULT_DIR     = PROJECT_ROOT / "results" / "trajectory"
FIGURE_DIR     = PROJECT_ROOT / "figures" / "trajectory"

SEED                 = 7
N_NEIGHBORS          = 15
DPI                  = 300
FIGSIZE              = (15, 8)
UMAP_LEGEND_FONTSIZE = 6
BACKGROUND           = "#12141a"
FOREGROUND           = "#e6e6e6"

CELL_TYPE_KEY        = "cell_type" 
DONOR_KEY            =  "donor"

EMBEDDINGS           = ("diffmap_pca", "diffmap_scvi")
MAIN_EMBEDDING       = "diffmap_pca"

ROOT_METHODS         =  ("diffusion", "marker")
MAIN_ROOT_METHODS    =  "diffusion"
ROOT_CELL_TYPE       =  "HSPC"
STEM_MARKERS         = ("CD34", "AVP", "CRHBP", "SPINK2")

RESULT_DIR.mkdir(parents = True, exist_ok = True)
FIGURE_DIR.mkdir(parents = True, exist_ok = True)


plt.rcParams["figure.figsize"] = FIGSIZE
plt.rcParams["figure.dpi"]     = 250
plt.rcParams["savefig.dpi"]    = DPI


def set_seed() -> None:
    """Set seed for python and NumPy"""
    random.seed(SEED)
    np.random.seed(SEED)

def embedding_key(name: str) -> str:
    """Return the ``obsm`` key for an embedding name"""
    return f"X_{name}"


def run_tag(embedding: str, root_method: str) -> str:
    """Return the output prefix for one embedding and root-method pair."""
    return f"{embedding}_{root_method}"

def run_grid() -> list[tuple[str, str]]:
    """Return every embedding and root-method pair, primary combination first."""
    pairs = list(product(EMBEDDINGS, ROOT_METHODS))
    pairs.sort(key = _grid_order) 
    return pairs # pyright: ignore[reportReturnType]

def _grid_order(pair: tuple[str, str]) -> tuple[int, int]:
    """
        Ordering key ranking the primary embedding and root method ahead of the rest 
    """
    embedding, root_method = pair
    return (embedding != MAIN_EMBEDDING, root_method != MAIN_ROOT_METHODS)

def load_annotated() -> AnnData:
    """Read the annotated expression object."""
    adata = sc.read_h5ad(ANNOTATED_H5AD)
    print(f"loaded {ANNOTATED_H5AD.name}: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata

def load_embeddings(name: str) -> np.ndarray:
    """Read precomputed embedding matrix by name"""
    path = EMBEDDING_DIR / f"{name}.npy"
    if not path.exists():
        raise FileNotFoundError(f"embedding not found: {path}")
    return np.load(path)

def attach_embeddings(adata: AnnData, names:tuple[str, ...] = EMBEDDINGS) -> AnnData:
    """
        Store the named embeddings in ``obsm``, rejecting any observation-count mismatch
    """
    for name in names:
        matrix = load_embeddings(name)
        if matrix.shape[0] != adata.n_obs:
            raise ValueError(
                f"{name} has {matrix.shape[0]} rows against {adata.n_obs} cells"
            )
        adata.obsm[embedding_key(name)] = matrix
        print(f"attached {name}: {matrix.shape[0]} x {matrix.shape[1]}")
    return adata

def root_by_diffusion_component(adata: AnnData, embedding: str = MAIN_EMBEDDING) -> str:
    """
        Return the progenitor barcode at the extreme of the first diffusion component, 
        the tail taken from the progenitor mean so the results is invariant to its sign.
    """

    values = adata.obsm[embedding_key(embedding)][:, 0]
    positions = np.flatnonzero(adata.obs[CELL_TYPE_KEY].to_numpy() == ROOT_CELL_TYPE)
    if positions.size == 0:
        raise ValueError(f"no cells labelled {ROOT_CELL_TYPE}")
    subset = values[positions]
    extreme = subset.argmax() if subset.mean() >= values.mean() else subset.argmin()
    return str(adata.obs_names[positions[extreme]])

def root_by_stem_marker(adata: AnnData, embedding: str = MAIN_EMBEDDING) -> str:
    """
        Return the progenitor barcode with the highest stem marker score; ``embedding`` 
        is unused and present for interface symmetry.
    """
    genes = [gene for gene in STEM_MARKERS if gene in adata.var_names]
    if not genes:
        raise ValueError("none of the stem markers are present in the dataset")
    if len(genes) < len(STEM_MARKERS):
        missing = sorted(set(STEM_MARKERS) - set(genes))
        print(f"stem markers absent from the dataset: {', '.join(missing)}")
    sc.tl.score_genes(adata, gene_list = genes, score_name = "stem_score", 
                        random_state = SEED, use_raw = False)
    scores = adata.obs["stem_score"].to_numpy()
    positions = np.flatnonzero(adata.obs[CELL_TYPE_KEY].to_numpy() == ROOT_CELL_TYPE)
    if positions.size == 0:
        raise ValueError(f"no cells labelled {ROOT_CELL_TYPE}")
    return str(adata.obs_names[positions[scores[positions].argmax()]])

ROOT_SELECTORS = {
    "diffusion" : root_by_diffusion_component,
    "marker"    : root_by_stem_marker,
}

def select_root(adata: AnnData, method: str, embedding: str = MAIN_EMBEDDING) -> str:
    """Return the root barcode produced by the named selection method."""
    if method not in ROOT_SELECTORS:
        raise ValueError(f"unknown root method: {method}")
    barcode = ROOT_SELECTORS[method](adata, embedding)
    print(f"root cell by {method} on {embedding}: {barcode}")
    return barcode

def save_table(frame: pd.DataFrame, name: str, index: bool = True) -> Path:
    """Write a table as CSV to the trajectroy results directory."""
    path = RESULT_DIR / f"{name}.csv"
    frame.to_csv(path, index = index)
    print(f"wrote {path} ({len(frame)} rows)")
    return path

def save_matrix(matrix: np.ndarray, name:str, features: list[str]) -> Path:
    """Write a per-cell matrix as ``.npy`` with its features names, in observation order."""
    if matrix.shape[1] != len(features):
        raise ValueError(
            f"{matrix.shape[1]} columns against {len(features)} feature names"
        )
    path = RESULT_DIR / f"{name}.npy"
    np.save(path, matrix)
    pd.DataFrame({"feature": features}).to_csv(
        RESULT_DIR / f"{name}_features.csv", index = False
    )
    print(f"wrote {path} {matrix.shape}")
    return path

def save_figure(figure: Figure, name: str) -> Path:
    """Generate figure for trajectory  as PNG and store under figure  dir."""
    path = FIGURE_DIR / f"{name}.png"
    figure.savefig(path, dpi = DPI, bbox_inches = "tight")
    plt.close(figure)
    print(f"wrote {path}")
    return path

def style_dark(figure: Figure, axis: Axes) -> None:
    """Apply a dark canvas with ligt text, ticks and spines"""
    figure.patch.set_facecolor(BACKGROUND)
    axis.set_facecolor(BACKGROUND)
    axis.title.set_color(FOREGROUND)
    axis.tick_params(colors=FOREGROUND)
    for spine in axis.spines.values():
        spine.set_color(FOREGROUND)

def style_bar(bar: Colorbar) -> None:
    """Apply light text and ticks to a color bar."""
    bar.ax.yaxis.label.set_color(FOREGROUND)
    bar.ax.tick_params(colors = FOREGROUND)

