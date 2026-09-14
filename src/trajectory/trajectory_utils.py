"""Shared paths, settings and helper funtions for bone marrow trajectroy analysis."""

from __future__ import annotations

from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc 
from anndata import AnnData

PROJECT_ROOT   = Path(__file__).resolve().parents[2]
ANNOTATED_H5AD = PROJECT_ROOT / "data" / "processed" / "bm_annotated.h5ad"
EMBEDDING_DIR  = PROJECT_ROOT / "results" / "benchmark" / "embeddings"
RESULT_DIR    = PROJECT_ROOT / "results" / "trajectory"
FIGURE_DIR     = PROJECT_ROOT / "figures" / "trajectory"

SEED                 = 7
N_NEIGHBORS          = 15
DPI                  = 300
FIGSIZE              = (15, 8)
UMAP_LEGEND_FONTSIZE = 6

CELL_TYPE_KEY        = "cell_type" 
DONOR_KEY            =  "donor"

EMBEDDINGS           = ("diffmap_pca", "diffmap_scvi")
MAIN_EMBEDDING       = "diffmap_pca"

ROOT_METHODS         =  ("diffusion")
MAIN_ROOT_METHODS    =  
ROOT_CELL_TYPE       =  
STEM_MARKERS         = 

RESULT_DIR.mkdir(parents = True, exist_ok = True)
FIGURE_DIR.mkdir(parents = True, exist_ok = True)


plt.rcParams["figure.figsize"] = FIGSIZE
plt.rcParams["figure.dpi"]     = 250
plt.rcParams["savefig.dpi"]    = DPI


def set_seed() -> None:
    """Set seed for reproduciblity"""
    random.seed(SEED)
    np.random.seed(SEED)

def embedding_key(name : str) -> str:
    """Return the ``obsm`` """
    return f"X_{name}"


def run_tag(embedding: str)