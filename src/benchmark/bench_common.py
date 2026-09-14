"""Shared paths, settings, and embedding storages for the dimensionality-reduiction benchmark"""

from pathlib import Path

import matplotlib.pyplot as plt 
import numpy as  np 
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
ANNOTATED = ROOT / "data" / "processed" / "bm_annotated.h5ad"
EMB_DIR = ROOT / "results" / "benchmark" / "embeddings"
TAB_DIR = ROOT / "results" / "benchmark"
FIG_DIR = ROOT / "figures" / "benchmark"

for _d in (EMB_DIR, TAB_DIR, FIG_DIR):
    _d.mkdir(parents = True, exist_ok = True)

N_COMPS = 30
N_NEIGHBORS = 15
RESOLUTION = 1.0
SEED = 7
BATCH_KEY = "donor"
LABEL_KEY = "cell_type"

DPI  = 300
FIGSIZE = (15, 8)
LEGEND_FONTSIZE = 6

def load_annotated_data():
    """Load annotated dataset, fail early if variable genes or raw counts are missing"""

    adata = sc.read_h5ad(ANNOTATED)

    if "highly_variable" not in adata.var:
        raise KeyError("hhighly_variable missing from the bm_annotated.h5ad")
    
    if "counts" not in adata.layers:
        raise KeyError("layers['counts'] missing from bm_annotated.h5ad")

    return adata 

def save_embedding(name, X):
    """Write one embedding to disk as float32"""
    X = np.asarray(X, dtype = np.float32)
    np.save(EMB_DIR / f"{name}.npy", X)
    print(f"saved {name}: {X.shape[0]} cells x {X.shape[1]} dims")

def load_embedding(name):
    """Load a previous saved embeddings by name."""
    return np.load(EMB_DIR / f"{name}.npy")

def list_embeddings():
    """Return the names of all saved embeddings"""
    return sorted(p.stem for p in EMB_DIR.glob("*.npy"))

def save_figure(path, fig = None):
    """Save a figures"""
    path = Path(path)
    path.parent.mkdir(parents = True, exist_ok = True)
    fig = plt.gcf() if fig is None else fig
    fig.savefig(path, dpi = DPI, bbox_inches = "tight")
    plt.close(fig)
    print(f"wrote {path}")

