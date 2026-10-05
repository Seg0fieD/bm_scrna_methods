"""Shared paths, settings and helper functions for pathway modules"""

from contextlib import contextmanager
from pathlib import Path
import re
import time 

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import scanpy as sc 
import gseapy

ROOT = Path(__file__).resolve().parents[2]

ANNOTATED      = ROOT / "data" / "processed" / "bm_annotated.h5ad"
CNMF_TOP_GENES = ROOT / "results" / "benchmark" / "emb_cnmf" / "cnmf_top_genes.csv"
GENE_SET_DIR   = ROOT / "data" / "external" / "gene_sets"
RESULTS        = ROOT / "results" / "pathway"
FIGURES        = ROOT / "figures" / "pathway" 


for _folder in (GENE_SET_DIR, RESULTS, FIGURES):
    _folder.mkdir(parents = True, exist_ok = True)

SEED            = 7
LABEL_KEY       = "cell_type"
EMBEDDING_KEY   = "X_pca"
COUNTS_LAYER    = "counts"
 
FDR_CUTOFF      = 0.05
TOP_N_MARKER    = 200
 
DPI             = 300
FIGSIZE         = (15, 8)
LEGEND_FONTSIZE = 6

GENE_SET_LIBRARIES = {
    "GO_BP"    : "GO_Biological_Process_2023",
    "KEGG"     : "KEGG_2021_Human",
    "Reactome" : "Reactome_Pathways_2024",
}

CONFOUNDING_PATTERNS = (r"^RP[SL]", r"^MT-", r"^HB[ABDEGMQZ][0-9]?$", r"^(FAU|UBA52|RACK1)$")

def find_result(name):
    """Saved result file by name: the pathway results root first, then one
    level of per-script subfolder."""
    direct = RESULTS / name
    if direct.exists():
        return direct
    nested = sorted(RESULTS.glob(f"*/{name}"))
    if nested:
        return nested[0]
    raise FileNotFoundError(f"result not found under {RESULTS}: {name}")

def setup_plots():
    """Apply the project-wide figure size and resolution to matplotlib and scanpy"""

    mpl.rcParams["figure.figsize"] = FIGSIZE
    mpl.rcParams["figure.dpi"]     = 200
    mpl.rcParams["savefig.dpi"]    = DPI
    mpl.rcParams["savefig.bbox"]   = "tight"
    sc.settings.verbosity          = 1
    np.random.seed(SEED)


def load_annotated():
    """Load the annotated dataset, fail early when no labels, raw counts or PCA not found"""
    if not ANNOTATED.exists():
        raise FileNotFoundError(f"annotated dataset not found: {ANNOTATED}")
    
    adata = sc.read_h5ad(ANNOTATED)

    if LABEL_KEY not in adata.obs:
        raise KeyError(f"obs is missing '{LABEL_KEY}'")
    if COUNTS_LAYER not in adata.layers:
        raise KeyError(f"layers is missing '{COUNTS_LAYER}'")
    if EMBEDDING_KEY not in adata.obsm:
        raise KeyError(f"obsm is missing '{EMBEDDING_KEY}'")


    print(f"loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes , "
          f"{adata.obs[LABEL_KEY].nunique()} labels")
    return adata

def drop_confounding_genes(genes):
    """
        Remove ribosomal, mitochondrial and haemoglobin genes from the 
        list of gene names.
        These genes reflect cell state rather than cell identity and 
        otherwise dominate ranked gene lists.
    """
    pattern = re.compile("|".join(CONFOUNDING_PATTERNS))
    kept = [gene for gene in genes if not pattern.match(gene)]
    print(f"dropped {len(genes) - len(kept)} ribosomal/mitochondrial/haemoglobin genes")
    return kept

def write_gmt(gene_sets, path):
    """Save gene sets to a GMT file: set name, a "na" description field, then genes, tab-seperated"""
    with open(path, "w") as f:
        for name, genes in gene_sets.items():
            f.write("\t".join([name, "na"] + list(genes)) + "\n")

def read_gmt(path):
    """Load gene sets from a GMT file as set name to gene list, ignoring the description field"""
    gene_sets = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2:
                gene_sets[fields[0]] = fields[2:]
    return gene_sets

def load_gene_sets(key):
    """
        Return one gene set collection, downloading and caching it on first use.
        `key` is one of the short names in GENE_SET_LIBRARIES.
    """

    library = GENE_SET_LIBRARIES[key]
    cached  = GENE_SET_DIR / f"{library}.gmt"

    if cached.exists():
        gene_sets = read_gmt(cached)
        print(f"{key} : {len(gene_sets)} sets from cache")
        return gene_sets

    gene_sets = gseapy.get_library(name = library, organism = "Human")
    write_gmt(gene_sets, cached)
    print(f"{key} : {len(gene_sets)} sets downloaded and cached to {cached.name}")
    return gene_sets

def check_libraries():
    """Check the configured gene set library names still exist on the Enrichr server"""

    available = set(gseapy.get_library_name(organism = "Human"))
    missing = {key: name for key, name in GENE_SET_LIBRARIES.items() if name not in available}

    if missing: 
        for key, name in missing.items():
            close = sorted(n for n in available if n.split("_20")[0] == name.split("_20")[0])
            print(f"{key}: '{name}' not found . Closest matches: {close}")
        raise ValueError("gene set library names are out of date, update GENE_SET_LIBRARIES")

    print("all gene set library names are valid")

def save_table(frame, script_name, name):
    """Results table written as CSV into the script's own results subfolder."""
    folder = RESULTS / script_name
    folder.mkdir(parents = True, exist_ok = True)
    path = folder / f"{name}.csv"
    frame.to_csv(path, index = False)
    print(f"wrote {path.relative_to(ROOT)} ({len(frame):,} rows)")

def save_figure(figure, name):
    """Save figure at the project resolution and close it"""
    path = FIGURES / f"{name}.png"
    figure.savefig(path, dpi = DPI, bbox_inches = "tight")
    plt.close(figure)
    print(f"wrote {path.relative_to(ROOT)}")

def save_cell_scores(frame, script_name, name):
    """
        Per-cell score matrix written as a NumPy array plus a CSV of its column
        names, in the script's own results subfolder. Row order matches the
        annotated dataset, so scores reattach by position.
    """
    folder = RESULTS / script_name
    folder.mkdir(parents = True, exist_ok = True)
    np.save(folder / f"{name}.npy", frame.to_numpy(dtype = np.float32))
    pd.Series(frame.columns, name = "feature").to_csv(
        folder / f"{name}_features.csv", index = False
    )
    print(f"wrote {name}.npy ({frame.shape[0]:,} cells x {frame.shape[1]} features)")

def load_cell_scores(name, obs_names = None):
    """Per-cell score matrix loaded by name from any script's results subfolder."""
    values = np.load(find_result(f"{name}.npy"))
    features = pd.read_csv(find_result(f"{name}_features.csv"))["feature"].tolist()
    return pd.DataFrame(values, index = obs_names, columns = features)

@contextmanager
def step(label):
    """Log the start and the runtime of a named block of work"""
    print(f"---{label}")
    started = time.perf_counter()
    yield
    print(f"---{label} finished in {time.perf_counter() - started:.1f}s")

