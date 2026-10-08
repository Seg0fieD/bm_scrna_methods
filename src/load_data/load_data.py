"""Load the 8 bone marrow files into one AnnData with a donor column."""

from pathlib import Path

import anndata as ad 
import scanpy as sc

import warnings

warnings.filterwarnings("ignore", message = "Variable names are not unique")

ROOT      = Path(__file__).resolve().parents[2]
RAW       = ROOT / "data" / "raw"
OUT       = ROOT / "data" / "interim" / "bm_merged.h5ad"
MIN_GENES = 200

def load_one(path: Path) -> ad.AnnData:
    """One 10x matrix read into AnnData, with low-gene barcodes dropped."""
    a = sc.read_10x_h5(path)
    a.var_names_make_unique()
    sc.pp.filter_cells(a, min_genes = MIN_GENES)
    a.obs_names_make_unique()
    return a


def main() -> None:
    """Eight donor matrices merged into one object written to disk."""
    files = sorted(RAW.glob("MantonBM[0-9]_HiSeq_1_raw_feature_bc_matrix.h5"))
    if len(files) != 8:
        raise SystemExit(f"expected 8 files in {RAW}, found {len(files)}")

    donors = [f.name.split("_")[0] for f in files]

    parts = []

    for donor, path in zip(donors, files):
        a = load_one(path)
        print(donor, a.n_obs, "cells")
        parts.append(a)

    merged = ad.concat(
        parts, label = "donor", keys = donors, index_unique = "-", 
        merge = "first"
    )
    merged.obs["donor"] = merged.obs["donor"].astype("category")

    OUT.parent.mkdir(parents = True, exist_ok = True)
    merged.write_h5ad(OUT, compression = "gzip")
    print(merged)


if __name__ == "__main__":
    main()
