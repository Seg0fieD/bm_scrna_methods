"""
    RNA velocity on the progenitor-side cell types only: refitted genes and
    graph, dynamical model with differential kinetics, HSPC-rooted orderings.
"""

from __future__ import annotations
import warnings
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scvelo as scv
from anndata import AnnData

from trajectory_utils import (CELL_TYPE_KEY, PROJECT_ROOT, SEED, find_result,
                              save_table, set_seed)
from velocity_scvelo import (LANE_KEY, LINEAGES,
                             MIN_LIKELIHOOD, MIN_SHARED_COUNTS, MODELS,
                             N_JOBS, N_NEIGHBORS, N_TOP_GENES,
                             PALANTIR_COLUMN, VELOCITY_H5AD, VKEYS,
                             cell_table, report_outputs)


FIT_H5AD        = (PROJECT_ROOT / "data" / "interim"
                   / "bm_velocity_progenitor_fit.h5ad")
SCRIPT_NAME     = Path(__file__).stem
WHOLE_SAMPLE    = "velocity_scvelo_paga_edges.csv"

ROOT_KEY        = "hspc_root"
N_PCS           = 30

PROGENITOR_TYPES = (
    "HSPC", "Granulocyte precursor", "Monocyte classical",
    "Monocyte non-classical", "Dendritic cell", "Plasmacytoid DC",
    "Erythroid early", "Erythroid late", "Platelet", "Pre-B large",
    "Pre-B small", "B naive",
)
PROGENITOR_LINEAGES = {name: labels for name, labels in LINEAGES.items()
                       if set(labels) <= set(PROGENITOR_TYPES)}
KNOWN_TRANSITIONS = (
    ("HSPC", "Pre-B large"),
    ("Pre-B large", "Pre-B small"),
    ("Pre-B small", "B naive"),
    ("HSPC", "Erythroid early"),
    ("Erythroid early", "Erythroid late"),
    ("HSPC", "Granulocyte precursor"),
    ("Granulocyte precursor", "Monocyte classical"),
    ("Monocyte classical", "Monocyte non-classical"),
)


def root_cell() -> str:
    """ 
        Barcode of the Palantir pseudotime origin, an HSPC cell, used as the
        velocity root.
    """
    palantir = pd.read_csv(find_result("pt_palantir_pseudotime.csv"),
                           index_col=0)[PALANTIR_COLUMN]
    return str(palantir.idxmin())


def prepared_subset() -> AnnData:
    """Progenitor-side cells with genes, principal components, neighbour
    graph and moments recomputed on the subset alone."""
    adata = ad.read_h5ad(VELOCITY_H5AD)
    adata = adata[adata.obs[CELL_TYPE_KEY].isin(PROGENITOR_TYPES)].copy()
    adata.obs[CELL_TYPE_KEY] = (adata.obs[CELL_TYPE_KEY]
                                .cat.remove_unused_categories())
    del adata.layers["ambiguous"]
    print(f"progenitor subset: {adata.n_obs:,} cells")
    scv.pp.filter_and_normalize(adata, min_shared_counts=MIN_SHARED_COUNTS,
                                layers_normalize=["X", "spliced", "unspliced"])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES,
                                batch_key=LANE_KEY, subset=True)
    sc.pp.pca(adata, n_comps=N_PCS, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca",
                    random_state=SEED)
    scv.pp.moments(adata, n_neighbors=None, n_pcs=None)
    root = root_cell()
    adata.obs[ROOT_KEY] = (adata.obs_names == root).astype(float)
    adata.uns[ROOT_KEY] = int(adata.obs_names.get_loc(root))
    print(f"root cell: {root}")
    return adata


def fitted_subset() -> AnnData:
    """Progenitor velocity object with both models fitted; the saved fit is
    reused when present."""
    if FIT_H5AD.exists():
        print(f"reusing {FIT_H5AD}; delete it to refit")
        return ad.read_h5ad(FIT_H5AD)
    adata = prepared_subset()
    stochastic, dynamical = VKEYS["stochastic"], VKEYS["dynamical"]

    print("fitting stochastic model")
    scv.tl.velocity(adata, mode="stochastic", vkey=stochastic)
    scv.tl.velocity_graph(adata, vkey=stochastic, n_jobs=N_JOBS)

    print("fitting dynamical model on all variable genes")
    scv.tl.recover_dynamics(adata, var_names="all", n_jobs=N_JOBS)
    scv.tl.differential_kinetic_test(adata, var_names="all",
                                     groupby=CELL_TYPE_KEY)
    scv.tl.velocity(adata, mode="dynamical", vkey=dynamical,
                    diff_kinetics=True)
    scv.tl.velocity_graph(adata, vkey=dynamical, n_jobs=N_JOBS)

    for vkey in (stochastic, dynamical):
        scv.tl.velocity_pseudotime(adata, vkey=vkey, root_key=ROOT_KEY)
        scv.tl.velocity_confidence(adata, vkey=vkey)
    scv.tl.latent_time(adata, vkey=dynamical, root_key=ROOT_KEY,
                       min_likelihood=MIN_LIKELIHOOD)
    adata.write_h5ad(FIT_H5AD, compression="gzip")
    print(f"wrote {FIT_H5AD}")
    return adata


def edge_comparison(progenitor: pd.DataFrame) -> pd.DataFrame:
    """Transition confidence per cell-type pair and model in the
    whole-sample and progenitor runs, restricted to progenitor types."""
    whole = pd.read_csv(find_result(WHOLE_SAMPLE))
    whole = whole[whole["source"].isin(PROGENITOR_TYPES)
                  & whole["target"].isin(PROGENITOR_TYPES)]
    keys = ["model", "source", "target"]
    return (whole.merge(progenitor, on=keys, how="outer",
                        suffixes=("_whole_sample", "_progenitor"))
            .fillna(0.0)
            .sort_values(keys))


def direction_check(edges: pd.DataFrame, run: str) -> pd.DataFrame:
    """Forward and reverse transition confidence for each known
    differentiation step, with the direction the run supports."""
    lookup = edges.set_index(["model", "source", "target"])["confidence"]
    rows = []
    for model in MODELS:
        for source, target in KNOWN_TRANSITIONS:
            forward = lookup.get((model, source, target), 0.0)
            reverse = lookup.get((model, target, source), 0.0)
            verdict = ("correct" if forward > reverse
                       else "reversed" if reverse > forward else "absent")
            rows.append({"run": run, "model": model, "source": source,
                         "target": target, "forward": forward,
                         "reverse": reverse, "direction": verdict})
    return pd.DataFrame(rows)


def main() -> None:
    """Progenitor velocity fit, the shared outputs and the comparison with
    the whole-sample run."""
    set_seed()
    warnings.filterwarnings("ignore", category=FutureWarning)
    scv.settings.verbosity = 2
    adata = fitted_subset()
    edges = report_outputs(adata, cell_table(adata), SCRIPT_NAME,
                           PROGENITOR_LINEAGES)

    whole = pd.read_csv(find_result(WHOLE_SAMPLE))
    directions = pd.concat([direction_check(whole, "whole sample"),
                            direction_check(edges, "progenitor")])
    kinetics = adata.var.get("fit_diff_kinetics")
    tables = {
        "edge_comparison": edge_comparison(edges),
        "direction_check": directions,
    }
    for name, table in tables.items():
        save_table(table, SCRIPT_NAME, f"{SCRIPT_NAME}_{name}", index=False)

    print(directions.round(3).to_string(index=False))
    if kinetics is not None:
        flagged = kinetics.fillna("").astype(str).str.len() > 0
        print(f"\ngenes with differential kinetics: {flagged.sum()} of "
              f"{adata.n_vars}")
    print(f"latent time and pseudotime rooted at "
          f"{adata.obs_names[adata.uns[ROOT_KEY]]}")
    print(f"cells per type:\n"
          f"{adata.obs[CELL_TYPE_KEY].value_counts().to_string()}")


if __name__ == "__main__":
    main()