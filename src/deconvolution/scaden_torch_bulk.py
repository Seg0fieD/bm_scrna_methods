"""Scaden ensemble in PyTorch trained on the genes shared with the sorted
bulk data (GSE74246) and applied to its pure samples and mixtures."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

import deconvolution_utils as du
import scaden_torch_fit as stf

SCRIPT      = "scaden_torch_bulk"
TRAIN_FILE  = du.INTERIM_DIR / "pseudobulk_train" / "training_raw.h5ad"
BULK_DIR    = du.INTERIM_DIR / "corces_bulk"
VARIANTS    = ("counts", "length")
SEEDS       = (7, 8, 9, 10, 11)
ORDER       = "full"
CHECKPOINTS = (5000, 20000)
VAR_CUTOFF  = 0.1


def load_inputs(bulk_file: Path) -> tuple[torch.Tensor, torch.Tensor,
                                          torch.Tensor, pd.Index, list[str]]:
    """Scaled training mixtures, fractions and bulk samples on the signature
    genes, on the compute device, with bulk sample names and cell types."""
    training = ad.read_h5ad(TRAIN_FILE)
    bulk = pd.read_table(bulk_file, index_col=0)
    shared = bulk.index.intersection(training.var_names)
    genes = du.signature_genes(bulk.loc[shared], VAR_CUTOFF)
    cell_types = list(training.uns["cell_types"])
    x = stf.scale_samples(pd.DataFrame(training[:, genes].X))
    y = training.obs[cell_types].to_numpy(dtype=np.float32)
    bulk_x = stf.scale_samples(bulk.loc[genes].T)
    print(f"{len(genes):,} of {len(shared):,} shared genes selected, "
          f"{len(x):,} training mixtures, {bulk.shape[1]:,} bulk samples, "
          f"device {stf.DEVICE}")
    return (torch.from_numpy(x).to(stf.DEVICE),
            torch.from_numpy(y).to(stf.DEVICE),
            torch.from_numpy(bulk_x).to(stf.DEVICE), bulk.columns, cell_types)


def main() -> None:
    """Training per seed and sub-model with predictions for the chosen bulk
    variant at every checkpoint; ensembles and per-pass loss."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant", choices=VARIANTS,
                        help="raw counts or gene-length corrected bulk")
    variant = parser.parse_args().variant
    du.set_seed()
    x, y, bulk_x, samples, cell_types = load_inputs(
        BULK_DIR / f"corces_bulk_{variant}.txt")

    started = time.perf_counter()
    members, losses = {}, {}
    for seed in SEEDS:
        for name in du.SCADEN_ARCHITECTURES:
            predictions, losses[seed, name] = stf.train_checkpoints(
                x, y, bulk_x, name, seed, ORDER, CHECKPOINTS)
            for steps, values in predictions.items():
                members[seed, steps, name] = pd.DataFrame(
                    values, index=samples, columns=cell_types)
        print(f"seed {seed} done, "
              f"{(time.perf_counter() - started) / 60:.1f} min")
    minutes = (time.perf_counter() - started) / 60

    ensembles = pd.concat(
        members, names=["seed", "steps", "model", "sample"],
    ).groupby(level=["seed", "steps", "sample"], sort=False).mean()
    summary = pd.Series({
        "device": str(stf.DEVICE),
        "bulk variant": variant,
        "genes": x.shape[1],
        "training mixtures": len(x),
        "bulk samples": len(samples),
        "seeds": len(SEEDS),
        "shuffle order": ORDER,
        "checkpoints": ", ".join(map(str, CHECKPOINTS)),
        "training minutes": round(minutes, 1),
        "largest deviation of a row sum from 1":
            float((ensembles.sum(axis=1) - 1).abs().max()),
    }, name="value", dtype=object)
    prefix = f"{SCRIPT}_{variant}"
    du.save_table(ensembles, SCRIPT, f"{prefix}_predictions")
    du.save_table(stf.pass_means(losses, len(x) / stf.BATCH_SIZE,
                                 ["seed", "model"]).to_frame(),
                  SCRIPT, f"{prefix}_loss")
    du.save_table(summary.to_frame(), SCRIPT, f"{prefix}_summary")


if __name__ == "__main__":
    main()