"""Scaden training variants in PyTorch: shuffle order of the training mixtures
and training length, over several seeds."""

from __future__ import annotations

import time

import pandas as pd

import deconvolution_utils as du
import scaden_torch_fit as stf

SCRIPT      = "scaden_torch_variants"
ORDERS      = ("buffer", "full")
SEEDS       = (7, 8, 9, 10, 11)
CHECKPOINTS = (5000, 10000, 20000)


def main() -> None:
    """Training per shuffle order, seed and sub-model with test predictions at
    every checkpoint; ensembles, per-pass loss and a reproduction check."""
    du.set_seed()
    x, y, bulk_x, samples, cell_types = stf.load_inputs()

    started = time.perf_counter()
    members, losses = {}, {}
    for order in ORDERS:
        for seed in SEEDS:
            for name in du.SCADEN_ARCHITECTURES:
                predictions, losses[order, seed, name] = stf.train_checkpoints(
                    x, y, bulk_x, name, seed, order, CHECKPOINTS)
                for steps, values in predictions.items():
                    members[order, seed, steps, name] = pd.DataFrame(
                        values, index=samples, columns=cell_types)
            print(f"{order} seed {seed} done, "
                  f"{(time.perf_counter() - started) / 60:.1f} min")
    minutes = (time.perf_counter() - started) / 60

    ensembles = pd.concat(
        members, names=["order", "seed", "steps", "model", "sample"],
    ).groupby(level=["order", "seed", "steps", "sample"], sort=False).mean()
    pass_loss = stf.pass_means(losses, len(x) / stf.BATCH_SIZE,
                               ["order", "seed", "model"])
    reference = pd.read_csv(
        du.RESULTS_DIR / stf.SCRIPT / "scaden_torch_fit_seeds.csv",
        index_col=[0, 1],
    )
    baseline = ensembles.loc["buffer"].xs(stf.N_STEPS, level="steps")
    summary = pd.Series({
        "device": str(stf.DEVICE),
        "seeds": len(SEEDS),
        "checkpoints": ", ".join(map(str, CHECKPOINTS)),
        "training minutes": round(minutes, 1),
        "largest difference from scaden_torch_fit, buffer order at "
        f"{stf.N_STEPS} steps":
            float((baseline - reference.loc[baseline.index, cell_types])
                  .abs().to_numpy().max()),
        "largest deviation of a row sum from 1":
            float((ensembles.sum(axis=1) - 1).abs().max()),
    }, name="value", dtype=object)
    du.save_table(ensembles, SCRIPT, "scaden_torch_variants_predictions")
    du.save_table(pass_loss.to_frame(), SCRIPT, "scaden_torch_variants_loss")
    du.save_table(summary.to_frame(), SCRIPT, "scaden_torch_variants_summary")


if __name__ == "__main__":
    main()