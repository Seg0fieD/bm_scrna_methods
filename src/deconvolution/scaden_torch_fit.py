"""PyTorch reimplementation of the Scaden ensemble: weight-transfer check
against TensorFlow Scaden, multi-seed training and prediction."""

from __future__ import annotations

import time
from collections.abc import Iterator
from itertools import batched, combinations, pairwise
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from torch import nn

import deconvolution_utils as du

SCRIPT         = "scaden_torch_fit"
FIT_SCRIPT     = "scaden_fit"
WEIGHTS_SCRIPT = "scaden_weights"
TEST_SCRIPT    = "pseudobulk_test"
SEEDS          = (7, 8, 9, 10, 11)
N_STEPS        = 5000
BATCH_SIZE     = 128
LEARNING_RATE  = 1e-4
ADAM_EPSILON   = 1e-7
SHUFFLE_BUFFER = 1000
DEVICE         = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class ScadenNet(nn.Sequential):
    """Scaden sub-model: ReLU dense layers with dropout and a softmax output
    over cell types, Glorot-uniform weights and zero biases."""

    def __init__(self, n_genes: int, n_types: int, hidden: tuple[int, ...],
                 dropout: tuple[float, ...]) -> None:
        layers: list[nn.Module] = []
        for (n_in, n_out), rate in zip(pairwise((n_genes, *hidden)), dropout,
                                       strict=True):
            layers += [nn.Linear(n_in, n_out), nn.ReLU(), nn.Dropout(rate)]
        super().__init__(*layers, nn.Linear(hidden[-1], n_types),
                         nn.Softmax(dim=1))
        for layer in self.dense_layers():
            nn.init.xavier_uniform_(layer.weight)
            nn.init.zeros_(layer.bias)

    def dense_layers(self) -> list[nn.Linear]:
        """Linear layers in forward order."""
        return [module for module in self if isinstance(module, nn.Linear)]

    @torch.no_grad()
    def load_keras(self, path: Path) -> None:
        """Kernels and biases copied from a Keras weight archive, kernels
        transposed to the PyTorch layout."""
        archive = np.load(path)
        arrays = [archive[f"arr_{i}"] for i in range(len(archive.files))]
        for layer, kernel, bias in zip(self.dense_layers(), arrays[::2],
                                       arrays[1::2], strict=True):
            layer.weight.copy_(torch.from_numpy(kernel.T))
            layer.bias.copy_(torch.from_numpy(bias))


def scale_samples(bulk: pd.DataFrame) -> np.ndarray:
    """Scaden scaling of a samples-by-genes matrix: log2(x + 1), then min-max
    to the unit interval within each sample."""
    logged = np.log2(bulk.to_numpy(dtype=np.float64) + 1)
    low = logged.min(axis=1, keepdims=True)
    span = logged.max(axis=1, keepdims=True) - low
    return ((logged - low) / np.where(span > 0, span, 1.0)).astype(np.float32)


def buffered_order(n_rows: int, rng: np.random.Generator) -> Iterator[int]:
    """Endless row order of tf.data shuffle with a fixed buffer, reshuffled
    on every pass over the data."""
    while True:
        buffer = list(range(min(SHUFFLE_BUFFER, n_rows)))
        for incoming in range(len(buffer), n_rows):
            slot = int(rng.integers(len(buffer)))
            yield buffer[slot]
            buffer[slot] = incoming
        yield from rng.permutation(buffer).tolist()


def full_order(n_rows: int, rng: np.random.Generator) -> Iterator[int]:
    """Endless row order with a complete reshuffle on every pass over the
    data."""
    while True:
        yield from rng.permutation(n_rows).tolist()


ORDERS = {"buffer": buffered_order, "full": full_order}


def batch_stream(n_rows: int, seed: int,
                 order: str = "buffer") -> Iterator[tuple[int, ...]]:
    """Endless mini-batches of training-row indices in the given shuffle
    order."""
    return batched(ORDERS[order](n_rows, np.random.default_rng(seed)),
                   BATCH_SIZE)


def build(name: str, n_genes: int, n_types: int,
          seed: int) -> tuple[ScadenNet, torch.optim.Adam]:
    """Seeded sub-model on the compute device and its Adam optimiser."""
    torch.manual_seed(seed)
    hidden, dropout = du.SCADEN_ARCHITECTURES[name]
    model = ScadenNet(n_genes, n_types, hidden, dropout).to(DEVICE)
    return model, torch.optim.Adam(model.parameters(), lr=LEARNING_RATE,
                                   eps=ADAM_EPSILON)


def run_steps(model: ScadenNet, optimiser: torch.optim.Adam,
              batches: Iterator[tuple[int, ...]], x: torch.Tensor,
              y: torch.Tensor, n_steps: int) -> np.ndarray:
    """Training steps on mean squared error; per-step training loss."""
    losses = torch.empty(n_steps, device=DEVICE)
    model.train()
    for step in range(n_steps):
        index = torch.tensor(next(batches), device=DEVICE)
        loss = nn.functional.mse_loss(model(x[index]), y[index])
        optimiser.zero_grad()
        loss.backward()
        optimiser.step()
        losses[step] = loss.detach()
    return losses.cpu().numpy()


@torch.no_grad()
def predict(model: ScadenNet, x: torch.Tensor) -> np.ndarray:
    """Cell-type fractions predicted in evaluation mode."""
    model.eval()
    return model(x).cpu().numpy()

def train_checkpoints(x: torch.Tensor, y: torch.Tensor, bulk_x: torch.Tensor,
                      name: str, seed: int, order: str,
                      checkpoints: tuple[int, ...]
                      ) -> tuple[dict[int, np.ndarray], np.ndarray]:
    """Sub-model trained to the last checkpoint; bulk predictions at every
    checkpoint and per-step training loss."""
    model, optimiser = build(name, x.shape[1], y.shape[1], seed)
    batches = batch_stream(len(x), seed, order)
    predictions, losses, done = {}, [], 0
    for steps in checkpoints:
        losses.append(run_steps(model, optimiser, batches, x, y, steps - done))
        predictions[steps] = predict(model, bulk_x)
        done = steps
    return predictions, np.concatenate(losses)


def pass_means(losses: dict[tuple, np.ndarray], steps_per_pass: float,
               names: list[str]) -> pd.Series:
    """Mean training loss per pass over the training data for every run,
    indexed by the run keys and the pass number."""
    return pd.concat(
        {key: pd.Series(values).groupby(
            (np.arange(len(values)) // steps_per_pass).astype(int)).mean()
         for key, values in losses.items()},
        names=[*names, "pass"],
    ).rename("loss")

def load_inputs() -> tuple[torch.Tensor, torch.Tensor, torch.Tensor,
                           pd.Index, list[str]]:
    """Scaled training mixtures, their fractions and the scaled test bulk on
    the compute device, with test sample names and cell-type order."""
    training = ad.read_h5ad(du.INTERIM_DIR / FIT_SCRIPT / "training.h5ad")
    cell_types, genes = list(training.uns["cell_types"]), training.var_names
    x = torch.from_numpy(np.asarray(training.X, dtype=np.float32)).to(DEVICE)
    y = torch.from_numpy(
        training.obs[cell_types].to_numpy(dtype=np.float32)).to(DEVICE)
    bulk = pd.read_table(du.INTERIM_DIR / TEST_SCRIPT / "test_bulk.txt",
                         index_col=0)
    bulk_x = torch.from_numpy(scale_samples(bulk.loc[genes].T)).to(DEVICE)
    print(f"{len(x):,} training mixtures, {bulk.shape[1]:,} test mixtures, "
          f"{len(genes):,} genes, device {DEVICE}")
    return x, y, bulk_x, bulk.columns, cell_types


def main() -> None:
    """Weight-transfer check, ensemble training per seed and prediction of the
    test mixtures."""
    du.set_seed()
    x, y, bulk_x, samples, cell_types = load_inputs()

    fit_results = du.RESULTS_DIR / FIT_SCRIPT
    reference = pd.read_csv(fit_results / "scaden_fit_sub_models.csv",
                            index_col=[0, 1])
    transfer = {}
    for name, (hidden, dropout) in du.SCADEN_ARCHITECTURES.items():
        model = ScadenNet(x.shape[1], len(cell_types), hidden, dropout)
        model.load_keras(du.INTERIM_DIR / WEIGHTS_SCRIPT / f"{name}.npz")
        expected = reference.loc[name].reindex(index=samples,
                                               columns=cell_types)
        transfer[name] = float(np.abs(predict(model.to(DEVICE), bulk_x)
                                      - expected.to_numpy()).max())
        print(f"{name} with Scaden weights: largest difference "
              f"{transfer[name]:.2e}")

    started = time.perf_counter()
    members, losses = {}, {}
    for seed in SEEDS:
        for name in du.SCADEN_ARCHITECTURES:
            model, optimiser = build(name, x.shape[1], y.shape[1], seed)
            losses[seed, name] = run_steps(model, optimiser,
                                           batch_stream(len(x), seed), x, y,
                                           N_STEPS)
            members[seed, name] = pd.DataFrame(
                predict(model, bulk_x), index=samples, columns=cell_types)
            print(f"seed {seed} {name}: final loss "
                  f"{losses[seed, name][-1]:.4f}")
    minutes = (time.perf_counter() - started) / 60

    sub_models = pd.concat(members, names=["seed", "model", "sample"])
    ensembles = sub_models.groupby(level=["seed", "sample"], sort=False).mean()
    loss_table = pd.concat({key: pd.Series(values) for key, values
                            in losses.items()},
                           names=["seed", "model", "step"]).rename("loss")
    scaden = pd.read_csv(fit_results / "scaden_fit_predictions.csv",
                         index_col=0).reindex(index=samples,
                                              columns=cell_types)
    seed_gaps = [du.rmse(ensembles.loc[a].to_numpy(),
                         ensembles.loc[b].to_numpy())
                 for a, b in combinations(SEEDS, 2)]

    summary = pd.Series({
        "device": str(DEVICE),
        "training mixtures": len(x),
        "test mixtures": len(samples),
        "genes": x.shape[1],
        "seeds": len(SEEDS),
        "training steps per sub-model": N_STEPS,
        "training minutes": round(minutes, 1),
        **{f"largest difference from Scaden with Scaden weights, {name}":
           value for name, value in transfer.items()},
        "largest deviation of a row sum from 1":
            float((ensembles.sum(axis=1) - 1).abs().max()),
        **{f"RMSE to the Scaden ensemble, seed {seed}":
           du.rmse(ensembles.loc[seed].to_numpy(), scaden.to_numpy())
           for seed in SEEDS},
        "mean RMSE between seed ensembles": float(np.mean(seed_gaps)),
    }, name="value", dtype=object)
    du.save_table(ensembles.loc[du.SEED], SCRIPT,
                  "scaden_torch_fit_predictions")
    du.save_table(ensembles, SCRIPT, "scaden_torch_fit_seeds")
    du.save_table(sub_models, SCRIPT, "scaden_torch_fit_sub_models")
    du.save_table(loss_table.to_frame(), SCRIPT, "scaden_torch_fit_loss")
    du.save_table(summary.to_frame(), SCRIPT, "scaden_torch_fit_summary")


if __name__ == "__main__":
    main()