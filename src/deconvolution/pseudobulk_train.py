"""Scaden training mixtures from the training donors: random and partial
cell-type compositions with realised fractions as labels."""

from __future__ import annotations

import numpy as np
import pandas as pd

import deconvolution_utils as du

SCRIPT           = "pseudobulk_train"
N_SAMPLES        = 8000
CELLS_PER_SAMPLE = 500
KIND_COLORS      = {"random": "#5aa9e6", "partial": "#f4a259"}


def main() -> None:
    """Simulation of the training mixtures and their composition summary."""
    du.set_seed()
    rng = np.random.default_rng(du.SEED)
    train = du.load_counts(du.TRAIN_DONORS)
    cell_types = list(train.obs[du.CELL_TYPE_KEY].cat.categories)
    half = N_SAMPLES // 2
    targets = np.vstack([
        du.random_fractions(half, len(cell_types), rng),
        du.random_fractions(half, len(cell_types), rng, partial=True),
    ])
    expression, fractions = du.pseudobulk(train, cell_types, targets,
                                          CELLS_PER_SAMPLE, rng)
    kinds = pd.Series(pd.Categorical(np.repeat(list(KIND_COLORS), half),
                                     categories=list(KIND_COLORS)),
                      index=fractions.index, name="kind")

    path = du.script_dir(du.INTERIM_DIR, SCRIPT) / "training_raw.h5ad"
    du.write_training(expression, fractions, SCRIPT, path)
    print(f"wrote {path.relative_to(du.PROJECT_ROOT)} "
          f"({expression.shape[0]:,} x {expression.shape[1]:,})")

    source_cells = (train.obs[du.CELL_TYPE_KEY].value_counts()
                    .reindex(cell_types).rename("source_cells"))
    du.save_table(pd.concat([kinds, fractions], axis=1), SCRIPT,
                  "pseudobulk_train_fractions")
    du.save_table(pd.concat([source_cells,
                             du.fraction_summary(fractions, kinds)], axis=1),
                  SCRIPT, "pseudobulk_train_cell_types")
    title = (f"Training mixtures: {len(fractions):,} samples of "
             f"{CELLS_PER_SAMPLE} cells")
    du.save_figure(du.plot_fraction_boxes(fractions, kinds, KIND_COLORS,
                                          title),
                   SCRIPT, "pseudobulk_train_fractions")


if __name__ == "__main__":
    main()