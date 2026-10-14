"""Held-out test mixtures with known composition: random, partial and
donor-composition sets from the test donors in one bulk matrix."""

from __future__ import annotations

import numpy as np
import pandas as pd

import deconvolution_utils as du

SCRIPT           = "pseudobulk_test"
N_PER_SET        = 500
CELLS_PER_SAMPLE = 500
CONCENTRATION    = 100.0
SET_COLORS       = {"random": "#5aa9e6", "partial": "#f4a259",
                    "donor": "#7bc96f"}


def donor_fractions(composition: np.ndarray, n_samples: int,
                    rng: np.random.Generator) -> np.ndarray:
    """Dirichlet draws centred on a donor's observed composition at fixed
    concentration; cell types absent from the donor remain at zero."""
    fractions = np.zeros((n_samples, composition.size))
    present = composition > 0
    fractions[:, present] = rng.dirichlet(
        CONCENTRATION * composition[present], size=n_samples,
    )
    return fractions


def main() -> None:
    """Simulation of the three test sets, their true fractions and summary."""
    du.set_seed()
    rng = np.random.default_rng(du.SEED)
    test = du.load_counts(du.TEST_DONORS)
    labels = test.obs[du.CELL_TYPE_KEY]
    cell_types = list(labels.cat.categories)
    n_types = len(cell_types)

    parts = [
        ("random", "pooled", *du.pseudobulk(
            test, cell_types, du.random_fractions(N_PER_SET, n_types, rng),
            CELLS_PER_SAMPLE, rng)),
        ("partial", "pooled", *du.pseudobulk(
            test, cell_types,
            du.random_fractions(N_PER_SET, n_types, rng, partial=True),
            CELLS_PER_SAMPLE, rng)),
    ]
    for donor in du.TEST_DONORS:
        cells = test[test.obs[du.DONOR_KEY] == donor]
        composition = (cells.obs[du.CELL_TYPE_KEY]
                       .value_counts(normalize=True).reindex(cell_types)
                       .to_numpy())
        targets = donor_fractions(composition,
                                  N_PER_SET // len(du.TEST_DONORS), rng)
        parts.append(("donor", donor, *du.pseudobulk(
            cells, cell_types, targets, CELLS_PER_SAMPLE, rng)))

    expressions, fractions, meta = [], [], []
    for set_name, donor, expression, truth in parts:
        index = [f"{set_name}_{donor}_{i}" for i in range(len(truth))]
        expressions.append(expression.set_axis(index))
        fractions.append(truth.set_axis(index))
        meta.append(pd.DataFrame({"set": set_name, "donor": donor},
                                 index=index))
    expression = pd.concat(expressions)
    fractions = pd.concat(fractions)
    meta = pd.concat(meta)
    meta["set"] = pd.Categorical(meta["set"], categories=list(SET_COLORS))

    path = du.script_dir(du.INTERIM_DIR, SCRIPT) / "test_bulk.txt"
    expression.T.to_csv(path, sep="\t", float_format="%.0f")
    print(f"wrote {path.relative_to(du.PROJECT_ROOT)} "
          f"({expression.shape[1]:,} genes x {expression.shape[0]:,} samples)")

    source_cells = (pd.crosstab(labels, test.obs[du.DONOR_KEY])
                    .reindex(cell_types).add_prefix("cells_"))
    du.save_table(pd.concat([meta, fractions], axis=1), SCRIPT,
                  "pseudobulk_test_fractions")
    du.save_table(pd.concat([source_cells,
                             du.fraction_summary(fractions, meta["set"])],
                            axis=1),
                  SCRIPT, "pseudobulk_test_cell_types")
    title = (f"Test mixtures from {' and '.join(du.TEST_DONORS)}: "
             f"{len(fractions):,} samples of {CELLS_PER_SAMPLE} cells")
    du.save_figure(du.plot_fraction_boxes(fractions, meta["set"], SET_COLORS,
                                          title),
                   SCRIPT, "pseudobulk_test_fractions")


if __name__ == "__main__":
    main()