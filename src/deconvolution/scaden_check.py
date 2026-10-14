"""End-to-end Scaden check on a small pseudobulk set: model persistence,
prediction validity and accuracy against a mean-fraction baseline."""

from __future__ import annotations

import shutil

import numpy as np
import pandas as pd
import tensorflow as tf
from scaden.model.functions import preprocess_h5ad_data
from scaden.predict import prediction
from scaden.train import training

import deconvolution_utils as du

SCRIPT        = "scaden_check"
N_TRAIN       = 400
N_TEST        = 20
TRAIN_CELLS   = 500
TEST_CELLS    = 500
N_STEPS       = 300
BATCH_SIZE    = 128
LEARNING_RATE = 1e-4
VAR_CUTOFF    = 0.1
SUB_MODELS    = ("m256", "m512", "m1024")


def main() -> None:
    """Small-scale simulation, training, prediction and validity checks."""
    du.require_legacy_keras()
    du.set_seed()
    rng = np.random.default_rng(du.SEED)
    work = du.script_dir(du.INTERIM_DIR, SCRIPT)
    raw, processed = work / "training_raw.h5ad", work / "training.h5ad"
    bulk, model = work / "test_bulk.txt", work / "model"
    predictions = work / "predictions.txt"

    train = du.load_counts(du.TRAIN_DONORS)
    test = du.load_counts(du.TEST_DONORS)
    cell_types = list(train.obs[du.CELL_TYPE_KEY].cat.categories)
    n_types = len(cell_types)

    train_fractions = np.vstack([
        du.random_fractions(N_TRAIN // 2, n_types, rng),
        du.random_fractions(N_TRAIN // 2, n_types, rng, partial=True),
    ])
    train_x, train_y = du.pseudobulk(train, cell_types, train_fractions,
                                     TRAIN_CELLS, rng)
    test_x, test_y = du.pseudobulk(test, cell_types,
                                   du.random_fractions(N_TEST, n_types, rng),
                                   TEST_CELLS, rng)

    du.write_training(train_x, train_y, SCRIPT, raw)
    test_x.T.to_csv(bulk, sep="\t")
    genes = du.signature_genes(test_x.T, VAR_CUTOFF)
    preprocess_h5ad_data(raw_input_path=str(raw),
                         processed_path=str(processed), sig_genes=genes)
    shutil.rmtree(model, ignore_errors=True)
    training(data_path=str(processed), train_datasets="",
             model_dir=str(model), batch_size=BATCH_SIZE,
             learning_rate=LEARNING_RATE, num_steps=N_STEPS, seed=du.SEED)
    prediction(model_dir=str(model), data_path=str(bulk),
               out_name=str(predictions), seed=du.SEED)

    reloaded = [tf.keras.models.load_model(model / name, compile=False)
                for name in SUB_MODELS]
    predicted = pd.read_table(predictions, index_col=0)[cell_types]
    truth = test_y.to_numpy()
    error = du.rmse(predicted.to_numpy(), truth)
    error_baseline = du.rmse(train_y.mean().to_numpy(), truth)
    pearson = float(np.corrcoef(predicted.to_numpy().ravel(),
                                truth.ravel())[0, 1])
    row_error = float((predicted.sum(axis=1) - 1).abs().max())

    summary = pd.DataFrame.from_records([
        ("signature genes", len(genes), len(genes) > 0),
        ("sub-models reloaded", len(reloaded), len(reloaded) == 3),
        ("largest deviation of a row sum from 1", row_error,
         row_error < 1e-3),
        ("RMSE, mean-fraction baseline", error_baseline, True),
        ("RMSE, Scaden", error, error < error_baseline),
        ("Pearson r, all fractions", pearson, pearson > 0),
    ], columns=["check", "value", "passed"]).set_index("check")
    du.save_table(summary, SCRIPT, "scaden_check_summary")
    print(summary.to_string(float_format="{:.4g}".format))
    failed = summary.index[~summary["passed"]].tolist()
    print(f"failed checks: {failed}" if failed else "all checks passed")


if __name__ == "__main__":
    main()