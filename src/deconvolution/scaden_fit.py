"""Scaden ensemble trained on the training mixtures and applied to the test
mixtures: ensemble and per-sub-model cell-type fraction predictions."""

from __future__ import annotations

import shutil
import time

import anndata as ad
import pandas as pd
import tensorflow as tf
from scaden.model.architectures import architectures
from scaden.model.functions import preprocess_h5ad_data
from scaden.model.scaden import Scaden
from scaden.train import training

import deconvolution_utils as du

SCRIPT        = "scaden_fit"
TRAIN_SCRIPT  = "pseudobulk_train"
TEST_SCRIPT   = "pseudobulk_test"
N_STEPS       = 5000
BATCH_SIZE    = 128
LEARNING_RATE = 1e-4
VAR_CUTOFF    = 0.1
SUB_MODELS    = ("m256", "m512", "m1024")


def main() -> None:
    """Gene selection, training, model reload check and prediction."""
    du.require_legacy_keras()
    du.set_seed()
    work = du.script_dir(du.INTERIM_DIR, SCRIPT)
    raw = du.INTERIM_DIR / TRAIN_SCRIPT / "training_raw.h5ad"
    bulk_path = du.INTERIM_DIR / TEST_SCRIPT / "test_bulk.txt"
    processed, model = work / "training.h5ad", work / "model"

    training_genes = ad.read_h5ad(raw, backed="r").var_names
    bulk = pd.read_table(bulk_path, index_col=0)
    shared = bulk.index.intersection(training_genes)
    genes = du.signature_genes(bulk.loc[shared], VAR_CUTOFF)
    print(f"{len(genes):,} of {len(bulk):,} bulk genes selected")
    preprocess_h5ad_data(raw_input_path=str(raw),
                         processed_path=str(processed), sig_genes=genes)

    shutil.rmtree(model, ignore_errors=True)
    started = time.perf_counter()
    training(data_path=str(processed), train_datasets="",
             model_dir=str(model), batch_size=BATCH_SIZE,
             learning_rate=LEARNING_RATE, num_steps=N_STEPS, seed=du.SEED)
    minutes = (time.perf_counter() - started) / 60

    predictions = {}
    for name in SUB_MODELS:
        tf.keras.models.load_model(model / name, compile=False)
        hidden_units, dropout_rates = architectures[name]
        predictions[name] = Scaden(
            model_dir=str(model / name), model_name=name, seed=du.SEED,
            hidden_units=hidden_units, do_rates=dropout_rates,
        ).predict(input_path=str(bulk_path))
    ensemble = sum(predictions.values()) / len(predictions)
    sub_models = pd.concat(predictions, names=["model", "sample"])

    summary = pd.Series({
        "training mixtures": ad.read_h5ad(processed, backed="r").n_obs,
        "test mixtures": len(ensemble),
        "genes": len(genes),
        "training steps per sub-model": N_STEPS,
        "training minutes": round(minutes, 1),
        "largest deviation of a row sum from 1":
            float((ensemble.sum(axis=1) - 1).abs().max()),
    }, name="value", dtype=object)
    du.save_table(ensemble, SCRIPT, "scaden_fit_predictions")
    du.save_table(sub_models, SCRIPT, "scaden_fit_sub_models")
    du.save_table(summary.to_frame(), SCRIPT, "scaden_fit_summary")


if __name__ == "__main__":
    main()