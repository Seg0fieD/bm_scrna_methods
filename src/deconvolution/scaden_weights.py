""" 
    Dense-layer weights of the trained Scaden sub-models, exported as NumPy
    archives for use outside TensorFlow.
"""

from __future__ import annotations

import numpy as np
import tensorflow as tf

import deconvolution_utils as du

SCRIPT     = "scaden_weights"
FIT_SCRIPT = "scaden_fit"


def main() -> None:
    """Kernels and biases of each sub-model in layer order, one archive per
    sub-model."""
    du.require_legacy_keras()
    model_dir = du.INTERIM_DIR / FIT_SCRIPT / "model"
    out_dir = du.script_dir(du.INTERIM_DIR, SCRIPT)
    for name in du.SCADEN_ARCHITECTURES:
        model = tf.keras.models.load_model(model_dir / name, compile=False)
        weights = model.get_weights()
        np.savez(out_dir / f"{name}.npz", *weights)
        shapes = ", ".join(str(array.shape) for array in weights[::2])
        print(f"{name}: kernels {shapes}")


if __name__ == "__main__":
    main()