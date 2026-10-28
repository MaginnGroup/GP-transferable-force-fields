"""Check final saved GP forms and record actual noise/mean settings.

Uses the GP environment and trusted repository pickles; does not refit models.
"""
import argparse
import hashlib
import json
from pathlib import Path

import gpflow
import numpy as np

from local_estimability import REPO, MOLS, PROPS, load_model
from fffit.fffit.models import MaskedLinear


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for molecule in MOLS:
        for prop in PROPS:
            model, source = load_model(molecule, prop)
            kernels = getattr(model.kernel, "kernels", [model.kernel])
            white = [k for k in kernels if isinstance(k, gpflow.kernels.White)]
            if prop == "liq_density":
                valid = (isinstance(model.kernel, gpflow.kernels.RationalQuadratic)
                         and not white and model.likelihood.variance.trainable
                         and isinstance(model.mean_function, gpflow.mean_functions.Linear))
            else:
                expected_mean = gpflow.mean_functions.Linear if molecule == "DMF" else MaskedLinear
                valid = (len(kernels) == 2 and isinstance(kernels[0], gpflow.kernels.Matern32)
                         and len(white) == 1 and not white[0].variance.trainable
                         and not model.likelihood.variance.trainable
                         and np.isclose(float(model.likelihood.variance.numpy()), 1e-5, rtol=0, atol=1e-12)
                         and isinstance(model.mean_function, expected_mean))
                if isinstance(model.mean_function, MaskedLinear):
                    columns = model.data[0].shape[1]
                    count = (columns - 1) // 2
                    np.testing.assert_array_equal(model.mean_function.mask,
                                                  np.r_[np.zeros(count), np.ones(columns - count)])
            if not valid:
                raise ValueError(f"Saved-model form differs from the manuscript: {molecule} {prop}")
            rows.append(dict(molecule=molecule, property=prop,
                source=source.relative_to(REPO).as_posix(), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                components=[type(k).__name__ for k in kernels], mean=type(model.mean_function).__name__,
                white_variances=[float(k.variance.numpy()) for k in white],
                white_trainable=[k.variance.trainable for k in white],
                likelihood_variance=float(model.likelihood.variance.numpy()),
                likelihood_trainable=model.likelihood.variance.trainable))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "saved_model_settings.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("Verified final kernel, mean and noise settings for all 12 saved GPs. No refits or simulations run.")


if __name__ == "__main__":
    main()
