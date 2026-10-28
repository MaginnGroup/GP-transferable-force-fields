"""Compare newly regenerated GP results against the fixed manuscript references.

This comparison uses NumPy/pandas, independently of TensorFlow or the fit code.
Runtime seconds and optimizer prose are not numerical reproduction targets.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REFERENCE = Path(__file__).resolve().parent / "reference"


def compare(reference, regenerated, columns=None, keys=None):
    old, new = pd.read_csv(reference), pd.read_csv(regenerated)
    if keys:
        old = old.sort_values(keys).reset_index(drop=True)
        new = new.sort_values(keys).reset_index(drop=True)
        if not old[keys].equals(new[keys]):
            raise ValueError(f"Row identities changed: {reference.name}")
    if old.shape != new.shape or old.columns.tolist() != new.columns.tolist():
        raise ValueError(f"Result schema changed: {reference.name}")
    columns = columns or old.select_dtypes("number").columns.tolist()
    a, b = old[columns].to_numpy(), new[columns].to_numpy()
    np.testing.assert_allclose(a, b, rtol=1e-8, atol=1e-10, equal_nan=True)
    finite = np.isfinite(a) & np.isfinite(b)
    return dict(file=reference.name, rows=len(old), columns=columns,
                max_absolute_difference=float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else 0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-results", type=Path)
    parser.add_argument("--grouped-results", type=Path)
    parser.add_argument("--global-results", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not (args.local_results or args.grouped_results or args.global_results):
        parser.error("Supply at least one regenerated result directory")
    checks = []
    if args.local_results:
        checks.append(compare(REFERENCE / "local_GP_predictions.csv", args.local_results / "local_GP_predictions.csv",
                              keys=["molecule", "variant", "property", "temperature"]))
        checks.append(compare(REFERENCE / "local_rank_summary.csv", args.local_results / "local_rank_summary.csv",
                              keys=["molecule", "variant", "block", "weighting"]))
        for path in sorted(REFERENCE.glob("*_weighted_jacobian.csv")):
            checks.append(compare(path, args.local_results / path.name))
    if args.grouped_results:
        name = "grouped_validation_predictions_3starts.csv"
        checks.append(compare(REFERENCE / name, args.grouped_results / name, keys=["molecule", "fold", "row"]))
        name = "grouped_validation_metrics_3starts.csv"
        checks.append(compare(REFERENCE / name, args.grouped_results / name,
                              columns=["train_rows", "test_rows", "iterations", "MAE", "RMSE", "coverage95",
                                       "mean_predictive_sd", "training_loss", "selected_restart", "white_noise_variance_scaled"],
                              keys=["molecule", "fold"]))
        for path in sorted((REFERENCE / "grouped_validation_3starts").glob("*_parameters.json")):
            other = args.grouped_results / "grouped_validation_3starts" / path.name
            old, new = json.loads(path.read_text()), json.loads(other.read_text())
            if old.keys() != new.keys():
                raise ValueError(f"Fitted parameter names changed: {path.name}")
            maximum = 0.0
            for key in old:
                a, b = np.asarray(old[key]), np.asarray(new[key])
                np.testing.assert_allclose(a, b, rtol=1e-8, atol=1e-10)
                maximum = max(maximum, float(np.max(np.abs(a - b))))
            checks.append(dict(file=path.name, max_absolute_difference=maximum))
    if args.global_results:
        root = Path(__file__).resolve().parents[2]
        for path in sorted((root / "Build_GPs/analysis").glob("*/*_iters/sens_approx/sens_approx*/sing_val_*.csv")):
            relative = path.relative_to(root)
            checks.append(compare(path, args.global_results / relative))
            ranking = path.with_name(path.name.replace("sing_val_", "param_rank_"))
            # Singular-vector signs are arbitrary; compare the named q_j scores
            # rather than raw signed basis vectors.
            checks.append(compare(ranking, args.global_results / ranking.relative_to(root),
                                  columns=["sensitivity_score"], keys=["param_rank"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = dict(checks=checks, tolerances=dict(rtol=1e-8, atol=1e-10),
                  maximum_absolute_difference=max(c["max_absolute_difference"] for c in checks),
                  limits="Compares regenerated analysis results, not molecular simulations. Runtime seconds and optimizer messages excluded.")
    (args.output_dir / "gp_rerun_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Compared {len(checks)} regenerated result files; maximum absolute difference {report['maximum_absolute_difference']:.6g}.")


if __name__ == "__main__":
    main()
