"""Independently check saved grouped predictions and local/global SVD arithmetic.

Requires NumPy and pandas only. Does not load models, fit GPs or run simulations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = Path(__file__).resolve().parent / "reference"
MOLECULES = ("MeOH", "EG", "Gly", "DMF", "DMSO", "DEC")
RTOL, ATOL = 1e-10, 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def same(a, b):
    np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)


def read_csv(path, inputs, root):
    inputs[path.relative_to(root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)


def verify(root, reference):
    inputs, grouped_results = {}, []
    specs_path = reference / "model_specs.json"
    specs_document = json.loads(specs_path.read_text())
    specs = specs_document["molecules"]
    require(hashlib.sha256((root / specs_document["source"]).read_bytes()).hexdigest()
            == specs_document["source_sha256"], "Model constant source changed; review exported specifications")
    inputs[specs_path.relative_to(root).as_posix()] = hashlib.sha256(specs_path.read_bytes()).hexdigest()
    predictions = read_csv(reference / "grouped_validation_predictions_3starts.csv", inputs, root)
    expected_summary = read_csv(reference / "grouped_validation_summary.csv", inputs, root).set_index("molecule")
    all_metrics = read_csv(reference / "grouped_validation_metrics_3starts.csv", inputs, root)
    require(len(predictions) == 745 and set(predictions.molecule) == set(MOLECULES),
            "Expected all 745 grouped predictions across six molecules")
    require(len(all_metrics) == 36 and not all_metrics.duplicated(["molecule", "fold"]).any(),
            "Expected six original tests and 30 distinct grouped-fold metrics")
    require(set(expected_summary.index) == set(MOLECULES) and expected_summary.index.is_unique,
            "Expected one summary per molecule")
    folds_dir = reference / "grouped_validation_3starts"
    protocol = json.loads((folds_dir / "protocol.json").read_text())
    require(protocol["seed"] == 42 and protocol["folds"] == 5 and protocol["maxiter"] == 2000,
            "Unexpected saved grouped-validation protocol")
    for molecule in MOLECULES:
        spec = specs[molecule]
        iterdir = max((root / "Build_GPs/analysis" / molecule / "vle_iters").glob("iter-*/gp_models.pkl"),
                      key=lambda p: int(p.parent.name.split("-")[1])).parent
        names = spec["param_names"] + ["temperature"]
        xs, ys = [], []
        for partition in ("train", "test"):
            xs.append(read_csv(iterdir / "train_test_sets" / f"sim_surf_tens_x_{partition}.csv", inputs, root)[names].to_numpy())
            ys.append(read_csv(iterdir / "train_test_sets" / f"sim_surf_tens_y_{partition}.csv", inputs, root)["sim_surf_tens"].to_numpy())
        x, y = np.vstack(xs), np.concatenate(ys)
        keys = np.round(x[:, :-1], 12)
        vectors, group_ids = np.unique(keys, axis=0, return_inverse=True)
        eligible = 120 if molecule == "DEC" else 125
        require(len(x) == eligible and len(vectors) == eligible // 5, f"Unexpected eligible design: {molecule}")
        require(set(np.unique(group_ids, return_counts=True)[1]) == {5}, f"Incomplete ST vectors: {molecule}")
        membership = read_csv(folds_dir / f"{molecule}_group_membership.csv", inputs, root)
        same(membership["row"], np.arange(len(x)))
        same(membership["group_id"], group_ids)
        saved_vectors = read_csv(folds_dir / f"{molecule}_group_vectors_scaled.csv", inputs, root)
        same(saved_vectors[names[:-1]].to_numpy(), vectors)
        p = predictions[predictions.molecule == molecule].sort_values("row")
        require(len(p) == eligible and p.row.is_unique, f"Predictions are not once per eligible row: {molecule}")
        same(p.row, np.arange(eligible))
        same(p.group_id, group_ids)
        same(p.scaled_target, y)
        lower, upper = spec["surf_tens_bounds"]
        # The original constants compute bounds/spans in float32; preserve the
        # actual span rather than subtracting JSON endpoints in float64.
        span = spec["property_spans"]["surf_tens"]
        same(p.physical_target, y * span + lower)
        same(p.physical_mean, p.scaled_mean * span + lower)
        require(np.isfinite(p.scaled_predictive_variance).all() and (p.scaled_predictive_variance >= 0).all(),
                f"Invalid saved predictive variances: {molecule}")
        same(p.physical_predictive_sd, np.sqrt(p.scaled_predictive_variance) * span)
        assigned = np.array_split(np.random.RandomState(42).permutation(len(vectors)), 5)
        for fold, test_groups in enumerate(assigned):
            heldout_rows = np.flatnonzero(np.isin(group_ids, test_groups))
            train_groups = set(group_ids) - set(test_groups)
            require(not train_groups & set(test_groups), f"Shared train/test vectors: {molecule} fold {fold}")
            pr = p[p.fold == fold].sort_values("row")
            same(pr.row, heldout_rows)
            per_fold = read_csv(folds_dir / f"{molecule}_fold{fold}_predictions.csv", inputs, root).sort_values("row")
            same(pr.drop(columns="molecule").to_numpy(), per_fold.drop(columns="molecule").to_numpy())
            receipts = json.loads((folds_dir / f"{molecule}_fold{fold}_restart_receipts.json").read_text())
            require(len(receipts) == 3, f"Expected three starts: {molecule} fold {fold}")
            converged = [r for r in receipts if r["success"]]
            best = min(converged or receipts, key=lambda r: r["loss"])
            metrics = all_metrics[(all_metrics.molecule == molecule) & (all_metrics.fold == fold)].iloc[0]
            require(int(metrics.selected_restart) == best["restart"], f"Wrong selected restart: {molecule} fold {fold}")
            residual = pr.physical_mean.to_numpy() - pr.physical_target.to_numpy()
            same(metrics.MAE, np.mean(np.abs(residual)))
            same(metrics.RMSE, np.sqrt(np.mean(residual**2)))
            same(metrics.coverage95, np.mean(np.abs(residual) <= 1.959963984540054 * pr.physical_predictive_sd))
        residual = p.physical_mean.to_numpy() - p.physical_target.to_numpy()
        result = dict(molecule=molecule, grouped_rows=len(p), grouped_MAE=float(np.mean(np.abs(residual))),
                      grouped_RMSE=float(np.sqrt(np.mean(residual**2))),
                      grouped_coverage95=float(np.mean(np.abs(residual) <= 1.959963984540054 * p.physical_predictive_sd)))
        for name, value in result.items():
            if name != "molecule":
                same(value, expected_summary.loc[molecule, name])
        grouped_results.append(result)

    # Recalculate every saved local rank from its actual Jacobian, including
    # row-count null directions. Temperature is a row locator, not a column.
    ranks = read_csv(reference / "local_rank_summary.csv", inputs, root)
    require(len(ranks) == 48, "Expected 48 local property/combined diagnostic blocks")
    for row in ranks.itertuples():
        stem = f"{row.molecule}_{row.variant.replace('-', '')}_{row.block}"
        if row.block == "combined":
            stem += "_" + row.weighting
        frame = read_csv(reference / (stem + "_weighted_jacobian.csv"), inputs, root)
        jac = frame.drop(columns=["temperature"], errors="ignore").to_numpy()
        require(jac.shape == (row.rows, row.free_parameters), f"Jacobian dimensions differ: {stem}")
        require(np.isfinite(jac).all(), f"Non-finite Jacobian: {stem}")
        singular = np.linalg.svd(jac, compute_uv=False)
        singular = np.pad(singular, (0, row.free_parameters - len(singular)))
        for threshold, column in [(1e-4, "rank_1e-04"), (1e-6, "rank_1e-06"), (1e-8, "rank_1e-08")]:
            actual = int(np.sum(singular / singular[0] > threshold))
            expected = int(ranks.loc[row.Index, column])
            require(actual == expected, f"Local rank differs: {stem} threshold {threshold}")
        same(singular[0], row.s_max)
        same(singular[-1], row.s_min)

    # Check the public sampled-domain decompositions independently of the local
    # analysis: singular values, orthonormal right vectors and q_j rankings.
    global_blocks = 0
    for molecule in MOLECULES:
        for prop, iterations in [("liq_density", "ld_iters"), ("surf_tens", "vle_iters")]:
            for suffix in ("", "_wo_temp"):
                directory = root / "Build_GPs/analysis" / molecule / iterations / "sens_approx" / ("sens_approx" + suffix)
                s = read_csv(directory / f"sing_val_{prop}.csv", inputs, root)["sing_value"].to_numpy()
                vectors = read_csv(directory / f"basis_vec_{prop}.csv", inputs, root)
                v = vectors.to_numpy()
                require(np.isfinite(s).all() and (s >= 0).all() and (np.diff(s) <= 1e-9).all(), "Invalid global singular values")
                same(v @ v.T, np.eye(len(v)))
                q = s @ np.abs(v)
                order = np.argsort(q)[::-1]
                ranking = read_csv(directory / f"param_rank_{prop}.csv", inputs, root)
                require(ranking.param_rank.tolist() == vectors.columns[order].tolist(), f"Global ranking differs: {directory}")
                same(ranking.sensitivity_score, q[order])
                global_blocks += 1
    return dict(grouped_validation=grouped_results, local_blocks_checked=len(ranks),
                global_decompositions_checked=global_blocks, source_sha256=inputs,
                tolerances=dict(rtol=RTOL, atol=ATOL),
                limits="Saved-result arithmetic and complete-vector fold exclusion only. Does not regenerate GP predictions/gradients, refit models or run molecular simulations.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = verify(ROOT, REFERENCE)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "saved_result_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    pd.DataFrame(report["grouped_validation"]).to_csv(args.output_dir / "grouped_summary.csv", index=False)
    print(f"Verified 30 grouped folds, {report['local_blocks_checked']} local blocks and {report['global_decompositions_checked']} global decompositions. No models loaded or simulations run.")


if __name__ == "__main__":
    main()
