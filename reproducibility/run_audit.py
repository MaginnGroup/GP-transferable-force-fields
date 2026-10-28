"""Run retained-data checks, optionally including GP regeneration/refits.

Requires the lightweight environment. --gp-python selects a separate compatible
GP interpreter. All outputs stay in --output-dir; no simulation engines run.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gp-python", type=Path)
    parser.add_argument("--refit", action="store_true", help="Refit all 30 grouped folds with three starts")
    parser.add_argument("--revision", help="Supply the tested git revision when auditing a git archive/export")
    args = parser.parse_args()
    if args.refit and not args.gp_python:
        parser.error("--refit requires --gp-python")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "audit.json").exists():
        parser.error("Choose a fresh audit output directory")
    revision = args.revision
    if not revision:
        result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
        revision = result.stdout.strip() if result.returncode == 0 else "unrecorded-export"
    # Preserve all original simulation/result/model inputs, not just those whose
    # filenames occur in one verification script.
    immutable = sorted(p for directory in ["Build_GPs/analysis", "Opt_ES/analysis"]
                       for p in (ROOT / directory).rglob("*") if p.suffix in (".csv", ".pkl"))
    before = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in immutable}
    jobs = [
        ("inputs", sys.executable, ["reproducibility/check_inputs.py"]),
        ("validation", sys.executable, ["reproducibility/validation_counts.py", "--output-dir", str(output / "validation")]),
        ("dmf", sys.executable, ["scripts/plot_dmf_corrected.py", "--output-dir", str(output / "dmf")]),
        ("saved_results", sys.executable, ["reproducibility/analysis/verify_saved_results.py", "--output-dir", str(output / "saved")]),
    ]
    if args.gp_python:
        gp = str(args.gp_python.resolve())
        jobs += [
            ("model_settings", gp, ["reproducibility/analysis/inspect_saved_models.py", "--output-dir", str(output / "models")]),
            ("inventory", gp, ["reproducibility/analysis/local_estimability.py", "inventory", "--output-dir", str(output / "inventory")]),
            ("local", gp, ["reproducibility/analysis/local_estimability.py", "local", "--output-dir", str(output / "local")]),
            ("global", gp, ["make_svd_sensitivity.py", "--output-dir", str(output / "global")]),
        ]
        if args.refit:
            jobs.append(("grouped_refit", gp, ["reproducibility/analysis/grouped_GP_validation.py", "--restarts", "3", "--maxiter", "2000", "--output-dir", str(output / "grouped")]))
        compare = ["reproducibility/analysis/verify_gp_reruns.py", "--local-results", str(output / "local"), "--global-results", str(output / "global"), "--output-dir", str(output / "comparison")]
        if args.refit:
            compare += ["--grouped-results", str(output / "grouped")]
        jobs.append(("rerun_comparison", sys.executable, compare))
    env = dict(os.environ, MPLCONFIGDIR=str(output / "matplotlib-cache"))
    receipts = []
    def scrub(text):
        return text.replace(str(output), "<OUTPUT>").replace(str(ROOT), "<REPO>").replace(str(Path.home()), "<HOME>")
    for name, interpreter, command in jobs:
        start = time.monotonic()
        result = subprocess.run([interpreter] + command, cwd=ROOT, env=env, capture_output=True, text=True)
        (output / (name + ".txt")).write_text(scrub(result.stdout + result.stderr))
        receipts.append(dict(check=name, command=["gp-python" if interpreter != sys.executable else "python"] + [scrub(x) for x in command],
                             exit_code=result.returncode, seconds=time.monotonic() - start, log=name + ".txt"))
        print(f"{name}: {'PASS' if result.returncode == 0 else 'FAIL'}", flush=True)
    changed = [name for name, digest in before.items() if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest]
    passed = not changed and all(r["exit_code"] == 0 for r in receipts)
    report = dict(tested_revision=revision, timestamp_utc=datetime.now(timezone.utc).isoformat(),
                  python=sys.version, platform=platform.platform(), machine=platform.machine(),
                  packages={n: importlib.metadata.version(n) for n in ["numpy", "pandas", "matplotlib"]},
                  checks=receipts, original_data_models_checked=len(before), changed_original_inputs=changed,
                  original_input_sha256=before, passed=passed,
                  not_run=["GROMACS preprocessing/MD", "Cassandra/GEMC", "Trajectory reanalysis", "Recovery of missing simulations"],
                  qualification="Retained-data, model and analysis audit; not molecular-simulation reproduction.")
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
