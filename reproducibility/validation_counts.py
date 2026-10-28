"""Audit retained final-validation rows without importing simulation software.

Run from any directory. Outputs go to a separate directory; source CSVs are read
only. A usable result is a finite value for that property, not a claim that every
property from the same replicate was available or that its trajectory survives.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ("MeOH", "EG", "Gly", "DMF", "DMSO", "DEC")
VARIANTS = {"Base": "ms_val_no_opt", "GP-Optimized": "ms_val_opt"}
PROPERTIES = ("liq_density", "surf_tens", "vap_density", "Pvap", "Hvap")


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def key(row):
    return number(row["temperature"]), number(row["restart"])


def audit(root):
    counts, summaries, discrepancies, inputs = [], [], [], {}
    for molecule in MOLECULES:
        for variant, folder in VARIANTS.items():
            directory = root / "Opt_ES/analysis/at_00" / molecule / "ExpVal/opt_res" / folder
            raw_path, processed_path = directory / "raw_data.csv", directory / "ms_data.csv"
            raw, processed = read(raw_path), read(processed_path)
            if not raw or not processed:
                raise ValueError(f"Empty validation data: {molecule} {variant}")
            for rows, prefix in ((raw, ""), (processed, "sim_")):
                required = {"temperature", "restart"} | {prefix + prop for prop in PROPERTIES}
                if not required <= rows[0].keys():
                    raise ValueError(f"Missing validation columns: {molecule} {variant} {required - rows[0].keys()}")
            for path in (raw_path, processed_path):
                inputs[path.relative_to(root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
            for rows, name in ((raw, "raw"), (processed, "processed")):
                keys = [key(row) for row in rows]
                if any(not all(math.isfinite(x) for x in k) for k in keys):
                    raise ValueError(f"Non-finite temperature/replicate identifier: {molecule} {variant} {name}")
                if len(keys) != len(set(keys)):
                    raise ValueError(f"Duplicate temperature/replicate identifier: {molecule} {variant} {name}")
            processed_by_key = {key(row): row for row in processed}
            if set(map(key, raw)) != set(processed_by_key):
                discrepancies.append({"molecule": molecule, "variant": variant, "issue": "raw/processed row keys differ"})
            temperatures = sorted({key(row)[0] for row in raw})
            for temperature in temperatures:
                rows = [row for row in raw if number(row["temperature"]) == temperature]
                for prop in PROPERTIES:
                    if any(prop not in row for row in rows):
                        raise ValueError(f"Missing property column: {raw_path} {prop}")
                    usable = [row for row in rows if math.isfinite(number(row[prop]))]
                    values = [number(row[prop]) for row in usable]
                    base = dict(molecule=molecule, variant=variant, temperature_K=temperature, property=prop)
                    counts.append(dict(**base, usable_replicates=len(values),
                        replicate_ids=",".join(format(number(row["restart"]), "g") for row in usable),
                        source=raw_path.relative_to(root).as_posix()))
                    mean = statistics.fmean(values) if values else math.nan
                    sd = statistics.pstdev(values) if values else math.nan
                    summaries.append(dict(**base, usable_replicates=len(values), mean=mean,
                        population_sd=sd, errorbar_halfwidth=1.96 * sd))
                    for row in rows:
                        other = processed_by_key.get(key(row), {})
                        a, b = number(row[prop]), number(other.get("sim_" + prop))
                        same = (math.isnan(a) and math.isnan(b)) or math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
                        if not same:
                            discrepancies.append(dict(**base, replicate=number(row["restart"]),
                                raw_value=str(a), processed_value=str(b), issue="raw/processed property differs"))
    return counts, summaries, {"source_sha256": inputs, "discrepancies": discrepancies,
        "checked_variants": len(MOLECULES) * len(VARIANTS), "property_state_points": len(counts),
        "summary_definition": "Finite retained values; arithmetic mean, population SD (ddof=0), 1.96*SD bars. Not confidence intervals for means.",
        "comparison": "Each retained raw property compared with corresponding sim_* processed row, rtol=atol=1e-12; summaries are independently recomputed, not stored-summary comparisons.",
        "not_checked": ["Trajectory equilibration and physical validity", "GROMACS/Cassandra execution", "Missing result causes", "Critical-point fitting"]}


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if isinstance(v, float) and not math.isfinite(v) else v) for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    counts, summaries, report = audit(args.repo_root.resolve())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "replicate_counts.csv", counts)
    write_csv(args.output_dir / "replicate_summaries.csv", summaries)
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Checked {report['checked_variants']} variants and {len(counts)} property/state points; {len(report['discrepancies'])} raw/processed discrepancies. No simulations run.")
    if report["discrepancies"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
