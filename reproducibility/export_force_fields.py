"""Export both archived parameter vectors without running simulations.

The GP-Optimized export supplies bonded topology and charges. Base changes
only LJ values, using the archived validation parameter vector. Existing
simulation inputs, trajectories and results are never changed.
"""
import csv
import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
R_KJ_MOL_K = 0.00831446261815324
MOLECULES = ("MeOH", "EG", "Gly", "DMF", "DMSO", "DEC")


def parameter_vector(molecule, variant):
    kind = "ms_val_opt" if variant == "GP-Optimized" else "ms_val_no_opt"
    path = ROOT / f"Opt_ES/analysis/at_00/{molecule}/ExpVal/opt_res/{kind}/ms_data.csv"
    rows = list(csv.DictReader(path.open()))
    columns = [key for key in rows[0] if key.startswith(("sigma_", "epsilon_"))]
    vectors = {tuple(float(row[key]) for key in columns) for row in rows}
    if len(vectors) != 1:
        raise ValueError(f"{path}: validation rows do not share one parameter vector")
    return dict(zip(columns, vectors.pop())), path.relative_to(ROOT).as_posix()


def replace_topology(text, vector):
    lines = text.splitlines(keepends=True)
    section = ""
    for index, line in enumerate(lines):
        if line.lstrip().startswith("["):
            section = line.strip().strip("[] ")
        elif section == "atomtypes" and line.strip() and not line.lstrip().startswith(";"):
            cells = line.split()
            atom = cells[0]
            cells[-2:] = [f"{vector['sigma_' + atom]:.6f}", f"{vector['epsilon_' + atom]:.6f}"]
            lines[index] = "  ".join(cells) + "\n"
    return "".join(lines)


def replace_mcf(text, vector):
    lines = text.splitlines(keepends=True)
    section = ""
    for index, line in enumerate(lines):
        if line.startswith("#"):
            section = line.strip()
        elif section == "# Atom_Info" and " LJ " in line:
            cells = line.split()
            atom = cells[1].removesuffix("_s1")
            # Foyer XML is serialized to six decimal places before conversion.
            sigma = float(f"{vector['sigma_' + atom]:.6f}")
            epsilon = float(f"{vector['epsilon_' + atom]:.6f}")
            cells[6:] = [f"{epsilon / R_KJ_MOL_K:.5f}", f"{sigma * 10:.5f}"]
            lines[index] = "    ".join(cells) + "\n"
    return "".join(lines)


def main():
    inventory = []
    for molecule in MOLECULES:
        optimized = ROOT / "Opt_ES/FF_params" / molecule / "GP-Optimized"
        original_xml = (optimized / "ff.xml").read_text()
        tree = ET.fromstring(original_xml)
        base, source = parameter_vector(molecule, "Base")
        for atom in tree.findall(".//NonbondedForce/Atom"):
            name = atom.attrib["type"]
            for property_name in ("sigma", "epsilon"):
                atom.set(property_name, f"{base[property_name + '_' + name]:.6f}")
        destination = optimized.parent / "Base"
        destination.mkdir(exist_ok=True)
        (destination / "ff.xml").write_text(ET.tostring(tree, encoding="unicode") + "\n")
        (destination / "system.top").write_text(replace_topology((optimized / "system.top").read_text(), base))
        (destination / "species1.mcf").write_text(replace_mcf((optimized / "species1.mcf").read_text(), base))
        (destination / "species1.pdb").write_bytes((optimized / "species1.pdb").read_bytes())
        for variant in ("Base", "GP-Optimized"):
            vector, source = parameter_vector(molecule, variant)
            inventory.append({"molecule": molecule, "variant": variant,
                              "recommended": variant == ("GP-Optimized" if molecule in ("MeOH", "EG", "DEC") else "Base"),
                              "parameter_source": source, "parameters_nm_kJ_per_mol": vector,
                              "export_origin": "reconstructed LJ-only Base export" if variant == "Base" else "archived GP-Optimized export from ES-FFO a41e279",
                              "directory": f"Opt_ES/FF_params/{molecule}/{variant}"})
    (ROOT / "Opt_ES/FF_params/manifest.json").write_text(json.dumps(inventory, indent=2) + "\n")


if __name__ == "__main__":
    main()
