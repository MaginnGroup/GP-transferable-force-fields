"""Check exported LJ values, charges, scaling and archived MDP settings.

This is a file consistency check; it does not invoke simulation software.
"""
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
R = 0.00831446261815324


def section_lines(text, name):
    current = ""
    for line in text.splitlines():
        if line.lstrip().startswith("["):
            current = line.strip().strip("[] ")
        elif current == name and line.strip() and not line.lstrip().startswith(";"):
            yield line.split(";")[0].split()


def mdp(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split(";", 1)[0]
        if "=" in line:
            name, value = line.split("=", 1)
            result[name.strip().lower()] = value.strip().lower()
    return result


def main():
    count = 0
    manifest = json.loads((ROOT / "Opt_ES/FF_params/manifest.json").read_text())
    for entry in manifest:
        directory = ROOT / entry["directory"]
        xml = ET.parse(directory / "ff.xml").getroot()
        atoms = {atom.attrib["type"]: atom.attrib for atom in xml.findall(".//NonbondedForce/Atom")}
        vector = entry["parameters_nm_kJ_per_mol"]
        for name, atom in atoms.items():
            for parameter in ("sigma", "epsilon"):
                assert abs(float(atom[parameter]) - vector[parameter + "_" + name]) <= 5.1e-7, (entry, name, parameter)
                count += 1
        topology = (directory / "system.top").read_text()
        for cells in section_lines(topology, "atomtypes"):
            atom = atoms[cells[0]]
            assert abs(float(cells[-2]) - float(atom["sigma"])) <= 1e-10
            assert abs(float(cells[-1]) - float(atom["epsilon"])) <= 1e-10
        charge = 0
        for cells in section_lines(topology, "atoms"):
            assert abs(float(cells[6]) - float(atoms[cells[1]]["charge"])) <= 1e-10
            charge += float(cells[6])
        assert abs(charge) <= 3e-6, (entry, "net charge", charge)
        defaults = next(section_lines(topology, "defaults"))
        assert defaults[0:3] == ["1", "2", "yes"]
        assert abs(float(defaults[3]) - 0.5) < 1e-10
        assert abs(float(defaults[4]) - 5/6) < 1e-6
        mcf = (directory / "species1.mcf").read_text()
        info = mcf.split("# Atom_Info", 1)[1].split("# Bond_Info", 1)[0]
        for line in info.splitlines():
            cells = line.split()
            if len(cells) >= 8 and cells[5] == "LJ":
                atom = atoms[cells[1].removesuffix("_s1")]
                assert abs(float(cells[4]) - float(atom["charge"])) < 1e-10
                assert abs(float(cells[6]) * R - float(atom["epsilon"])) < 1e-6
                assert abs(float(cells[7]) / 10 - float(atom["sigma"])) < 1e-6
        scales = mcf.split("# Intra_Scaling", 1)[1].strip().splitlines()
        lj_scale = [float(v) for v in scales[0].split()]
        assert lj_scale == [0, 0, 0.5, 1] or (entry["molecule"] == "MeOH" and lj_scale == [0, 0, 0, 1])
        if entry["molecule"] == "MeOH":
            # Every 1-4 endpoint pair in methanol contains hydroxyl H2,
            # whose LJ epsilon is zero. Thus the archived zero 1-4 LJ
            # scaling and XML/topology scaling 0.5 yield the same energy.
            assert float(atoms["H2"]["epsilon"]) == 0
            dihedrals = mcf.split("# Dihedral_Info", 1)[1].split("# Improper_Info", 1)[0]
            atom_rows = {int(line.split()[0]): line.split()[1].removesuffix("_s1") for line in info.splitlines() if "LJ" in line}
            for line in dihedrals.splitlines():
                cells = line.split()
                if len(cells) > 5 and cells[5] == "CHARMM":
                    assert "H2" in (atom_rows[int(cells[1])], atom_rows[int(cells[4])])
        assert abs(float(scales[1].split()[2]) - 5/6) < 4e-5
    for example in (ROOT / "reproducibility/examples").iterdir():
        if not example.is_dir():
            continue
        for stage in ("inter_eq", "inter_prod"):
            settings = mdp(example / stage / (stage + ".mdp"))
            assert settings["coulombtype"] == "pme"
            assert settings["vdwtype"] == "cut-off"
            assert settings["dispcorr"] == "no"
            assert float(settings["rvdw"]) == float(settings["rcoulomb"]) == 2.4
            processed = mdp(example / stage / "mdout.mdp")
            assert processed["ewald-geometry"] == "3d"
    print(f"Verified {len(manifest)} exports, {count} LJ values, neutral charges, topology/MCF unit conversion, 1-4 scaling and two archived MDP examples. No simulations run.")


if __name__ == "__main__":
    main()
