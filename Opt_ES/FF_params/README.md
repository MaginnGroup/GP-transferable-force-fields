# Force-field files for both fitted models

Each molecule has separate `Base/` and `GP-Optimized/` directories. The recommended variants in the manuscript are **GP-Optimized for MeOH, EG and DEC** and **Base for Gly, DMF and DMSO**. These names refer to the six individual molecule models; atom labels do not imply transferable parameters between molecules.

| File | Role | LJ and charge units |
|---|---|---|
| `ff.xml` | Foyer force field: atom definitions, bonded terms, RESP charges and LJ values | sigma: nm; epsilon: kJ/mol; charge: elementary charge |
| `system.top` | Standalone GROMACS topology for an archived molecule count | sigma: nm; epsilon: kJ/mol; charge: elementary charge |
| `species1.mcf` | Cassandra molecular connectivity, fixed bond lengths, charges and LJ values | sigma: angstrom; epsilon/kB: K; charge: elementary charge |
| `species1.pdb` | Single-molecule coordinate/atom-order reference for Cassandra | Coordinates: angstrom |

The `GP-Optimized` files come from the saved `Opt_ES/FF_params` exports (commit `a41e279`). All six original exports correspond to GP-Optimized, including molecules for which Base is recommended. The `Base` files are reconstructed from the invariant parameter vector in each public `ms_val_no_opt/ms_data.csv`; only the LJ values are changed. They are explicitly derived parameter files, rather than separately recovered original simulation input files. The bonded terms, RESP charges, atom ordering and molecule counts are preserved. The single-molecule PDB is copied because LJ optimization did not change connectivity. `manifest.json` identifies the exact archived source and unrounded parameter vector for every variant.

GROMACS `[defaults]` uses combination rule 2 (arithmetic sigma and geometric epsilon), LJ 1-4 scaling 0.5 and charge 1-4 scaling 0.8333333. The XML and Cassandra files store approximately 5/6 charge scaling with their original serialization precision (usually 0.8333). Cassandra's `Intra_Scaling` excludes 1-2 and 1-3 nonbonded terms and leaves more distant intramolecular terms unscaled. Its LJ 1-4 factor is 0.5 except in the archived methanol MCF, where it is zero. All three methanol 1-4 pairs contain hydroxyl H2 with epsilon=0, so this exception yields the same zero LJ contribution as scaling by 0.5. The input checker verifies those endpoints rather than silently changing the archive. Cassandra mixing is selected by `mixing_rule='lb'` in the project generator.

The XML exports for DMF and DEC originally contained unescaped ampersands in SMARTS definitions. Those characters are encoded as `&amp;` here so the documents are valid XML; the parsed SMARTS and physical parameters are unchanged. GROMACS and MCF files retain their original formatting/precision for GP-Optimized. Base reconstruction uses the XML writer's six-decimal sigma/epsilon serialization; MCF conversion uses R=0.00831446261815324 kJ mol^-1 K^-1. Original archived MCFs may reflect slightly different library constants within serialization precision. Use the full-precision archived vector when comparing with the displayed parameter table.

`system.top` includes an archived bulk system molecule count. Before using it with a different coordinate configuration, make `[molecules]` match that configuration. A single-molecule PDB does not supply a bulk/slab configuration, Cassandra fragment library, restart or trajectory. See [reproducibility instructions](../../reproducibility/README.md) for stage inputs, their provenance and missing assets.

From the repository root, validate the package without simulation software:

```sh
python3 reproducibility/check_inputs.py
```

To recreate the Base exports and manifest from the included optimized templates and archived vectors:

```sh
python3 reproducibility/export_force_fields.py
python3 reproducibility/check_inputs.py
```
