"""Rebuild DMF candidate figure without changing historical simulation data.

Parser: Build_GPs/vle_iters/project.py, std=0.05*sqrt(max Err.Est).
Hence corrected maximum error estimate (mN/m) = stored_unc**2 / 0.05.
Bars retain the original plotting multiplier 1.96; no CI coverage is asserted.
"""
from pathlib import Path
import argparse
import importlib.metadata
import ast
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
ROOT = args.repo_root.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
DATA = ROOT/'Build_GPs/analysis/DMF/vle_iters'
results = pd.read_csv(DATA/'all_results.csv', index_col=0)
pareto = pd.read_csv(DATA/'iter-1/pareto-params.csv', index_col=0)
params = [c for c in pareto if c.startswith(('sigma_', 'epsilon_'))]
if list(pareto.index) != [11, 13, 15, 17, 22]:
    raise ValueError('Unexpected DMF candidate identifiers')
temps = np.array([308.15, 313.15, 318.15, 323.15, 328.15])

# Read the exact RD2 reference dictionaries, avoiding imports with side effects.
tree = ast.parse((ROOT/'utils/molec_class_files/esolvs.py').read_text())
refs = {}
in_dmf = False
for node in tree.body:
    if isinstance(node, ast.Assign):
        names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if 'name' in names:
            try:
                in_dmf = ast.literal_eval(node.value) == 'DMF'
            except (ValueError, TypeError):
                pass
        if in_dmf:
            for key in ['expt_liq_density', 'expt_surf_tens']:
                if key in names:
                    refs[key] = ast.literal_eval(node.value)
if len(refs) != 2:
    raise ValueError('Missing DMF reference dictionaries')

plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':21,
                     'axes.linewidth':1.4, 'savefig.facecolor':'white'})
fig = plt.figure(figsize=(12, 8.4), layout='constrained')
grid = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], wspace=.12)
ax_a = fig.add_subplot(grid[0, 0])
ax_b = fig.add_subplot(grid[1, 0], sharex=ax_a)
ax_c = fig.add_subplot(grid[:, 1])
rows = []
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for k, (candidate, par) in enumerate(pareto.iterrows()):
    mask = np.isclose(results[params].to_numpy(), par[params].to_numpy(dtype=float),
                      atol=1e-12, rtol=0).all(axis=1)
    group = results.loc[mask].sort_values('temperature').copy()
    if len(group) != 5:
        raise ValueError(f'Candidate {candidate} does not have five temperatures')
    np.testing.assert_array_equal(group.temperature.to_numpy(), temps)
    corrected = 20 * group.surf_tens_unc.to_numpy()**2
    expected = {11:[1.1,1.35,.85,.95,1.05], 13:[.95,1.5,.65,1.2,.95],
                15:[.55,1.2,.85,.6,1.25], 17:[.8,1.4,.8,1,1],
                22:[.75,.8,1.1,1.25,1.3]}
    np.testing.assert_allclose(corrected, expected[candidate], atol=1e-10)
    color = colors[k]
    for ax, prop, unc in [(ax_a,'liq_density',group.liq_density_unc.to_numpy()),
                          (ax_b,'surf_tens',corrected)]:
        ax.errorbar(temps, group[prop], yerr=1.96*unc, fmt='o',
                    color=color, markersize=6, capsize=3, elinewidth=1.2,
                    alpha=.85, label=f'Candidate {candidate}')
    ax_c.scatter(par.mapd_liq_density, par.mapd_surf_tens,
                 color=color, s=70, label=f'Candidate {candidate}')
    for j, (_, row) in enumerate(group.iterrows()):
        rows.append({'candidate':candidate, 'temperature_K':row.temperature,
                     'density_kg_m3':row.liq_density,
                     'density_errorbar_halfwidth_kg_m3':1.96*row.liq_density_unc,
                     'surface_tension_mN_m':row.surf_tens,
                     'stored_transformed_uncertainty':row.surf_tens_unc,
                     'recovered_max_error_estimate_mN_m':corrected[j],
                     'surface_tension_errorbar_halfwidth_mN_m':1.96*corrected[j]})

for ax, key in [(ax_a,'expt_liq_density'), (ax_b,'expt_surf_tens')]:
    ax.scatter(temps, [refs[key][t] for t in temps], marker='x', color='black',
               linewidths=1.8, s=55, zorder=10, label='RD2 reference')
    ax.set_xticks([310,315,320,325])
    ax.set_xlim(307.1,329.2)
    ax.tick_params(direction='in', top=True, right=True)
ax_a.tick_params(labelbottom=False)
ax_a.set_ylabel(r'$\rho_l$ / kg m$^{-3}$')
ax_b.set_ylabel(r'$\gamma$ / mN m$^{-1}$')
ax_b.set_xlabel(r'$T$ / K')
ax_c.set_xlabel(r'MAPD $\rho_l$ / %')
ax_c.set_ylabel(r'MAPD $\gamma$ / %')
ax_c.tick_params(direction='in', top=True, right=True)
ax_c.set_xlim(.04,.71)
ax_c.set_ylim(-1,39)
for ax, label, x, y in [(ax_a,'A',.04,.10),(ax_b,'B',.97,.95),(ax_c,'C',.04,.10)]:
    ax.text(x,y,f'({label})',transform=ax.transAxes,va='top',ha='right' if label == 'B' else 'left',fontweight='bold',fontsize=25)
ax_a.text(.97,.95,'DMF',transform=ax_a.transAxes,va='top',ha='right',fontsize=25)
handles, labels = ax_a.get_legend_handles_labels()
fig.legend(handles, labels, loc='outside upper center',ncol=3,frameon=False,fontsize=17)
fig.savefig(OUT/'DMF_pareto_corrected.png',dpi=240)
fig.savefig(OUT/'DMF_pareto_corrected.svg')
actual = pd.DataFrame(rows)
expected_table = pd.read_csv(ROOT/'reproducibility/dmf/reference_figure_data.csv')
if list(actual.columns) != list(expected_table.columns) or actual.shape != expected_table.shape:
    raise ValueError('Plotted data schema differs from approved reference')
np.testing.assert_allclose(actual.to_numpy(), expected_table.to_numpy(), rtol=1e-12, atol=1e-12)
actual.to_csv(OUT/'corrected_figure_data.csv',index=False)
manifest = {'sources':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in
             [DATA/'all_results.csv', DATA/'iter-1/pareto-params.csv',
              ROOT/'utils/molec_class_files/esolvs.py']},
            'formula':'corrected_error_estimate_mN_m = 20 * stored_uncertainty**2',
            'bar_multiplier':1.96, 'points':len(rows),
            'environment':{name:importlib.metadata.version(name) for name in ['numpy','pandas','matplotlib']},
            'verification':'All 25 plotted rows agree with approved reference_figure_data.csv at rtol=atol=1e-12',
            'note':'Only maximum GROMACS error estimate is recoverable; not calibrated 95% CI.'}
(OUT/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Verified and regenerated all 25 DMF candidate points.')
