"""Regenerate final-force-field diagnostics from committed public GPs and data.

Use the documented GP environment. Outputs require --output-dir; no simulations run.
The scaled-parameter Jacobian below is NOT an objective Hessian or a claim of
structural identifiability. Five fitting temperatures give ten mean observations.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = None  # Set by the required --output-dir CLI argument.
sys.path.insert(0, str(REPO))
from utils.molec_class_files import esolvs

MOLS = ['MeOH', 'EG', 'Gly', 'DMF', 'DMSO', 'DEC']
PROPS = ['liq_density', 'surf_tens']
RECOMMENDED = dict(MeOH='GP-Optimized', EG='GP-Optimized', Gly='Base',
                   DMF='Base', DMSO='Base', DEC='GP-Optimized')
SOURCE = dict(liq_density=('ld_iters', 'RQ'), surf_tens=('vle_iters', 'Matern32'))


def file_sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()


def latest(mol, prop):
    iters, kernel = SOURCE[prop]
    paths = list((REPO/'Build_GPs'/'analysis'/mol/iters).glob('iter-*/gp_models.pkl'))
    path = max(paths, key=lambda p: int(p.parent.name.split('-')[1]))
    return path, kernel


def theta_vectors(mol, data):
    optdir = REPO/'Opt_ES'/'analysis'/'at_00'/mol/'ExpVal'/'opt_res'
    res = {}
    for variant, folder in [('Base', 'ms_val_no_opt'), ('GP-Optimized', 'ms_val_opt')]:
        path = optdir/folder/'raw_data.csv'
        df = pd.read_csv(path)
        vectors = df[list(data.param_names)].drop_duplicates().to_numpy()
        if len(vectors) != 1:
            raise ValueError((mol, variant, 'expected unique validation vector', len(vectors)))
        res[variant] = (vectors[0], path)
    return res


def preferred_units(theta, data):
    # Constants defines bounds_sigma in nm and bounds_epsilon in kJ/mol.
    conv = np.array([10 if n.startswith('sigma') else 1/0.008314462145468951
                     for n in data.param_names])
    return np.asarray(theta)*conv


def manuscript_vectors():
    path = Path(__file__).resolve().parent/'reference/manuscript_parameters.csv'
    table = pd.read_csv(path)
    result = {}
    for mol, data in esolvs.make_dict(MOLS).items():
        result[mol] = {}
        for variant in ['Base', 'GP-Optimized']:
            rows = table[(table.molecule == mol) & (table.variant == variant)]
            result[mol][variant] = rows.set_index('parameter').loc[list(data.param_names), 'table6_rounded_value'].to_numpy()
    return result, path


def inventory():
    data_dict = esolvs.make_dict(MOLS)
    manuscript, ms_path = manuscript_vectors()
    rows, split_rows, starts_rows, oxygen_rows, criterion_rows = [], [], [], [], []
    provenance = {ms_path.relative_to(REPO).as_posix(): file_sha(ms_path)}
    for mol, data in data_dict.items():
        bounds = np.asarray(data.param_bounds)
        free = bounds[:, 1]-bounds[:, 0] > 1e-12
        st = pd.read_csv(REPO/'Build_GPs'/'analysis'/mol/'vle_iters'/'all_results.csv')
        fit_temps = sorted(st['temperature'].unique())
        for variant, (theta, path) in theta_vectors(mol, data).items():
            provenance[path.relative_to(REPO).as_posix()] = file_sha(path)
            pref = preferred_units(theta, data)
            if not np.all(np.abs(pref-manuscript[mol][variant]) <= 0.0050001):
                raise AssertionError((mol, variant, pref, manuscript[mol][variant]))
            for j, n in enumerate(data.param_names):
                width = bounds[j, 1]-bounds[j, 0]
                z = (theta[j]-bounds[j, 0])/width if free[j] else 0.0
                rows.append(dict(molecule=mol, variant=variant, recommended=RECOMMENDED[mol]==variant,
                                 parameter=n, free=bool(free[j]), real_value=theta[j],
                                 unit='nm' if n.startswith('sigma') else 'kJ/mol',
                                 preferred_value=pref[j], lower_bound=bounds[j, 0], upper_bound=bounds[j, 1],
                                 scaled_value=z, distance_to_bound_fraction=min(z,1-z) if free[j] else np.nan,
                                 inside_bounds=bool(theta[j]>=bounds[j,0]-1e-12 and theta[j]<=bounds[j,1]+1e-12),
                                 table6_rounded_value=manuscript[mol][variant][j], source=path.relative_to(REPO).as_posix()))
                if 'O' in n:
                    oxygen_rows.append(dict(molecule=mol, variant=variant, parameter=n,
                                            value=pref[j], recommended=RECOMMENDED[mol]==variant))
        for prop in PROPS:
            model_path, _ = latest(mol, prop)
            provenance[model_path.relative_to(REPO).as_posix()] = file_sha(model_path)
            tdir = model_path.parent/'train_test_sets'
            tx = pd.read_csv(tdir/f'sim_{prop}_x_train.csv')[list(data.param_names)]
            vx = pd.read_csv(tdir/f'sim_{prop}_x_test.csv')[list(data.param_names)]
            tr = set(map(tuple, np.round(tx.to_numpy(), 12)))
            te = set(map(tuple, np.round(vx.to_numpy(), 12)))
            split_rows.append(dict(molecule=mol, property=prop, train_rows=len(tx), test_rows=len(vx),
                                   train_parameter_sets=len(tr), test_parameter_sets=len(te),
                                   shared_parameter_sets=len(tr&te), test_sets_unseen=len(te-tr),
                                   model_source=model_path.relative_to(REPO).as_posix()))
        optdir = REPO/'Opt_ES'/'analysis'/'at_00'/mol/'ExpVal'/'opt_res'
        starts = pd.read_csv(optdir/'best_per_run.csv')
        provenance[(optdir/'best_per_run.csv').relative_to(REPO).as_posix()] = file_sha(optdir/'best_per_run.csv')
        obj_col = 'Min Obj Cum.'
        obj = starts[obj_col].to_numpy()
        for j, n in enumerate(data.param_names):
            vals = starts[n+'_cum'].to_numpy()
            for threshold in [0.01, 0.05]:
                selected = vals[obj <= np.min(obj)*(1+threshold)]
                if not len(selected): continue
                conv = 10 if n.startswith('sigma') else 1/0.008314462145468951
                width = (bounds[j, 1]-bounds[j, 0])*conv
                starts_rows.append(dict(molecule=mol, parameter=n, relative_objective_threshold=threshold,
                                        selected_runs=len(selected), total_runs=len(vals), minimum_objective=float(np.min(obj)),
                                        free=bool(free[j]), minimum=float(np.min(selected)),
                                        maximum=float(np.max(selected)), spread=float(np.ptp(selected)),
                                        spread_fraction_of_bound_width=float(np.ptp(selected)/width) if free[j] else 0,
                                        successful_selected_runs=int(np.sum((obj<=np.min(obj)*(1+threshold))&(starts['Term Status'].to_numpy()==0)))))
        full = REPO/'Build_GPs'/'analysis'/mol/'ld_iters'/'mse-less10-full.csv'
        c = pd.read_csv(full)
        criterion_rows.append(dict(molecule=mol, selected_rows=len(c),
                                   minimum_rmse=float(np.sqrt(c.mse.min())), maximum_rmse=float(np.sqrt(c.mse.max())),
                                   max_mse=float(c.mse.max()), temperatures=fit_temps,
                                   n_free=int(free.sum()), n_listed=len(free)))
    for name, content in [('parameters',rows),('saved_split',split_rows),('multistart_variability',starts_rows),
                          ('oxygen_parameters',oxygen_rows),('density_criterion',criterion_rows)]:
        pd.DataFrame(content).to_csv(OUT/f'{name}.csv', index=False)
    (OUT/'inventory_provenance.json').write_text(json.dumps(provenance, indent=2))
    print(pd.DataFrame(criterion_rows).to_string(index=False), flush=True)
    print(pd.DataFrame(split_rows).to_string(index=False), flush=True)


def load_model(mol, prop):
    path, kernel = latest(mol, prop)
    with open(path,'rb') as f:
        models, labels = pickle.load(f)
    g = models[f'sim_{prop}'][kernel]
    return g, path


def svd_info(jac, names):
    _, s, vt = np.linalg.svd(jac, full_matrices=True)
    # Include exact row-count null directions so P>N is visible.
    s_padded = np.pad(s, (0, len(names)-len(s)))
    relative = s_padded/s_padded[0] if s_padded[0] else np.zeros(len(names))
    ranks = {f'rank_{tol:.0e}': int(np.sum(relative>tol)) for tol in [1e-4,1e-6,1e-8]}
    retained = s[relative[:len(s)]>1e-6]
    info = dict(rows=len(jac), free_parameters=len(names), **ranks,
                s_max=float(s[0]), s_min=float(s_padded[-1]),
                condition_full=float(s[0]/s_padded[-1]) if s_padded[-1]>0 else None,
                condition_retained_1e6=float(retained[0]/retained[-1]) if len(retained) else None)
    return info, s_padded, vt


def local_analysis(molecules):
    import tensorflow as tf
    import gpflow
    from gpflow.utilities import parameter_dict
    version = dict(python=sys.version, tensorflow=tf.__version__, gpflow=gpflow.__version__,
                   numpy=np.__version__, pandas=pd.__version__)
    (OUT/'environment.json').write_text(json.dumps(version, indent=2))
    data_dict = esolvs.make_dict(molecules)
    all_rows, prediction_rows, derivative_rows, hyper_rows = [], [], [], []
    for mol, data in data_dict.items():
        bounds = np.asarray(data.param_bounds)
        widths = bounds[:,1]-bounds[:,0]
        free = widths>1e-12
        names = np.asarray(data.param_names)[free].tolist()
        vectors = theta_vectors(mol,data)
        st = pd.read_csv(REPO/'Build_GPs'/'analysis'/mol/'vle_iters'/'all_results.csv')
        temps = np.sort(st.temperature.unique())
        blocks, aug_blocks = {k:[] for k in vectors}, {k:[] for k in vectors}
        for prop in PROPS:
            model, model_path = load_model(mol, prop)
            # Cached posterior computes the training solve once, keeping gradients
            # with respect to test inputs valid and original archived GP immutable.
            post = model.posterior(precompute_cache=gpflow.posteriors.PrecomputeCacheType.TENSOR)
            for name, val in parameter_dict(model).items():
                hyper_rows.append(dict(molecule=mol, property=prop, name=name, values=np.asarray(val).tolist(),
                                       source=model_path.relative_to(REPO).as_posix()))
            tb = data.temperature_bounds('expt_'+prop)
            pb = getattr(data,prop+'_bounds')
            exp = getattr(data,'expt_'+prop)
            rel_unc = max(float(data.uncertainty['expt_'+prop]),0.02)
            measurement_sd = np.array([rel_unc*exp[float(t)] for t in temps])
            span = float(pb[1]-pb[0])
            for variant, (theta, _) in vectors.items():
                z = np.divide(theta-bounds[:,0], widths, out=np.zeros_like(theta), where=free)
                x = np.column_stack([np.tile(z,(len(temps),1)), (temps-tb[0])/(tb[1]-tb[0])])
                xt = tf.Variable(x,dtype=tf.float64)
                with tf.GradientTape() as tape:
                    mean,var = post.predict_f(xt)
                jac_scaled = tape.gradient(mean,xt).numpy()[:, :-1][:,free]
                jac_physical_output = jac_scaled*span
                jac_weighted = jac_physical_output/measurement_sd[:,None]
                sd_gp = np.sqrt(np.maximum(var.numpy().ravel(),0))*span
                jac_aug = jac_physical_output/np.sqrt(measurement_sd**2+sd_gp**2)[:,None]
                blocks[variant].append(jac_weighted)
                aug_blocks[variant].append(jac_aug)
                label = f'{mol}_{variant.replace("-","")}_{prop}'
                pd.DataFrame(jac_weighted,index=temps,columns=names).to_csv(OUT/f'{label}_weighted_jacobian.csv',index_label='temperature')
                info,s,vt = svd_info(jac_weighted,names)
                all_rows.append(dict(molecule=mol,variant=variant,block=prop,weighting='objective_diagonal',**info))
                pd.DataFrame({'mode':np.arange(1,len(s)+1),'singular_value':s,'relative_singular_value':s/s[0]}).to_csv(OUT/f'{label}_singular_values.csv',index=False)
                pd.DataFrame(vt,columns=names).to_csv(OUT/f'{label}_right_vectors.csv',index_label='mode_0based')
                # Central differences on all free columns, at two step sizes.
                for step in [1e-3,1e-4,1e-5,1e-6]:
                    finite=np.empty_like(jac_scaled)
                    for k,j in enumerate(np.flatnonzero(free)):
                        plus=x.copy();minus=x.copy();plus[:,j]+=step;minus[:,j]-=step
                        finite[:,k]=(post.predict_f(plus)[0].numpy().ravel()-post.predict_f(minus)[0].numpy().ravel())/(2*step)
                    delta=finite-jac_scaled
                    derivative_rows.append(dict(molecule=mol,variant=variant,property=prop,step_scaled=step,
                                                max_absolute_discrepancy=float(np.max(np.abs(delta))),
                                                relative_frobenius_discrepancy=float(np.linalg.norm(delta)/max(np.linalg.norm(jac_scaled),1e-30))))
                for i,t in enumerate(temps):
                    prediction_rows.append(dict(molecule=mol,variant=variant,property=prop,temperature=float(t),
                                                target=float(exp[float(t)]),GP_mean=float(mean.numpy()[i,0]*span+pb[0]),
                                                GP_latent_sd=float(sd_gp[i]),objective_sd_floor=float(measurement_sd[i]),
                                                relative_uncertainty_used=rel_unc))
            del model,post
        for variant in vectors:
            for weighting,bs in [('objective_diagonal',blocks),('measurement_plus_GP_diagonal',aug_blocks)]:
                jac=np.vstack(bs[variant]); info,s,vt=svd_info(jac,names)
                all_rows.append(dict(molecule=mol,variant=variant,block='combined',weighting=weighting,**info))
                label=f'{mol}_{variant.replace("-","")}_combined_{weighting}'
                pd.DataFrame(jac,columns=names).to_csv(OUT/f'{label}_weighted_jacobian.csv',index=False)
                pd.DataFrame({'mode':np.arange(1,len(s)+1),'singular_value':s,'relative_singular_value':s/s[0]}).to_csv(OUT/f'{label}_singular_values.csv',index=False)
                pd.DataFrame(vt,columns=names).to_csv(OUT/f'{label}_right_vectors.csv',index_label='mode_0based')
        # Checkpoint after every molecule; restartable without losing completed data.
        pd.DataFrame(all_rows).to_csv(OUT/'local_rank_summary.csv',index=False)
        pd.DataFrame(prediction_rows).to_csv(OUT/'local_GP_predictions.csv',index=False)
        pd.DataFrame(derivative_rows).to_csv(OUT/'derivative_validation.csv',index=False)
        (OUT/'local_hyperparameters.json').write_text(json.dumps(hyper_rows,indent=2))
        print(pd.DataFrame([r for r in all_rows if r['molecule']==mol]).to_string(index=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['inventory','local'])
    parser.add_argument('--molecules',nargs='+',default=MOLS)
    parser.add_argument('--output-dir', type=Path, required=True)
    args=parser.parse_args()
    OUT = args.output_dir.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.mode=='inventory':inventory()
    else:local_analysis(args.molecules)
