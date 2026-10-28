"""Additional five-fold ST-GP validation holding complete LJ vectors out.

Reads immutable archived train/test CSVs and density GPs. Fits fresh Matern32
models, with the archived mean-function form, fixed Gaussian 1e-5 variance and
training-only Eotvos noise where the archive has fixed White variance.
No held-out ST values initialize or select fitted hyperparameters.
Three deterministic initializations/fold, finite fit limit: additional validation.
"""
from __future__ import annotations
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tensorflow as tf
import gpflow
from gpflow.utilities import parameter_dict
from sklearn.linear_model import LinearRegression

from local_estimability import REPO, MOLS, latest, load_model
OUT = None  # Set by --output-dir.
sys.path.insert(0,str(REPO))
from utils.molec_class_files import esolvs
from fffit.fffit.models import MaskedLinear

SEED=42


def metrics(y, mean, variance, span):
    residual=(np.asarray(mean).ravel()-np.asarray(y).ravel())*span
    sigma=np.sqrt(np.maximum(np.asarray(variance).ravel(),0))*span
    return dict(MAE=float(np.mean(np.abs(residual))),RMSE=float(np.sqrt(np.mean(residual**2))),
                coverage95=float(np.mean(np.abs(residual)<=1.959963984540054*sigma)),
                mean_predictive_sd=float(np.mean(sigma)))


def fold_noise(x,y,data,rho_posterior,mol):
    gamma_span=float(np.diff(data.surf_tens_bounds)[0])
    gamma=np.asarray(y).ravel()*gamma_span+data.surf_tens_bounds[0]
    rho=rho_posterior.predict_f(x)[0].numpy().ravel()*float(np.diff(data.liq_density_bounds)[0])+data.liq_density_bounds[0]
    tb=data.temperature_bounds('expt_surf_tens')
    temps=x[:,-1]*float(tb[1]-tb[0])+tb[0]
    groupkeys=[tuple(np.round(row[:-1],12)) for row in x]
    mads=[]
    for group in sorted(set(groupkeys)):
        sel=np.array([g==group for g in groupkeys])
        volume=(data.molecular_weight/1000/rho[sel])**(2/3)
        t=temps[sel].reshape(-1,1)
        transformed=gamma[sel]*volume
        fit=LinearRegression().fit(t,transformed)
        predicted=fit.predict(t)/volume
        mads.append(float(np.mean(np.abs(predicted-gamma[sel]))))
    mad=float(np.mean(mads) if mol=='Gly' else np.max(mads))
    return (mad/gamma_span)**2, mads


def run(molecules,nfold,maxiter,restarts):
    np.random.seed(SEED);tf.random.set_seed(SEED)
    data_dict=esolvs.make_dict(molecules)
    rows,pred_rows=[],[]
    suffix='' if restarts==1 else f'_{restarts}starts'
    resultsdir=OUT/('grouped_validation'+suffix);resultsdir.mkdir(exist_ok=True)
    (resultsdir/'protocol.json').write_text(json.dumps(dict(seed=SEED,folds=nfold,maxiter=maxiter,
        initialization=f'{restarts} starts: ones first; independent seeded broad kernel draws next, all zero linear coefficients',
        selection='published Matern32 ST kernel, archived mean form; no held-out model selection',
        noise='fixed white Eotvos MAD recomputed using training ST groups only where original fixed; otherwise trained white initialized at 1; Gaussian fixed 1e-5',
        targets='archived scaled sim_surf_tens_y_train/y_test, recombined; no original replacement'),indent=2))
    for mol,data in data_dict.items():
        archived,path=load_model(mol,'surf_tens');tdir=path.parent/'train_test_sets'
        names=list(data.param_names)+['temperature']
        xs={s:pd.read_csv(tdir/f'sim_surf_tens_x_{s}.csv')[names].to_numpy() for s in ['train','test']}
        ys={s:pd.read_csv(tdir/f'sim_surf_tens_y_{s}.csv')['sim_surf_tens'].to_numpy().reshape(-1,1) for s in ['train','test']}
        span=float(np.diff(data.surf_tens_bounds)[0])
        am,av=archived.predict_y(xs['test'])
        rows.append(dict(molecule=mol,split='original_row_split',fold=-1,train_rows=len(xs['train']),test_rows=len(xs['test']),
                         success=True,iterations=-1,**metrics(ys['test'],am.numpy(),av.numpy(),span)))
        x=np.vstack([xs['train'],xs['test']]);y=np.vstack([ys['train'],ys['test']])
        keys=np.round(x[:,:-1],12)
        groups=np.unique(keys,axis=0)
        unique_count=len(groups)
        rng=np.random.RandomState(SEED)
        order=rng.permutation(unique_count)
        folds=np.array_split(order,nfold)
        group_ids=np.array([next(j for j,g in enumerate(groups) if np.array_equal(row,g)) for row in keys])
        pd.DataFrame(dict(row=np.arange(len(x)),group_id=group_ids)).to_csv(resultsdir/f'{mol}_group_membership.csv',index=False)
        pd.DataFrame(groups,columns=data.param_names).to_csv(resultsdir/f'{mol}_group_vectors_scaled.csv',index=False)
        rho_model,_=load_model(mol,'liq_density')
        rho_post=rho_model.posterior(precompute_cache=gpflow.posteriors.PrecomputeCacheType.TENSOR)
        fixed_white=not archived.kernel.kernels[-1].variance.trainable
        masked=isinstance(archived.mean_function,MaskedLinear)
        mask=np.asarray(archived.mean_function.mask) if masked else None
        for i,test_groups in enumerate(folds):
            if (resultsdir/f'{mol}_fold{i}_metrics.json').exists():
                raise FileExistsError('Choose a fresh output directory; fold outputs already exist')
            test=np.isin(group_ids,test_groups);train=~test
            assert not set(group_ids[train])&set(group_ids[test])
            noise,mads=fold_noise(x[train],y[train],data,rho_post,mol) if fixed_white else (1.,[])
            start=time.monotonic()
            candidates=[];restart_receipts=[]
            for restart in range(restarts):
                irng=np.random.RandomState(restart)
                lengthscales=np.ones(x.shape[1]) if restart==0 else irng.uniform(.1,100.,x.shape[1])
                variance=1. if restart==0 else irng.lognormal(0.,1.)
                linear=gpflow.mean_functions.Linear(A=np.zeros((x.shape[1],1)))
                mean=MaskedLinear(linear,mask) if masked else linear
                kernel=gpflow.kernels.Matern32(variance=variance,lengthscales=lengthscales)+gpflow.kernels.White(variance=max(noise,1e-12))
                gpflow.utilities.set_trainable(kernel.kernels[-1].variance,not fixed_white)
                candidate=gpflow.models.GPR((x[train],y[train]),kernel,mean_function=mean,noise_variance=1e-5)
                gpflow.utilities.set_trainable(candidate.likelihood.variance,False)
                fit=gpflow.optimizers.Scipy().minimize(candidate.training_loss,candidate.trainable_variables,
                       method='L-BFGS-B',options=dict(maxiter=maxiter,ftol=1e-9,maxls=50))
                loss=float(candidate.training_loss().numpy())
                candidates.append((candidate,fit,loss,restart))
                restart_receipts.append(dict(restart=restart,loss=loss,success=bool(fit.success),iterations=int(fit.nit),message=str(fit.message)))
            eligible=[c for c in candidates if c[1].success] or candidates
            model,opt,_,selected_restart=min(eligible,key=lambda c:c[2])
            pm,pv=model.predict_y(x[test]);elapsed=time.monotonic()-start
            result=dict(molecule=mol,split='grouped_set_holdout',fold=i,train_rows=int(train.sum()),test_rows=int(test.sum()),
                        train_sets=int(len(set(group_ids[train]))),test_sets=int(len(test_groups)),shared_sets=0,
                        success=bool(opt.success),status=int(opt.status),message=str(opt.message),iterations=int(opt.nit),
                        seconds=elapsed,training_loss=float(model.training_loss().numpy()),selected_restart=selected_restart,restarts=restarts,
                        fixed_white=fixed_white,white_noise_variance_scaled=float(model.kernel.kernels[-1].variance.numpy()),
                        **metrics(y[test],pm.numpy(),pv.numpy(),span))
            (resultsdir/f'{mol}_fold{i}_metrics.json').write_text(json.dumps(result,indent=2))
            (resultsdir/f'{mol}_fold{i}_parameters.json').write_text(json.dumps({k:np.asarray(v).tolist() for k,v in parameter_dict(model).items()},indent=2))
            (resultsdir/f'{mol}_fold{i}_noise_mads.json').write_text(json.dumps(mads))
            (resultsdir/f'{mol}_fold{i}_restart_receipts.json').write_text(json.dumps(restart_receipts,indent=2))
            with open(resultsdir/f'{mol}_fold{i}_model.pkl','wb') as f:pickle.dump(model,f)
            pr=[]
            for j,row in enumerate(np.flatnonzero(test)):
                pr.append(dict(molecule=mol,fold=i,row=int(row),group_id=int(group_ids[row]),
                               scaled_target=float(y[row,0]),scaled_mean=float(pm.numpy()[j,0]),scaled_predictive_variance=float(pv.numpy()[j,0]),
                               physical_target=float(y[row,0]*span+data.surf_tens_bounds[0]),
                               physical_mean=float(pm.numpy()[j,0]*span+data.surf_tens_bounds[0]),
                               physical_predictive_sd=float(np.sqrt(pv.numpy()[j,0])*span)))
            pd.DataFrame(pr).to_csv(resultsdir/f'{mol}_fold{i}_predictions.csv',index=False)
            pred_rows.extend(pr);rows.append(result)
            pd.DataFrame(rows).to_csv(OUT/f'grouped_validation_metrics{suffix}.csv',index=False)
            pd.DataFrame(pred_rows).to_csv(OUT/f'grouped_validation_predictions{suffix}.csv',index=False)
            print(result,flush=True)
        del archived,rho_model,rho_post
    pd.DataFrame(rows).to_csv(OUT/f'grouped_validation_metrics{suffix}.csv',index=False)
    pd.DataFrame(pred_rows).to_csv(OUT/f'grouped_validation_predictions{suffix}.csv',index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--molecules',nargs='+',default=MOLS)
    p.add_argument('--folds',type=int,default=5);p.add_argument('--maxiter',type=int,default=2000)
    p.add_argument('--restarts',type=int,default=3)
    p.add_argument('--output-dir', type=Path, required=True)
    args=p.parse_args()
    OUT = args.output_dir.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    run(args.molecules,args.folds,args.maxiter,args.restarts)
