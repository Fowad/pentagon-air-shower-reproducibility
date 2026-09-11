#!/usr/bin/env python3
"""Persist the Gaussian family-null maximum draws used by covariance-based source tests."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

def corr_from_cov(cov):
    s=np.sqrt(np.diag(cov)); return cov/np.outer(s,s)
def draws_max(corr,trials,seed):
    eig,vec=np.linalg.eigh((corr+corr.T)/2); fac=vec@np.diag(np.sqrt(np.clip(eig,0,None))); rng=np.random.default_rng(seed); out=np.empty(trials,dtype=np.float32)
    pos=0
    while pos<trials:
        n=min(50000,trials-pos); out[pos:pos+n]=(rng.standard_normal((n,corr.shape[0]))@fac.T).max(axis=1).astype(np.float32); pos+=n
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,required=True); ap.add_argument('--trials',type=int,default=5_000_000); a=ap.parse_args(); r=a.root.resolve()/'audit/results'
    # LHAASO 37 full archive
    d=np.load(r/'lhaaso_37_exact_full_covariance.npz',allow_pickle=False); obs=d['observed'].astype(float); mean=d['mean']; cov=d['covariance']; z=(obs-mean)/np.sqrt(np.diag(cov)); null=draws_max(corr_from_cov(cov),a.trials,370001); np.savez_compressed(r/'lhaaso_37_family_gaussian_max_5m.npz',max_z=null,observed_max_z=np.max(z),seed=np.int64(370001),trials=np.int64(a.trials))
    # 17x6 family and full 17-source family
    d=np.load(r/'p300_exact_family_covariance.npz',allow_pickle=False); cov=d['covariance']; obs=d['observed'].astype(float); mean=d['mean']; m=len(d['source']); z=(obs-mean)/np.sqrt(np.diag(cov)); interval_null=draws_max(corr_from_cov(cov),a.trials,300171)
    cov4=cov.reshape(6,m,6,m); full_cov=cov4.sum(axis=(0,2)); obs_full=obs.reshape(6,m).sum(0); mean_full=mean.reshape(6,m).sum(0); zfull=(obs_full-mean_full)/np.sqrt(np.diag(full_cov)); full_null=draws_max(corr_from_cov(full_cov),a.trials,300172)
    np.savez_compressed(r/'p300_family_gaussian_max_5m.npz',interval_max_z=interval_null,full_max_z=full_null,interval_observed_max_z=np.max(z),full_observed_max_z=np.max(zfull),interval_seed=np.int64(300171),full_seed=np.int64(300172),trials=np.int64(a.trials))
    print(json.dumps({'lhaaso_null':str(r/'lhaaso_37_family_gaussian_max_5m.npz'),'p300_null':str(r/'p300_family_gaussian_max_5m.npz'),'trials':a.trials},indent=2))
if __name__=='__main__': main()
