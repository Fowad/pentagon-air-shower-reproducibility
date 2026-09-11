#!/usr/bin/env python3
"""Paired 8°/11° exact-source robustness check using the rebuilt event archive."""
from __future__ import annotations
import argparse, importlib.util, io, json, multiprocessing as mp, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm

V12=None; STATE=None; CFG=None; SRC_RA=None; SRC_DEC=None; SS=None; CS=None; COSCUT=None; RMAX=None

def load_module(path: Path):
    spec=importlib.util.spec_from_file_location('v12_repro_ap',path); m=importlib.util.module_from_spec(spec); sys.modules[spec.name]=m; spec.loader.exec_module(m); return m

def init_worker(v12_path, events_path, runs_path):
    global V12,STATE,CFG,SRC_RA,SRC_DEC,SS,CS,COSCUT,RMAX
    V12=load_module(Path(v12_path)); events=pd.read_csv(events_path); runs=pd.read_csv(runs_path)
    CFG=V12.Config(data_dir=Path('.'),output_dir=Path('.'),event_cache=None,run_summary_cache=None,make_figures=False)
    interval_id,_,_=V12.equal_livetime_interval_definition(events,runs,6); STATE=V12.randomization_state(events,CFG,interval_id)
    lhaaso=pd.read_csv(io.StringIO(V12.LHAASO_CATALOG_CSV)); targets=pd.concat([lhaaso[['source_name','ra_deg','dec_deg']],pd.DataFrame([{'source_name':'Crab_exact_benchmark','ra_deg':83.6331,'dec_deg':22.0145}])],ignore_index=True)
    SRC_RA=np.radians(targets.ra_deg.to_numpy(float)); SRC_DEC=np.radians(targets.dec_deg.to_numpy(float)); SS=np.sin(SRC_DEC); CS=np.cos(SRC_DEC); COSCUT=np.cos(np.radians((8.,11.))); RMAX=np.radians(11.)

def counts_for_trial(trial):
    rng=V12.realization_rng(CFG.random_seed,12001,int(trial)); n=len(STATE.hour_angle_deg); randomized_lst=np.empty(n); randomized_pre=[np.empty(n) for _ in range(3)]; randomized_interval=np.empty(n,dtype=np.int16)
    for inds in STATE.blocks:
        reassigned=inds[rng.permutation(len(inds))]; randomized_lst[inds]=STATE.observed_lst_deg[reassigned]; randomized_interval[inds]=STATE.interval_id[reassigned]
        for d,s in zip(randomized_pre,STATE.precession): d[inds]=s[reassigned]
    rr,dd=V12.mean_of_date_to_icrs_with_angles((randomized_lst-STATE.hour_angle_deg)%360.,STATE.declination_mean_deg,*randomized_pre); rr=np.radians(rr); dd=np.radians(dd); sd=np.sin(dd); cd=np.cos(dd); out=np.zeros((2,len(SRC_RA),7),dtype=np.int32)
    for i in range(len(SRC_RA)):
        ix=np.flatnonzero(np.abs(dd-SRC_DEC[i])<=RMAX); c=sd[ix]*SS[i]+cd[ix]*CS[i]*np.cos(rr[ix]-SRC_RA[i])
        for q,cc in enumerate(COSCUT):
            hit=ix[c>=cc]; out[q,i,0]=len(hit)
            if len(hit): out[q,i,1:]=np.bincount(randomized_interval[hit],minlength=6)
    return int(trial),out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,required=True); ap.add_argument('--trials',type=int,default=2000); ap.add_argument('--workers',type=int,default=min(4,mp.cpu_count())); a=ap.parse_args(); root=a.root.resolve(); v12p=root/'audit/handoff/PENTAGON_V12_NEW_WORK_HANDOFF_PACKAGE/Pentagon_Array_Directional_Analysis_V12.py'; ep=root/'audit/v12/reconstructed_events.csv.gz'; rp=root/'audit/v12/run_summary.csv'; out=root/'audit/results/aperture_8deg_vs_11deg'; out.mkdir(parents=True,exist_ok=True)
    v12=load_module(v12p); events=pd.read_csv(ep); runs=pd.read_csv(rp); cfg=v12.Config(data_dir=Path('.'),output_dir=Path('.'),event_cache=None,run_summary_cache=None,make_figures=False); interval_id,_,_=v12.equal_livetime_interval_definition(events,runs,6); lhaaso=pd.read_csv(io.StringIO(v12.LHAASO_CATALOG_CSV)); targets=pd.concat([lhaaso[['source_name','ra_deg','dec_deg']],pd.DataFrame([{'source_name':'Crab_exact_benchmark','ra_deg':83.6331,'dec_deg':22.0145}])],ignore_index=True)
    sra=np.radians(targets.ra_deg.to_numpy(float)); sdec=np.radians(targets.dec_deg.to_numpy(float)); ss=np.sin(sdec); cs=np.cos(sdec); cc=np.cos(np.radians((8.,11.))); rmax=np.radians(11.); rr=np.radians(events.ra_deg.to_numpy(float)); dd=np.radians(events.dec_deg.to_numpy(float)); sd=np.sin(dd); cd=np.cos(dd); obs=np.zeros((2,len(targets),7),dtype=np.int32)
    for i in range(len(targets)):
        ix=np.flatnonzero(np.abs(dd-sdec[i])<=rmax); c=sd[ix]*ss[i]+cd[ix]*cs[i]*np.cos(rr[ix]-sra[i])
        for q,cut in enumerate(cc):
            hit=ix[c>=cut]; obs[q,i,0]=len(hit)
            if len(hit): obs[q,i,1:]=np.bincount(interval_id[hit],minlength=6)
    null=np.empty((a.trials,2,len(targets),7),dtype=np.int32); t0=time.time(); ctx=mp.get_context('fork')
    with ctx.Pool(a.workers,initializer=init_worker,initargs=(str(v12p),str(ep),str(rp))) as pool:
        for done,(trial,arr) in enumerate(pool.imap_unordered(counts_for_trial,range(a.trials),chunksize=2),1):
            null[trial]=arr
            if done%100==0: print(f'[8v11] {done}/{a.trials} elapsed={(time.time()-t0)/60:.1f} min',flush=True)
    split=a.trials//2; ref=null[:split].astype(float); cal=null[split:].astype(float); B=ref.mean(0); sB=ref.std(0,ddof=1); Z=np.divide(obs-B,sB,out=np.full_like(B,np.nan),where=sB>0); pemp=(1+(cal>=obs[None,...]).sum(0))/(len(cal)+1); zcal=np.divide(cal-B[None,...],sB[None,...],out=np.full_like(cal,np.nan),where=sB[None,...]>0)
    np.savez_compressed(out/f'paired_null_{a.trials}.npz',counts=null,observed=obs,radii_deg=np.array([8.,11.]),source_name=targets.source_name.to_numpy(str),ra_deg=targets.ra_deg,dec_deg=targets.dec_deg)
    periods=['full']+[f'interval_{i}' for i in range(1,7)]; rows=[]
    for q,rad in enumerate((8.,11.)):
      for si,row in targets.iterrows():
       for pi,period in enumerate(periods): rows.append({'radius_deg':rad,'source':row.source_name,'period':period,'N':int(obs[q,si,pi]),'B':float(B[q,si,pi]),'sB':float(sB[q,si,pi]),'Z':float(Z[q,si,pi]),'p_empirical':float(pemp[q,si,pi]),'p_normal':float(norm.sf(Z[q,si,pi]))})
    pd.DataFrame(rows).to_csv(out/'exact_source_8v11_results.csv',index=False); fam=[]; nsrc=len(lhaaso)
    for q,rad in enumerate((8.,11.)):
        o=float(np.nanmax(Z[q,:nsrc,0])); n=np.nanmax(zcal[:,q,:nsrc,0],axis=1); fam.append({'radius_deg':rad,'family':'37_1LHAASO_full','observed_max_Z':o,'empirical_p':float((1+(n>=o).sum())/(len(n)+1))})
        o=float(np.nanmax(Z[q,:nsrc,1:])); n=np.nanmax(zcal[:,q,:nsrc,1:],axis=(1,2)); fam.append({'radius_deg':rad,'family':'37_1LHAASO_x_6_intervals','observed_max_Z':o,'empirical_p':float((1+(n>=o).sum())/(len(n)+1))})
    pd.DataFrame(fam).to_csv(out/'family_8v11_results.csv',index=False); (out/'summary.json').write_text(json.dumps({'trials':a.trials,'reference':split,'calibration':a.trials-split,'event_fingerprint':v12.event_fingerprint(events),'families':fam},indent=2)+'\n')
if __name__=='__main__': main()
