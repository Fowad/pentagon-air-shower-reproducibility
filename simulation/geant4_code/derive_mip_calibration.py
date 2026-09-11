#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import moyal

def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--input',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args(); d=np.load(a.input,allow_pickle=False); x=d['deposited_energy_mev'].astype(float); loc,scale=moyal.fit(x); out={'input':str(a.input),'input_sha256':sha(a.input),'n_events':len(x),'moyal_mpv_mev':float(loc),'moyal_scale_mev':float(scale),'mean_mev':float(x.mean()),'median_mev':float(np.median(x)),'published_fold_reference_mip_mev':3.57852,'absolute_difference_from_published_mev':float(abs(loc-3.57852))}; a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2));
 if abs(loc-3.57852)>0.03: raise SystemExit('Regenerated MIP MPV differs unexpectedly from the frozen published response scale.')
if __name__=='__main__': main()
