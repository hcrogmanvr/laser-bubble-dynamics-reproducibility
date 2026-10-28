"""Export the arrays returned by the unmodified supplied Figure 5–6 scripts.

This is a source-execution check, not a corrected physical simulation or a
claim that the existing manuscript images are reproduced pixel for pixel.
"""
from pathlib import Path
import importlib.util,sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'returned_arrays';OUT.mkdir(exist_ok=True)
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
result=[]
for name,file,entry in [('figure5','figure5_two_reservoir_original.py','run_multibubble_2T_plasma'),('figure6','figure6_extended_original.py','simulate')]:
 src=ROOT/'source_originals'/file;m=module(name,src)
 settings=dict(N=8,t_end=1e-3,dt=1e-7,stride=5,seed=2)
 if name=='figure6':settings['cluster_sigma']=2e-6
 out=getattr(m,entry)(**settings)
 arrays={k:v for k,v in out.items() if isinstance(v,np.ndarray)}
 np.savez_compressed(OUT/(name+'_returned_arrays.npz'),**arrays)
 scalars={k:(v.item() if isinstance(v,np.generic) else v) for k,v in out.items() if not isinstance(v,np.ndarray)}
 metrics=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),settings=settings,returned_scalars=scalars,all_returned_arrays_finite=all(np.isfinite(v).all() for v in arrays.values()),stored_samples=len(out['t']),maximum_stored_mean_Te_K=float(np.max(out['Te_mean'])),minimum_stored_mean_radius_m=float(np.min(out['R_mean'])),maximum_stored_mean_radius_m=float(np.max(out['R_mean'])))
 if name=='figure6':metrics['final_ruptured_fraction']=float(np.mean(out['final_ruptured']))
 (OUT/(name+'_execution_summary.json')).write_text(json.dumps(metrics,indent=2))
 result.append({'source':name,**metrics})
(OUT/'execution_checks.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
