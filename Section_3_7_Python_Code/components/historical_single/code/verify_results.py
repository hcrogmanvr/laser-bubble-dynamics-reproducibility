"""Verify source identity, saved comparisons, and an independent equilibrium check."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'results'
read=lambda s:json.loads((P/(s+'.json')).read_text())
for x in read('provenance')['original_sources']:
 assert hashlib.sha256((ROOT/'source_originals'/x['file']).read_bytes()).hexdigest()==x['sha256']
for x in read('source_reproduction'):
 assert x['event_count']==577
 assert max(v['max_relative_to_peak'] for v in x['comparison'].values())<1e-10
for model,k,cap in [(2,5e7,1.5e9),(3,6e7,2e9)]:
 s=np.load(P/f'continuous_figure{model}_reference.npz')
 # Independent static force balance at the retained pressure ceiling.
 a=4*np.pi*(cap-1e5);r0=2e-6
 equilibrium=2*k*r0/(k+np.sqrt(k*k-4*a*k*r0))
 assert abs(s['R'].max()-equilibrium)*1e9<1e-5
 assert s['R'].min()>(1.5e-6 if model==2 else 1.2e-6)
 rows=[x for x in read('continuous_convergence') if x['model'].startswith(f'Figure {model}')]
 assert all(x['success'] and x['events']==0 for x in rows)
 assert rows[1]['max_radius_difference_from_finest_nm']<1e-6
for x in read('legacy_controls'):
 if x['test']=='timestep':assert x['radius_clip_fraction']==x['velocity_clip_fraction']==1.
 if x['test']=='optical_postprocessing_null':assert x['max_dynamic_state_change']==0.
audit=read('archived_coupling_audit');assert all(v==0 for v in audit['maximum_dynamic_state_differences'].values())
q=read('reconstructed_quenching_control');assert q[1]['absorbed_energy_nominal_J']<q[0]['absorbed_energy_nominal_J']
for row in read('new_dynamic_coupling_controls'):
 model=2 if row['model'].startswith('Figure 2') else 3
 s=np.load(P/f"new_coupling_figure{model}_{row['spacing_um']:g}um.npz")
 assert np.isclose(np.max(abs(s['R_on']-s['R_off']))*1e9,row['max_radius_on_off_difference_nm'],rtol=1e-12,atol=0)
print('PASS source hashes, original reproduction, independent equilibrium, refinement and saved controls')
