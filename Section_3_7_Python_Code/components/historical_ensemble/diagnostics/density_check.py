"""Observe the original Figure 6 Saha call without changing its returned values.
Run from any directory; source and archived baseline arrays remain unchanged.
"""
from pathlib import Path
import importlib.util, sys, json, hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent
PKG=ROOT/'package' if (ROOT/'package').exists() else ROOT.parent
src=PKG/'source_originals/figure6_extended_original.py'
spec=importlib.util.spec_from_file_location('figure6_density_source',src)
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
original=mod.SahaEOS.pressure
rows=[]
N_heavy=1e25*(4*np.pi/3)*(2e-7)**3
def observe(self,n_total,Ti,Te):
    pressure,x=original(self,n_total,Ti,Te)
    ne=x*n_total
    j=int(np.argmax(ne));R=(3*N_heavy/(4*np.pi*n_total[j]))**(1/3)
    rows.append([len(rows)*1e-7,float(ne.min()),float(ne.mean()),float(ne.max()),float(n_total.max()),float(x.max()),j,float(Te[j]),float(R)])
    return pressure,x
mod.SahaEOS.pressure=observe
out=mod.simulate()
baseline=np.load(PKG/'returned_arrays/figure6_returned_arrays.npz')
comparisons={}
for key in baseline.files:
    if key=='runtime_s':continue
    if key in out:
        comparisons[key]=bool(np.array_equal(baseline[key],np.asarray(out[key])))
assert comparisons and all(comparisons.values()),comparisons
values=np.array(rows)
epsilon0,me,c,e=8.8541878128e-12,9.1093837015e-31,299792458.,1.602176634e-19
nc=epsilon0*me*(2*np.pi*c/1.064e-6)**2/e**2
peak=values[int(np.argmax(values[:,3]))]
folder=ROOT/'diagnostics' if (ROOT/'package').exists() else ROOT
folder.mkdir(exist_ok=True)
np.savetxt(folder/'figure6_electron_density_history.csv',values,delimiter=',',comments='',header='grid_label_s,ne_min_m3,ne_mean_m3,ne_max_m3,nheavy_max_m3,xion_max,bubble_index_at_ne_max,Te_at_ne_max_K,R_at_ne_max_m')
report={'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'observation_stage':'Immediately after the original Saha pressure evaluation, before radial advance and event-temperature kick; grid labels preserve original timing convention.','N_heavy_per_bubble':N_heavy,'critical_density_m3':nc,'critical_density_cm3':nc/1e6,'critical_density_definition':'cold collisionless vacuum plasma-frequency equality at 1064 nm; not a dielectric/collisional breakdown or reflectivity calculation','analytic_source_ne_upper_bound_m3':1e25*(2e-7/5e-8)**3,'max_instantaneous_ne_m3':float(peak[3]),'max_instantaneous_ne_cm3':float(peak[3]/1e6),'max_ne_over_nc':float(peak[3]/nc),'max_grid_label_s':float(peak[0]),'bubble_index':int(peak[6]),'Te_at_max_ne_K':float(peak[7]),'R_at_max_ne_m':float(peak[8]),'pressure_evaluations':len(rows),'returned_arrays_bitwise_equal_to_archive':comparisons,'adiabatic_volume_ratio':(0.18/0.12)**3,'adiabatic_temperature_ratios':{'gamma_1_4':((0.18/0.12)**3)**.4,'gamma_5_3':((0.18/0.12)**3)**(2/3)},'laplace_pressure_at_0_15_um_Pa':2*.072/.15e-6}
(folder/'figure6_density_and_scale_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='returned_arrays_bitwise_equal_to_archive'},indent=2))
