"""Explicit scalar calculations cited in the manuscript; no fitted flash data."""
from pathlib import Path
import json,csv,math,importlib.util,sys
import numpy as np
from scipy.constants import epsilon_0,m_e,e,c

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/calculations';OUT.mkdir(parents=True,exist_ok=True)
def critical_density(wavelength_m):return epsilon_0*m_e*(2*math.pi*c/wavelength_m)**2/e**2
def ideal_gas_compression(initial_radius_m,final_radius_m,gamma):
    return (initial_radius_m/final_radius_m)**(3*(gamma-1))
def shell_volume(particle_radius_m,penetration_m):return 4*math.pi*particle_radius_m**2*penetration_m
def absorption_only(optical_depth):return math.exp(-optical_depth),-math.expm1(-optical_depth)
def diameter_squared_slope_to_radius_squared(slope):return slope/4
def rayleigh_scale(radius_m,rho,pressure_difference_Pa):return .9146813565*radius_m*math.sqrt(rho/pressure_difference_Pa)
def load(path,name):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sys.modules[name]=m;sp.loader.exec_module(m);return m
def main():
    nc=critical_density(1064e-9);rows=[]
    def add(name,value,unit,source):rows.append(dict(quantity=name,value=float(value),unit=unit,definition_or_source=source))
    add('cold_vacuum_critical_density_1064nm',nc,'m^-3','Eq. 16m; vacuum plasma-frequency scale')
    add('cold_vacuum_critical_density_1064nm_cm3',nc*1e-6,'cm^-3','Eq. 16m')
    nmax=1e25*(.2/.05)**3
    add('historical_inventory_full_ionization_density_bound',nmax,'m^-3','n_ref=1e25 m^-3 at R_ref=0.2 µm; R_min=0.05 µm')
    add('historical_full_ionization_density_bound_over_nc',nmax/nc,'1','Inventory bound; not a calculated or measured ionization state')
    for gamma in [1.4,5/3]:
        factor=ideal_gas_compression(.18e-6,.12e-6,gamma)
        add(f'ideal_gas_temperature_factor_gamma_{gamma:g}',factor,'1','Illustrative fixed-mass ideal gas; Section 4.3')
        add(f'ideal_gas_final_temperature_gamma_{gamma:g}',300*factor,'K','Illustration starting at 300 K, not a source EOS output')
    for model,mass,damping in [(2,1e-12,3e4),(3,4*math.pi/3*1000*(2e-6)**3,4e4)]:
        add(f'source{model}_damping_time',mass/damping,'s','m_eff / Gamma; local explicit damping diagnostic')
        add(f'source{model}_isolated_explicit_damping_dt_limit',2*mass/damping,'s','2 m_eff / Gamma; not a full nonlinear stability theorem')
    add('empty_cavity_Rayleigh_time_30um',rayleigh_scale(30e-6,1000,101325)*1e6,'µs','Full-collapse empty-cavity scale; finite-radius benchmark is separate')
    m=load(ROOT/'components/corrected_models/code/corrected_ensemble.py','scalar_corrected_model')
    with np.load(ROOT/'components/corrected_models/results/corrected_figure6.npz') as z:
        ne=m.NPART*z['x']/(4*np.pi/3*z['R']**3)
        add('corrected_reference_formal_max_ne',ne.max(),'m^-3','Paired per-bubble Saha fraction and density; effective fractional population')
        add('corrected_reference_formal_max_ne_over_nc',ne.max()/nc,'1','No resolved plasma shielding in the cold reference run')
    fit=ROOT/'outputs/plotted_data/growth_fit_summary.csv'
    if not fit.exists():
        from plot_paper_figures import growth_analysis
        growth_analysis()
    import pandas as pd
    for r in pd.read_csv(fit).to_dict('records'):
        add(f'{r["trace"]}_diameter_squared_growth',r['slope_um2_s'],'µm²/s','Measured diameter initial-expansion OLS; descriptive only')
        add(f'{r["trace"]}_radius_squared_equivalent',diameter_squared_slope_to_radius_squared(r['slope_um2_s']),'µm²/s','Unit convention conversion D=2R; not a separate fit')
    (OUT/'reference_quantities.json').write_text(json.dumps(rows,indent=2))
    with (OUT/'reference_quantities.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print('Calculated',len(rows),'reference quantities; unassigned conceptual parameters are not invented.')
if __name__=='__main__':main()
