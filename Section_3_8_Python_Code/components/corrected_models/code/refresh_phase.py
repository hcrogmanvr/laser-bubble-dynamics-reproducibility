"""Recompute phase cases without repeating completed ensemble controls."""
from run_suite import *
if __name__=='__main__':
    phases=[];sols=[]
    for tol,step in [(1e-6,.08),(1e-7,.04),(1e-8,.02)]:
        o,s=run_phase(PhaseParameters(rtol=tol,max_step_us=step));s['test']='refinement';phases.append(s);sols.append((o,s))
    fine,sm=sols[-1];np.savez_compressed(OUT/'phase_on.npz',**fine)
    conv=[]
    for o,s in sols:
        conv.append(dict(max_step_us=s['parameters']['max_step_us'],rtol=s['parameters']['rtol'],final_radius_error_nm=float(abs(o['R_m'][-1]-fine['R_m'][-1])*1e9),final_vapor_mass_relative_error=abs(s['vapor_mass_ratio_final']/sm['vapor_mass_ratio_final']-1),peak_wall_temperature_difference_K=abs(s['max_Twall_K']-sm['max_Twall_K']),energy_residual_relative=s['energy_residual_relative']))
    for name,kw in [('phase_off',dict(phase_on=False)),('accommodation_low',dict(accommodation=1e-4)),('accommodation_high',dict(accommodation=1e-2)),('absorption_low',dict(absorption_fraction=1e-5)),('absorption_high',dict(absorption_fraction=1e-3))]:
        o,s=run_phase(PhaseParameters(**kw));s['test']=name;phases.append(s);np.savez_compressed(OUT/f'{name}.npz',**o)
    dump('phase_controls',phases);dump('phase_convergence',conv)
    assert max(s['energy_residual_relative'] for s in phases)<1e-6
    print(json.dumps(phases,indent=2))
