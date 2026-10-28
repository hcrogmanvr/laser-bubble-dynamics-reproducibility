from pathlib import Path
import json,csv,time,sys,platform
from dataclasses import replace
import numpy as np
from scipy.integrate import solve_ivp
from corrected_ensemble import Parameters,run,saha,drude,evaluate,apply_event,geometry,NPART,EV
from phase_change import PhaseParameters,run_phase,rayleigh_benchmark,saturation_pressure
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';OUT.mkdir(exist_ok=True)

def dump(name,rows):
    (OUT/(name+'.json')).write_text(json.dumps(rows,indent=2))
    flat=[]
    for x in rows:
        flat.append({**x.get('parameters',{}),**{k:v for k,v in x.items() if k!='parameters'}})
    keys=list(dict.fromkeys(k for r in flat for k in r))
    with (OUT/(name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(flat)

def main():
    started=time.time();checks=[]
    for tol in [1e-5,1e-7,1e-9]:checks.append({'test':'Rayleigh finite-radius benchmark',**rayleigh_benchmark(tol)})
    # Saha x=1/2 iff A=1/2; verify a broad stable-root grid against high precision.
    from decimal import Decimal as D,getcontext
    getcontext().prec=70
    maxerr=0.
    for n in np.logspace(18,29,12):
        for T in np.logspace(2.5,6,24):
            x=saha(n,T,13.5984)[0]
            z=2*D('3.1415926535897932384626433832795028841971693993751')*D('9.1093837015e-31')*D('1.380649e-23')*D(str(T))/D('6.62607015e-34')**2
            A=2/D(str(n))*z*z.sqrt()*(-D('13.5984')*D('1.602176634e-19')/(D('1.380649e-23')*D(str(T)))).exp()
            exact=2*A/(A+(A*A+4*A).sqrt())
            maxerr=max(maxerr,abs(float(x)-float(exact)))
    checks.append(dict(test='Hydrogen Saha algebra against 70-digit arithmetic',max_absolute_fraction_error=maxerr,species='hydrogen benchmark only; ensemble retains effective 12.6 eV'))
    ne=np.r_[0.,np.logspace(10,29,200)];R=np.full_like(ne,2e-7)
    rr,tt,aa=drude(ne,R)
    checks.append(dict(test='Drude optical budget',max_fraction_sum_error=float(abs(rr+tt+aa-1).max()),min_fraction=float(min(rr.min(),tt.min(),aa.min())),vacuum_R=float(rr[0]),vacuum_T=float(tt[0]),vacuum_A=float(aa[0])))
    # Exact conservative finite-temperature exchange at fixed radius, no heat input/loss.
    p=Parameters(variant=5,N=1,power_multiplier=0,tau_e=float('inf'),tau_i=float('inf'),coupling=False)
    _,r,_,w=geometry(p);y=np.r_[r*1e6,0.,10000.,300.,0.,np.zeros(5)]
    def thermal(t,z):
        zz=y.copy();zz[2:4]=z;return evaluate(zz,p,w,0)[0][2:4]
    sol=solve_ivp(thermal,(0,20),y[2:4],method='Radau',rtol=1e-10,atol=1e-9)
    Tbar=(p.Ce*10000+p.Ci*300)/(p.Ce+p.Ci)
    delta=9700*np.exp(-(1+p.Ce/p.Ci)*20e-6/p.tau_ei)
    exact=np.array([Tbar+p.Ci/(p.Ce+p.Ci)*delta,Tbar-p.Ce/(p.Ce+p.Ci)*delta])
    checks.append(dict(test='Two-reservoir exchange exact solution',max_temperature_error_K=float(abs(sol.y[:,-1]-exact).max()),relative_internal_energy_drift=float(abs(p.Ce*sol.y[0]+p.Ci*sol.y[1]-(p.Ce*10000+p.Ci*300)).max()/(p.Ce*10000+p.Ci*300))))
    p=Parameters(N=1,variant=6,coupling=False);_,r,_,w=geometry(p)
    y=np.r_[r*1e6,-10.,20000.,5000.,0.,np.zeros(5)]
    _,loss,res=apply_event(y,0,p,w)
    checks.append(dict(test='Hot-state event energy removal',removed_J=loss,absolute_energy_residual_J=res,relative_residual=abs(res)/loss))
    checks.append(dict(test='IAPWS saturation endpoints',p_triple_Pa=float(saturation_pressure(273.16)),p_critical_Pa=float(saturation_pressure(647.096))))
    dump('closure_checks',checks)
    rows=[];traces={};convergence=[]
    for v in [5,6]:
        sols=[]
        for j,(tol,step) in enumerate([(1e-6,1.),(1e-7,.5),(1e-8,.25)]):
            p=Parameters(variant=v,rtol=tol,max_step_us=step)
            o,s=run(p);s['test']='refinement';s['name']=f'figure{v}_refine{j}';rows.append(s);sols.append((o,s))
            print(s['name'],s['events'],s['peak_Te_K'],s['energy_residual_relative'],flush=True)
        fine,sm=sols[-1];traces[v]=fine;np.savez_compressed(OUT/f'corrected_figure{v}.npz',**fine)
        # Compare final states and integrated outputs; do not compare across event jumps by interpolation.
        for o,s in sols:
            convergence.append(dict(variant=v,rtol=s['parameters']['rtol'],max_step_us=s['parameters']['max_step_us'],events=s['events'],final_radius_max_error_nm=float(abs(o['R'][-1]-fine['R'][-1]).max()*1e9),final_Te_max_error_K=float(abs(o['Te'][-1]-fine['Te'][-1]).max()),absorbed_energy_relative_error=abs(s['absorbed_J']/sm['absorbed_J']-1),max_event_time_difference_ns=float(abs(o['events'][:,0]-fine['events'][:,0]).max()*1e3) if len(o['events'])==len(fine['events']) and len(o['events']) else None,energy_residual_relative=s['energy_residual_relative']))
        cases=[('coupling_off',dict(coupling=False)),('one_temperature',dict(one_temperature=True)),('drude_off',dict(drude_on=False)),('forcing_half',dict(power_multiplier=.5)),('forcing_1p5',dict(power_multiplier=1.5)),('spacing_half',dict(spacing_multiplier=.5)),('spacing_double',dict(spacing_multiplier=2)),('radius_guard_10nm',dict(cutoff_um=.01)),('radius_guard_20nm',dict(cutoff_um=.02)),('mach_guard_0p075',dict(mach_limit=.075))]
        for name,kw in cases:
            p=Parameters(variant=v,max_step_us=.5,rtol=1e-7,**kw);o,s=run(p);s['test']=name;s['name']=f'figure{v}_{name}';rows.append(s)
            if name in ['coupling_off','one_temperature','drude_off']:np.savez_compressed(OUT/f'{s["name"]}.npz',**o)
            print(s['name'],s['events'],s['peak_Te_K'],flush=True)
    for c in convergence:
        candidate=next(s for s in rows if s['test']=='refinement' and s['parameters']['variant']==c['variant'] and s['parameters']['rtol']==c['rtol'])
        fine_summary=next(s for s in rows if s['test']=='refinement' and s['parameters']['variant']==c['variant'] and s['parameters']['rtol']==1e-8)
        c['sampled_minimum_radius_difference_nm']=abs(candidate['min_radius_um']-fine_summary['min_radius_um'])*1e3
        c['peak_Te_difference_K']=abs(candidate['peak_Te_K']-fine_summary['peak_Te_K'])
        c['sampled_peak_speed_relative_difference']=abs(candidate['peak_speed_m_s']/fine_summary['peak_speed_m_s']-1)
        c['event_removed_energy_relative_difference']=abs(candidate['removed_J']/fine_summary['removed_J']-1)
    dump('corrected_ensemble_controls',rows);dump('corrected_ensemble_convergence',convergence)
    phases=[];psols=[]
    for tol,step in [(1e-6,.08),(1e-7,.04),(1e-8,.02)]:
        o,s=run_phase(PhaseParameters(rtol=tol,max_step_us=step));s['test']='refinement';phases.append(s);psols.append((o,s));print('phase refinement',step,s['energy_residual_relative'],flush=True)
    fine,sm=psols[-1];np.savez_compressed(OUT/'phase_on.npz',**fine)
    conv=[]
    for o,s in psols:
        conv.append(dict(max_step_us=s['parameters']['max_step_us'],rtol=s['parameters']['rtol'],final_radius_error_nm=float(abs(o['R_m'][-1]-fine['R_m'][-1])*1e9),final_vapor_mass_relative_error=abs(s['vapor_mass_ratio_final']/sm['vapor_mass_ratio_final']-1),peak_wall_temperature_difference_K=abs(s['max_Twall_K']-sm['max_Twall_K']),energy_residual_relative=s['energy_residual_relative']))
    for name,kw in [('phase_off',dict(phase_on=False)),('accommodation_low',dict(accommodation=1e-4)),('accommodation_high',dict(accommodation=1e-2)),('absorption_low',dict(absorption_fraction=1e-5)),('absorption_high',dict(absorption_fraction=1e-3))]:
        o,s=run_phase(PhaseParameters(**kw));s['test']=name;phases.append(s);np.savez_compressed(OUT/f'{name}.npz',**o);print(name,s['vapor_mass_ratio_final'],flush=True)
    dump('phase_controls',phases);dump('phase_convergence',conv)
    assert max(x['energy_residual_relative'] for x in rows)<1e-6
    assert max(x['energy_residual_relative'] for x in phases)<1e-6
    assert all(x['terminated'] is None for x in rows)
    assert maxerr<1e-12
    (OUT/'run_manifest.json').write_text(json.dumps(dict(seconds=time.time()-started,python=platform.python_version(),ensemble_cases=len(rows),phase_cases=len(phases),checks=len(checks),experimental_flash_comparison='Outside study scope: prescribed optical proxies, with no calibrated flash trace, response function or spectrum in supplied records.'),indent=2))
    print('COMPLETE',time.time()-started,flush=True)

if __name__=='__main__':main()
