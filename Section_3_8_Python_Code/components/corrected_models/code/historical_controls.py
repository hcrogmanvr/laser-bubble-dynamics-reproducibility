"""Source-faithful numerical controls plus explicitly identified extensions.
Run from any directory: python run_controls.py
No claim of experimental validation is made. See README and result scope fields.
"""
from pathlib import Path
from math import pi, exp
import importlib.util,json,csv,hashlib,time,platform,sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)

def load_module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod);return mod

def cfg(model):
 common=dict(model=model,R0=2e-6,Rmax=30e-6,Pamb=1e5,c_sound=1500.,rho=1000.,cooldown=.5,eta_rad=1e-3,tau_flash=.002)
 if model==2:common.update(Rmin=.5e-6,Rc=1.5e-6,k=5e7,damping=3e4,mass=1e-12,vmax=200.,tau_E=15.,source=8.,Pcap=1.5e9,burst_fraction=.4,Ce=None)
 else:common.update(Rmin=.4e-6,Rc=1.2e-6,k=6e7,damping=4e4,mass=4/3*pi*1000*(2e-6)**3,vmax=250.,tau_E=20.,source=10.,Pcap=2e9,burst_fraction=.35,Ce=2e-9)
 return common

def legacy(model=2,dt=.01,duration=300.,forcing=1.,on=60.,off=30.,rmin_factor=1.,rmax_factor=1.,vmax_factor=1.,pcap_factor=1.,cooldown_factor=1.,growing_floor=None,tau0=1.,volume_gain=1.):
 """Only scalar-loop optimization; same source update order and sampled optics.
 Flash tails use the algebraically equivalent exponential recurrence; baseline
 comparison against untouched source verifies numerical differences explicitly.
 """
 p=cfg(model);growing_floor=(model==3) if growing_floor is None else growing_floor
 t=np.arange(0,duration,dt);N=len(t);gate=((t%(on+off))<on).astype(float)
 names=['R','v','E','P','E_burst','burst_flag','I_flash','m_np','T','A','floor','raw_v','raw_R']
 s={k:np.zeros(N) for k in names};R,v,E,P,Eb,B,I,M,T,A,F,RV,RR=[s[k] for k in names]
 R[0]=p['R0'];F[0]=p['Rmin']*rmin_factor;T[0]=exp(-tau0);A[0]=1-T[0]
 last=-np.inf;armed=True;vmax=p['vmax']*vmax_factor;rmax=p['Rmax']*rmax_factor
 rmin=p['Rmin']*rmin_factor;cap=p['Pcap']*pcap_factor
 re=exp(-dt/p['tau_flash']);vc=rc=pc=ec=0;input_energy=removed=0.
 for i in range(1,N):
  oldE=E[i-1];oldR=R[i-1];oldV=v[i-1];oldM=M[i-1]
  floor=rmin*(1+volume_gain*oldM)**(1/3) if growing_floor else rmin;F[i]=floor
  drive=gate[i] if model==2 else min((1+4*oldM)*gate[i],3.)
  src=p['source']*forcing*drive;input_energy+=src*dt
  en=max(oldE+(src-oldE/p['tau_E'])*dt,0.)
  rawp=2e8*en if model==2 else .25*en/(4/3*pi*max(oldR,floor)**3)
  pres=min(rawp,cap);pc+=rawp>cap;P[i]=pres
  force=-p['k']*(oldR-p['R0'])+4*pi*oldR**2*(pres-p['Pamb'])-p['damping']*oldV
  rv=oldV+force/p['mass']*dt;RV[i]=rv;vc+=abs(rv)>vmax
  vv=min(max(rv,-vmax),vmax);v[i]=vv
  rr=oldR+vv*dt;RR[i]=rr;rc+=(rr<floor or rr>rmax)
  R[i]=min(max(rr,floor),rmax)
  if R[i]<p['Rc'] and vv<0 and armed and t[i]-last>p['cooldown']*cooldown_factor:
   B[i]=1.;Eb[i]=p['burst_fraction']*en;en-=Eb[i];removed+=Eb[i];last=t[i];armed=False
  if R[i]>1.2*p['Rc']:armed=True
  E[i]=en
  I[i]=I[i-1]*re+p['eta_rad']*Eb[i]/p['tau_flash']
  if model==3:
   M[i]=max(oldM+(.005*(1-gate[i])-.002*I[i])*dt,0.)
   depth=tau0*(1+4*oldM);T[i]=exp(-depth);A[i]=-np.expm1(-depth)
 s.update(t=t,laser=gate)
 if model==2:
  w=max(1,int(.005/dt));s['I_cam']=np.convolve(I,np.ones(w)/w,mode='same')
 else:s['I_cam']=I
 events=t[B>0];offmins=[]
 for start in np.arange(on,duration,on+off):
  end=start+off
  if end>duration+1e-10:continue
  mask=(t>=start)&(t<end)
  if np.any(mask):offmins.append(float(R[mask].min()*1e6))
 met=dict(model=f'Figure {model} discrete source',dt_s=dt,duration_s=duration,forcing_factor=forcing,on_s=on,off_s=off,
  rmin_factor=rmin_factor,rmax_factor=rmax_factor,vmax_factor=vmax_factor,pcap_factor=pcap_factor,cooldown_factor=cooldown_factor,
  growing_floor=growing_floor,tau0=tau0,volume_gain=volume_gain,events=len(events),first_event_s=float(events[0]) if len(events) else None,
  R_min_um=float(R.min()*1e6),R_max_um=float(R.max()*1e6),v_peak_m_s=float(abs(v).max()),E_max_proxy=float(E.max()),E_final_proxy=float(E[-1]),
  pressure_peak_MPa=float(P.max()*1e-6),energy_input_proxy=input_energy,event_energy_proxy=removed,radiated_event_energy_proxy=p['eta_rad']*removed,
  sampled_flash_integral_proxy=float(np.trapezoid(I,t)),velocity_clip_fraction=vc/(N-1),radius_clip_fraction=rc/(N-1),pressure_cap_fraction=pc/(N-1),
  off_minima_um=offmins,final_floor_um=float(F[-1]*1e6))
 return s,met

def grid(duration):
 regular=np.linspace(0,duration,int(duration/.02)+1)
 # Sampling is independent of integration: capture startup and gate transitions.
 tiny=np.r_[0.,np.logspace(-18,0,240)]
 extras=[tiny]
 for switch in [60.,90.,150.,180.,240.,270.]:
  if switch<=duration:extras.extend([switch+np.r_[0.,np.logspace(-9,-1,80)],switch-np.logspace(-9,-1,80)])
 return np.unique(np.r_[regular,*extras])[(np.unique(np.r_[regular,*extras])>=0)&(np.unique(np.r_[regular,*extras])<=duration)]

def continuous(model=2,duration=300.,forcing=1.,on=60.,off=30.,rtol=1e-7,max_step=1.,pcap_factor=1.,n_bubbles=1,spacing_um=None,coupling=False):
 """Integrate the source spring/damping force balance using implicit Radau.
 NO radius/velocity clipping. Same pressure ceiling. Energy is integrated exactly
 on the no-event branch within each gate interval. If an event is detected, this
 reference STOPS rather than silently omitting the source jump. The runs reported
 here never reach the source trigger radius.
 For N>1 only: NEW quasi-instantaneous far-field pressure coupling, not an
 implementation recovered from Figures 4-6. See dynamic_coupling_equations.md.
 """
 p=cfg(model);rs,vs=1e-9,1e-6;n=n_bubbles
 if n>1 and (spacing_um is None or spacing_um<=0):raise ValueError('spacing required')
 pos=np.arange(n)*(spacing_um or 1.)*1e-6
 dist=np.abs(pos[:,None]-pos[None,:]);inv=np.zeros((n,n));np.divide(1.,dist,out=inv,where=dist>0)
 tt=grid(duration);ys=np.zeros((2*n,len(tt)));ener=np.zeros(len(tt));gates=np.zeros(len(tt));loads=np.zeros(len(tt))
 y0=np.zeros(2*n);e0=0.;m0=0.;start=0.;success=True;messages=[];nfev=0;trigger_times=[];minsep=float('inf')
 bounds=sorted(set([0.,duration]+[x for base in np.arange(0,duration,on+off) for x in [float(base),float(base+on)] if 0<x<duration]))
 for a,b in zip(bounds[:-1],bounds[1:]):
  gate=1. if ((a+min(1e-8,(b-a)/2))%(on+off))<on else 0.
  drive=gate*(1. if model==2 else min(1+4*m0,3.));source=p['source']*forcing*drive
  def energy(t):return e0*exp(-t/p['tau_E'])+source*p['tau_E']*(-np.expm1(-t/p['tau_E']))
  def pressures(R,t):
   en=energy(t)
   raw=2e8*en*np.ones(n) if model==2 else .25*en/(4/3*pi*R**3)
   return np.minimum(raw,p['Pcap']*pcap_factor)
  def fun(t,y):
   R=p['R0']+rs*y[:n];v=vs*y[n:];P=pressures(R,t)
   force=-p['k']*(rs*y[:n])+4*pi*R**2*(P-p['Pamb'])-p['damping']*v
   if coupling and n>1:
    area=4*pi*R**2
    mass=np.eye(n)*p['mass']+area[:,None]*p['rho']*inv*(R[None,:]**2)
    rhs=force-area*p['rho']*(inv@(2*R*v*v))
    acceleration=np.linalg.solve(mass,rhs)
   else:acceleration=force/p['mass']
   return np.r_[v/rs,acceleration/vs]
  def jac(t,y):
   R=p['R0']+rs*y[:n];P=pressures(R,t)
   raw=.25*energy(t)/(4/3*pi*R**3) if model==3 else np.zeros(n)
   dp=np.where(raw<p['Pcap']*pcap_factor,-3*raw/R,0.) if model==3 else np.zeros(n)
   J=np.zeros((2*n,2*n));J[:n,n:]=np.eye(n)*vs/rs
   J[n:,:n]=np.diag((-p['k']+8*pi*R*(P-p['Pamb'])+4*pi*R**2*dp)*rs/(p['mass']*vs));J[n:,n:]=-np.eye(n)*p['damping']/p['mass'];return J
  events=[]
  for j in range(n):
   def event(t,y,j=j):return p['R0']+rs*y[j]-p['Rc']
   event.direction=-1;event.terminal=True;events.append(event)
  opts=dict(method='Radau',rtol=rtol,atol=np.r_[np.full(n,rtol*.01),np.full(n,max(1e-6,rtol*10))],max_step=max_step,dense_output=True,events=events)
  if n==1:opts['jac']=jac
  sol=solve_ivp(fun,(0.,b-a),y0,**opts);nfev+=sol.nfev
  messages.append(sol.message)
  trigger_times.extend([float(ti+a) for x in sol.t_events for ti in x])
  if not sol.success or trigger_times:
   success=False;print('FAIL',a,b,sol.t[-1],sol.message,flush=True);break
  mask=(tt>=a)&(tt<=b);ys[:,mask]=sol.sol(tt[mask]-a);ener[mask]=[energy(ti-a) for ti in tt[mask]];gates[mask]=gate;loads[mask]=m0+(0.005*(tt[mask]-a) if model==3 and gate==0 else 0.)
  y0=sol.y[:,-1];e0=energy(b-a)
  if model==3 and gate==0:m0+=.005*(b-a)
 R=p['R0']+rs*ys[:n];v=vs*ys[n:]
 P=np.minimum(2e8*ener[None,:]*np.ones((n,1)) if model==2 else .25*ener[None,:]/(4/3*pi*R**3),p['Pcap']*pcap_factor)
 met=dict(model=f'Figure {model} continuous spring-damping equations'+(' + NEW pressure-coupling extension' if n>1 else ''),
  solver='Radau',rtol=rtol,radius_atol_nm=rtol*.01,velocity_atol_m_s=max(1e-6,rtol*10)*vs,max_step_s=max_step,duration_s=duration,forcing_factor=forcing,on_s=on,off_s=off,pcap_factor=pcap_factor,
  n_bubbles=n,spacing_um=spacing_um,coupling=coupling,success=success,nfev=nfev,events=len(trigger_times),
  R_min_um=float(R.min()*1e6),R_max_um=float(R.max()*1e6),v_peak_m_s=float(abs(v).max()),E_max_proxy=float(ener.max()),E_final_proxy=float(ener[-1]),
  pressure_peak_MPa=float(P.max()*1e-6),radius_clip_fraction=0.,velocity_clip_fraction=0.,R_max_excursion_nm=float((R-p['R0']).max()*1e9),
  R_min_excursion_nm=float((R-p['R0']).min()*1e9),minimum_radius_margin_to_trigger_um=float((R.min()-p['Rc'])*1e6),
  mechanical_relaxation_fast_s=p['mass']/p['damping'],mechanical_relaxation_slow_s=p['damping']/p['k'],
  eps_inertia=p['mass']*p['k']/p['damping']**2)
 return dict(t=tt,R=R,v=v,E=ener,P=P,gate=gates,m_np=loads),met

def dump_rows(name,rows):
 (OUT/(name+'.json')).write_text(json.dumps(rows,indent=2,allow_nan=False))
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with (OUT/(name+'.csv')).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader()
  for r in rows:wr.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()})

def main():
 start=time.time();manifest=[]
 for name in ['figure2_updated.py','figure3_updated.py','reconstructed_model.py']:
  path=ROOT/'source_originals'/name;manifest.append(dict(file=name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 (OUT/'provenance.json').write_text(json.dumps(dict(original_sources=manifest,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
  scopes={'legacy':'Parameter changes to the supplied discrete source equations; no physical convergence inferred.',
  'continuous':'Unclipped spring-damping force equations from supplied Figure 2/3 sources, not Rayleigh-Plesset dynamics.',
  'dynamic_coupling':'New explicitly documented far-field extension of source mechanical equations, not recovery of archived multi-bubble manuscript figures.',
  'reconstructed':'Archived reconstruction; README explicitly states original source was not available when created.'}),indent=2))
 # Verify optimized loop against the untouched uploaded source implementation.
 repro=[]
 for model in [2,3]:
  mod=load_module('original'+str(model),ROOT/'source_originals'/f'figure{model}_updated.py')
  src=mod.simulate() if model==2 else mod.simulate(True)
  s,m=legacy(model);np.savez_compressed(OUT/f'legacy_figure{model}_baseline.npz',**s)
  keys=['R','v','E','P','burst_flag','E_burst']+(['I_cam'] if model==2 else ['I_flash','m_np','T','A','R_floor'])
  errors={}
  for key in keys:
   target=s['floor'] if key=='R_floor' else s[key]
   err=float(np.max(np.abs(target-src[key])));den=float(np.max(np.abs(src[key])))
   errors[key]={'max_absolute_error':err,'max_relative_to_peak':err/max(den,1e-300)}
  assert all(z['max_relative_to_peak']<1e-10 for z in errors.values()),errors
  repro.append(dict(model=model,event_count=m['events'],comparison=errors));print('Source reproduced',model,m['events'],flush=True)
 dump_rows('source_reproduction',repro)
 # Legacy refinement and hard-bound sensitivity, explicitly not an integration repair.
 rows=[]
 for model in [2,3]:
  for dt in [.01,.005,.0025,.00125,.000625,.0003125]:
   _,m=legacy(model,dt=dt);m['test']='timestep';rows.append(m);print('Legacy dt',model,dt,m['events'],flush=True)
  for parameter in ['rmin_factor','rmax_factor','vmax_factor','pcap_factor','cooldown_factor']:
   for factor in [.5,1.,2.]:
    _,m=legacy(model,**{parameter:factor});m.update(test='bound_or_event_rule',parameter=parameter,factor=factor);rows.append(m)
  for factor in [.5,.75,1.,1.25,1.5]:
   _,m=legacy(model,forcing=factor);m.update(test='forcing');rows.append(m)
  for on,off in [(30.,60.),(45.,45.),(60.,30.),(75.,15.)]:
   _,m=legacy(model,on=on,off=off);m.update(test='duty_cycle_fixed_peak_source');rows.append(m)
 # Optical diagnostic null controls: these parameters do not enter original energy source.
 s0,m0=legacy(3)
 for tau in [.5,1.,2.]:
  sx,m=legacy(3,tau0=tau);m.update(test='optical_postprocessing_null',max_dynamic_state_change=max(float(np.max(abs(sx[k]-s0[k]))) for k in ['R','v','E','P','E_burst','m_np']),T_min=float(sx['T'].min()),T_max=float(sx['T'].max()));rows.append(m)
 dump_rows('legacy_controls',rows)
 # Genuine refinement of the source continuous force equations.
 refs={};rr=[]
 for model in [2,3]:
  sols=[]
  for rtol,ms in [(1e-5,5.),(1e-7,1.),(1e-9,.25)]:
   s,m=continuous(model,rtol=rtol,max_step=ms);assert m['success'],m;sols.append((s,m));rr.append(m);print('Continuous',model,rtol,m['events'],m['R_max_um'],flush=True)
  fine=sols[-1][0];refs[model]=fine
  for s,m in sols:
   m['max_radius_difference_from_finest_nm']=float(np.max(abs(s['R']-fine['R']))*1e9)
   m['max_velocity_difference_from_finest_m_s']=float(np.max(abs(s['v']-fine['v'])))
  np.savez_compressed(OUT/f'continuous_figure{model}_reference.npz',**fine)
 dump_rows('continuous_convergence',rr)
 # Controls on the stable no-event branch: forcing and physical pressure ceiling.
 rr=[]
 for model in [2,3]:
  for factor in [.5,.75,1.,1.25,1.5]:
   s,m=continuous(model,forcing=factor);assert m['success'];m['test']='forcing';rr.append(m)
  for factor in [.5,1.,2.]:
   s,m=continuous(model,pcap_factor=factor);assert m['success'];m['test']='pressure_cap';rr.append(m)
  for on,off in [(30.,60.),(45.,45.),(60.,30.),(75.,15.)]:
   s,m=continuous(model,on=on,off=off);assert m['success'];m['test']='duty_cycle_fixed_peak_source';rr.append(m)
 dump_rows('continuous_controls',rr)
 # NEW physical pressure-mediated extension: N fixed; spacing changes, not force coefficients.
 rows=[]
 for model in [2,3]:
  for spacing in [40.,100.,200.]:
   no,mno=continuous(model,duration=2.,n_bubbles=3,spacing_um=spacing,coupling=False,rtol=1e-8,max_step=.05)
   yes,my=continuous(model,duration=2.,n_bubbles=3,spacing_um=spacing,coupling=True,rtol=1e-8,max_step=.05)
   finer,mf=continuous(model,duration=2.,n_bubbles=3,spacing_um=spacing,coupling=True,rtol=1e-10,max_step=.01)
   assert mno['success'] and my['success'] and mf['success']
   my['comparison_window_s']=2.;my['max_radius_on_off_difference_nm']=float(np.max(abs(yes['R']-no['R']))*1e9)
   my['max_radius_refinement_difference_nm']=float(np.max(abs(yes['R']-finer['R']))*1e9)
   my['max_velocity_on_off_difference_m_s']=float(np.max(abs(yes['v']-no['v'])))
   my['max_radius_to_spacing_ratio']=float(yes['R'].max()/(spacing*1e-6));rows.append(my)
   np.savez_compressed(OUT/f'new_coupling_figure{model}_{spacing:g}um.npz',t=yes['t'],R_on=yes['R'],R_off=no['R'],v_on=yes['v'],v_off=no['v'])
   print('New dynamic coupling',model,spacing,my['max_radius_on_off_difference_nm'],flush=True)
 dump_rows('new_dynamic_coupling_controls',rows)
 # Archived reconstruction: actual implemented quenching and postprocessed coupling.
 mod=load_module('archived_reconstruction',ROOT/'source_originals/reconstructed_model.py')
 rows=[]
 for q in [0.,.65]:
  p=mod.Params(duration_s=1.,dt_s=1e-4,quench_strength=q)
  s=mod.simulate_single(p);g=s['gate'];T=s['temperature_K'];af=np.array([mod.absorbed_fraction(r,t,p) for r,t in zip(s['R_um'],T)])
  absorbed=p.power_W*g*af;rows.append(dict(scope='Archived reconstruction only',quench_strength=q,duration_s=1.,dt_s=p.dt_s,
   E_final_nominal_J=float(s['energy_J'][-1]),events=int(np.count_nonzero(s['burst_J'])),absorbed_energy_nominal_J=float(np.sum(absorbed[:-1])*p.dt_s),R_final_um=float(s['R_um'][-1]),
   pressure_cap_sample_fraction=float(np.mean(s['pressure_Pa']==p.pressure_cap_Pa))))
  np.savez_compressed(OUT/f'reconstructed_quench_{q:g}.npz',**s)
 dump_rows('reconstructed_quenching_control',rows)
 p=mod.Params(duration_s=.1,dt_s=1e-4);p0=mod.Params(**{**p.__dict__,'coupling_strength':0.})
 off=mod.simulate_multi(p0,n_bubbles=3,seed=7);on=mod.simulate_multi(p,n_bubbles=3,seed=7)
 physical=['R_um','velocity_m_s','energy_J','temperature_K','burst_J','camera_flash']
 audit=dict(scope='Archived reconstruction only',same_seed=7,n_bubbles=3,duration_s=.1,
  maximum_dynamic_state_differences={k:max(float(np.max(abs(x[k]-y[k]))) for x,y in zip(on,off)) for k in physical},
  maximum_postprocessed_pressure_difference_Pa=max(float(np.max(abs(x['pressure_coupled_Pa']-y['pressure_coupled_Pa']))) for x,y in zip(on,off)),
  spacing_parameter_present=False,coupling_enters_time_integration=False)
 (OUT/'archived_coupling_audit.json').write_text(json.dumps(audit,indent=2))
 (OUT/'completion_scope.json').write_text(json.dumps(dict(
  completed=['Source Figure 2/3 reproduction','Discrete-source timestep sweeps','Numerical-bound and cooldown sensitivity','Forcing and duty-cycle sweeps','Unclipped stiff integration of source force equations and refinement','Source Figure 3 optical-diagnostic null control','Quenching on/off in explicitly archived reconstruction','Archived coupling on/off and missing-spacing audit','New dynamic far-field pressure-coupling extension with spacing and refinement controls'],
  not_claimed=['Independent experimental flash validation','A converged repeat-collapse plasma regime','Quenching implemented in original Figure 2/3 energy equations','Recovery or validation of original manuscript Figures 4-6 dynamic coupling','A validated two-temperature EOS','Calibration of nominal energy coefficients to joules'],runtime_s=time.time()-start),indent=2))
 print('COMPLETE seconds',time.time()-start,flush=True)

if __name__=='__main__':main()
