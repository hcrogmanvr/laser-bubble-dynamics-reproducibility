"""Paired pressure-coupling controls for the original, dimensionally defective sources.
This audits algorithm sensitivity; it does not repair or validate bubble physics.
"""
from pathlib import Path
import importlib.util,sys,json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent
PKG=ROOT/'package' if (ROOT/'package').exists() else ROOT.parent
DEST=ROOT/'coupling_results' if (ROOT/'package').exists() else ROOT/'coupling_results'
DEST.mkdir(exist_ok=True)
DT=1e-7
names={5:'figure5_two_reservoir_original.py',6:'figure6_extended_original.py'}
summary=[];cases={}
for fig,filename in names.items():
 for enabled in [True,False]:
  tag=f'figure{fig}_coupling_'+('on' if enabled else 'off')
  src=PKG/'source_originals'/filename
  spec=importlib.util.spec_from_file_location(tag,src);mod=importlib.util.module_from_spec(spec);sys.modules[tag]=mod;spec.loader.exec_module(mod)
  cls=mod.CouplingModel if fig==5 else mod.Coupling
  if not enabled:
   cls.coupled_pressure=lambda self,W,P:np.zeros_like(P)
   if fig==6:mod.Shock.amplitude=lambda self,Eburst,r:0.0
  original_step=mod.BubbleDynamics.step
  masks=[];hist_R=[]
  def observe(self,R,v,*args,_original=original_step):
   newR,newv=_original(self,R,v,*args)
   masks.append((R>1.5e-7)&(newR<=1.5e-7)&(v<0))
   hist_R.append(np.array(newR,copy=True))
   return newR,newv
  mod.BubbleDynamics.step=observe
  out=mod.run_multibubble_2T_plasma() if fig==5 else mod.simulate()
  mask=np.asarray(masks);assert int(mask.sum())==int(out['total_bursts'])
  full_t=np.arange(len(mask))*DT
  ev=[full_t[mask[:,i]] for i in range(mask.shape[1])]
  if enabled:
   old=np.load(PKG/f'returned_arrays/figure{fig}_returned_arrays.npz')
   equal={k:bool(np.array_equal(old[k],np.asarray(out[k]))) for k in old.files if k!='runtime_s' and k in out}
   assert equal and all(equal.values()),equal
  else:equal={}
  arrays={k:np.asarray(v) for k,v in out.items() if k!='runtime_s'}
  arrays.update(full_grid_labels_s=full_t,per_bubble_event_mask=mask,per_bubble_radius_m=np.asarray(hist_R))
  np.savez_compressed(DEST/f'{tag}.npz',**arrays)
  pkey='P_mean_MPa' if 'P_mean_MPa' in out else 'P_over_MPa'
  row={'figure':fig,'coupling':'on' if enabled else 'off','events':int(out['total_bursts']),'event_counts_per_bubble':[len(x) for x in ev],'peak_stored_mean_Te_K':float(np.max(out['Te_mean'])),'peak_stored_mean_pressure_MPa':float(np.max(out[pkey])),'mean_radius_peak_um':float(np.max(out['R_mean'])*1e6),'final_flag_fraction':float(out['rupt_frac'][-1]) if fig==6 else None,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'default_returned_arrays_equal_to_archive':equal}
  summary.append(row);cases[(fig,enabled)]=(ev,row)

# Event phases interpolate between consecutive individual triggers. Compare on/off on
# exactly the same grid and the intersection of all valid bracketing intervals.
for fig in names:
 all_lists=[e for enabled in [True,False] for e in cases[(fig,enabled)][0]]
 valid=all(len(e)>=2 for e in all_lists)
 if valid:
  lo=max(e[0] for e in all_lists);hi=min(e[-1] for e in all_lists)
  grid=np.arange(np.ceil(lo/DT),np.floor(hi/DT)+1,dtype=np.int64)*DT
  valid=hi>lo and len(grid)>1
 if valid:
  for enabled in [True,False]:
   ev,row=cases[(fig,enabled)];phases=[]
   for events in ev:
    idx=np.searchsorted(events,grid,side='right')-1;idx=np.clip(idx,0,len(events)-2)
    phases.append(2*np.pi*(grid-events[idx])/(events[idx+1]-events[idx]))
   order=np.abs(np.mean(np.exp(1j*np.array(phases)),axis=0))
   row.update(event_phase_coherence_mean=float(order.mean()),phase_window_s=[float(grid[0]),float(grid[-1])],phase_grid_points=len(grid))
   np.savetxt(DEST/(f'figure{fig}_phase_'+('on' if enabled else 'off')+'.csv'),np.column_stack([grid,order]),delimiter=',',header='grid_label_s,event_phase_coherence',comments='')
 else:
  for enabled in [True,False]:cases[(fig,enabled)][1].update(event_phase_coherence_mean=None,phase_window_s=None,phase_grid_points=0)
 # A separate, explicitly restricted statistic for bubbles that trigger at least
 # twice in BOTH conditions. Never substitute this for the full-ensemble metric.
 active=[i for i in range(8) if all(len(cases[(fig,b)][0][i])>=2 for b in [True,False])]
 if len(active)>=2:
  selected=[cases[(fig,b)][0][i] for b in [True,False] for i in active]
  lo=max(x[0] for x in selected);hi=min(x[-1] for x in selected)
  grid=np.arange(np.ceil(lo/DT),np.floor(hi/DT)+1,dtype=np.int64)*DT
  if hi>lo and len(grid)>1:
   for enabled in [True,False]:
    ev,row=cases[(fig,enabled)];ph=[]
    for i in active:
     events=ev[i];idx=np.searchsorted(events,grid,side='right')-1;idx=np.clip(idx,0,len(events)-2)
     ph.append(2*np.pi*(grid-events[idx])/(events[idx+1]-events[idx]))
    order=np.abs(np.mean(np.exp(1j*np.array(ph)),axis=0))
    row.update(active_bubble_indices=active,active_phase_coherence_mean=float(order.mean()),active_phase_window_s=[float(grid[0]),float(grid[-1])],active_phase_grid_points=len(grid))
    np.savetxt(DEST/(f'figure{fig}_active_phase_'+('on' if enabled else 'off')+'.csv'),np.column_stack([grid,order]),delimiter=',',header='grid_label_s,active_subset_event_phase_coherence',comments='')
report={'scope':'Conditional implementation audit at unchanged default timestep/seed/forcing. Original radial dimensional and event-energy defects are retained. No convergence or physical validation.','off_switch':'Figure 5: zero instantaneous coupling correction. Figure 6: zero instantaneous coupling correction and zero delayed-pulse amplitudes. All other equations/defaults unchanged.','event_observer':'Passive observation of original radial step; exact original crossing rule (old R > threshold, new R <= threshold, old v < 0), threshold 0.15 micrometer. Total observed triggers must equal original returned total_bursts.','phase_definition':'Between consecutive per-bubble events, phi_i=2*pi*(t-t_i_prev)/(t_i_next-t_i_prev). R_phi=abs(mean(exp(i*phi_i))). On/off use the same 0.1 microsecond grid over the intersection of valid bracketing intervals of every bubble in both cases; no-event/insufficient-event cases are undefined, not zero. Shared gating can also create phase coherence.','cases':summary}
(DEST/'matched_coupling_controls.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(summary,indent=2))
