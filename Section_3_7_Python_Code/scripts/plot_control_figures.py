from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; R=ROOT/'components/historical_single/results'; F=ROOT/'outputs/figures';F.mkdir(parents=True,exist_ok=True)
load=lambda n:json.loads((R/(n+'.json')).read_text())
legacy=load('legacy_controls');conv=load('continuous_convergence');couple=load('new_dynamic_coupling_controls');quench=load('reconstructed_quenching_control')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
colors=['#225b8d','#b44d31']
fig,ax=plt.subplots(1,3,figsize=(10.5,3.2))
for model,col in zip([2,3],colors):
 rr=[x for x in legacy if x['model'].startswith(f'Figure {model}') and x['test']=='timestep'];dt=[x['dt_s'] for x in rr]
 ax[0].plot(dt,[x['events'] for x in rr],'o-',color=col,label=f'Figure {model}')
 ax[1].plot(dt,[x['sampled_flash_integral_proxy']/x['radiated_event_energy_proxy'] for x in rr],'o-',color=col)
 rr=[x for x in legacy if x['model'].startswith(f'Figure {model}') and x['test']=='bound_or_event_rule' and x['parameter']=='cooldown_factor']
 ax[2].plot([x['factor']*.5 for x in rr],[x['events'] for x in rr],'o-',color=col)
for a in ax[:2]:a.set_xscale('log');a.set_xlabel('Timestep (s)');a.invert_xaxis()
ax[0].set_ylabel('Event count in 300 s');ax[0].legend(frameon=False)
ax[1].set_ylabel('Sampled flash integral / event radiation');ax[1].axhline(1,color='gray',ls='--',lw=1)
ax[2].set(xlabel='Imposed waiting time (s)',ylabel='Event count in 300 s')
for i,a in enumerate(ax):a.text(.01,1.04,chr(97+i),transform=a.transAxes,fontweight='bold')
fig.tight_layout();[fig.savefig(F/f'Figure9.{ext}',dpi=220) for ext in ('png','pdf','svg')];plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(10.5,3.3))
for model,col in zip([2,3],colors):
 s=np.load(R/f'continuous_figure{model}_reference.npz');d=(s['R'][0]-2e-6)*1e9
 ax[0].plot(s['t'],d,color=col,label=f'Figure {model}')
 mask=(s['t']>0)&(s['t']<.1);ax[1].semilogx(s['t'][mask],d[mask],color=col)
ax[0].set(xlabel='Time (s)',ylabel='Radius displacement from 2 µm (nm)');ax[0].legend(frameon=False)
ax[1].set(xlabel='Time after start (s)',ylabel='Radius displacement from 2 µm (nm)')
for i,a in enumerate(ax):a.text(.01,1.04,chr(97+i),transform=a.transAxes,fontweight='bold')
fig.tight_layout();[fig.savefig(F/f'Figure10.{ext}',dpi=220) for ext in ('png','pdf','svg')];plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(10.5,3.3))
for q,label,c in [(0,'Quenching off',colors[0]),(.65,'Quenching on',colors[1])]:
 s=np.load(R/f'reconstructed_quench_{q:g}.npz');ax[0].plot(s['t'],s['energy_J'],label=label,color=c)
ax[0].set(xlabel='Time (s)',ylabel='Stored energy (nominal J)');ax[0].legend(frameon=False)
for model,col in zip([2,3],colors):
 rr=[x for x in couple if x['model'].startswith(f'Figure {model}')]
 ax[1].semilogy([x['spacing_um'] for x in rr],[x['max_radius_on_off_difference_nm'] for x in rr],'o-',color=col,label=f'Figure {model} on/off')
 ax[1].semilogy([x['spacing_um'] for x in rr],[x['max_radius_refinement_difference_nm'] for x in rr],'x--',color=col,label=f'Figure {model} refinement')
ax[1].set(xlabel='Center spacing (µm)',ylabel='Maximum radius difference (nm)');ax[1].legend(frameon=False,fontsize=8)
for i,a in enumerate(ax):a.text(.01,1.04,chr(97+i),transform=a.transAxes,fontweight='bold')
fig.tight_layout();[fig.savefig(F/f'Figure11.{ext}',dpi=220) for ext in ('png','pdf','svg')];plt.close(fig)
