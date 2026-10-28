from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];RES=ROOT/'results';FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def finish(fig,name):
    for i,a in enumerate(fig.axes):
        a.grid(alpha=.2);a.text(.01,.97,f'({chr(97+i)})',transform=a.transAxes,va='top',weight='bold')
    fig.savefig(FIG/(name+'.png'),dpi=200);fig.savefig(FIG/(name+'.pdf'));plt.close(fig)

fig,ax=plt.subplots(4,2,figsize=(10.5,11.5),layout='constrained')
for j,v in enumerate([5,6]):
    x=np.load(RES/f'corrected_figure{v}.npz');z=np.load(RES/f'figure{v}_one_temperature.npz');t=x['t_us']
    ax[0,j].semilogx(np.maximum(t,1e-5),x['R'].mean(axis=1)*1e6,label=f'Corrected source {v}')
    ax[0,j].set(ylabel='Mean radius (µm)',xlabel='Time (µs)');ax[0,j].legend(loc='upper right',fontsize=9)
    ax[1,j].plot(t,x['Te'].mean(axis=1),label='Electron-reservoir mean')
    ax[1,j].plot(t,x['Ti'].mean(axis=1),label='Heavy-reservoir mean')
    ax[1,j].plot(z['t_us'],z['Te'].mean(axis=1),'--',label='Matched one-temperature mean')
    ax[1,j].set(ylabel='Effective temperature (K)',xlabel='Time (µs)');ax[1,j].legend(fontsize=8)
    ev=x['events'];ax[2,j].step(np.r_[0,ev[:,0],t[-1]],np.r_[0,np.arange(1,len(ev)+1),len(ev)],where='post')
    ax[2,j].set(xlabel='Time (µs)',ylabel='Cumulative threshold events',xlim=(0,.02))
    ax[2,j].text(.98,.25,'Initial threshold crossings\nare not measured flashes',transform=ax[2,j].transAxes,ha='right',fontsize=9)
    ax[3,j].plot(t,x['energy_residual_J']/max(x['states'][-1,-5]*1e-9,1e-18))
    ax[3,j].set(xlabel='Time (µs)',ylabel='Normalized energy residual')
finish(fig,'Figure12_corrected_ensembles')

fig,ax=plt.subplots(2,2,figsize=(10.5,7.5),layout='constrained')
for name,label in [('phase_on','Phase change ON'),('phase_off','Phase change OFF')]:
    x=np.load(RES/f'{name}.npz');t=x['t_us']
    ax[0,0].plot(t,x['R_m']*1e6,label=label)
    ax[1,0].plot(t,x['vapor_kg']/x['vapor_kg'][0],label=label)
    ax[1,1].plot(t,x['energy_residual_J']/7.5e-10,label=label)
x=np.load(RES/'phase_on.npz');ax[0,1].plot(x['t_us'],x['Twall_K'],label='Heated interface/absorber node');ax[0,1].plot(x['t_us'],x['Tgas_K'],label='Gas–vapor mixture')
for a,yl in zip(ax.flat,['Radius (µm)','Temperature (K)','Vapor mass / initial vapor mass','Normalized energy residual']):
    a.set(xlabel='Time (µs)',ylabel=yl);a.legend(fontsize=8)
ax[0,0].legend(fontsize=8,loc='lower right')
finish(fig,'Figure13_phase_change_controls')

fig,ax=plt.subplots(1,2,figsize=(10.5,3.8),layout='constrained')
for v in [2,3]:
    d=np.load(RES/f'unclipped_source{v}.npz');ax[0].plot(d['t'],(d['R'][0]-2e-6)*1e9,label=f'Source {v}')
ax[0].set(xlabel='Time (s)',ylabel='Radius displacement (nm)');ax[0].legend()
from corrected_ensemble import drude
ne=np.logspace(15,29,250);r,t,a=drude(ne,np.full_like(ne,2e-7))
for z,l in [(r,'Reflectance'),(t,'Transmittance'),(a,'Absorptance')]:ax[1].semilogx(ne,z,label=l)
ax[1].set(xlabel='Assumed electron density (m⁻³)',ylabel='Optical power fraction',ylim=(-.02,1.02));ax[1].legend(fontsize=8)
finish(fig,'FigureS1_reference_and_optical_checks')
