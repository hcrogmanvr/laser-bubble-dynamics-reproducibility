from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;RES=ROOT/'results';FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def finish(fig,name):
    for i,a in enumerate(fig.axes):
        a.grid(alpha=.2);a.text(.01,.97,f'({chr(97+i)})',transform=a.transAxes,va='top',weight='bold')
    fig.savefig(FIG/(name+'.png'),dpi=200);fig.savefig(FIG/(name+'.pdf'));plt.close(fig)

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

