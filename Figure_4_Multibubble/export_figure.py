"""Run the unmodified model and export its original six-panel figure."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multibubble_fixed_full import Medium, Camera, simulate_multibubble

folder = os.path.dirname(os.path.abspath(__file__))
out, bursts = simulate_multibubble(medium=Medium(), cam=Camera(tau=8e-6, gain=2e-3, saturate=1.0))
t = out['t']; radius = out['R']*1e6; energy=out['E']; pressure=out['P']/1e6
T=out['Te']; laser=out['laser']; cam_sig=out['cam']; events=np.where(out['burst'].sum(axis=1)>0)[0]
fig=plt.figure(figsize=(12, 10)); gs=fig.add_gridspec(6,1,hspace=.25)
axes=[fig.add_subplot(gs[i]) for i in range(6)]
for ax,y,label in zip(axes,[radius.mean(axis=1),laser,energy.mean(axis=1),pressure.mean(axis=1),T.mean(axis=1),cam_sig],['Mean R (µm)','Laser (0/1)','Mean E (arb)','Mean P (MPa)','Mean Te (K)','I_cam (arb)']):
    ax.plot(t*1e3,y); ax.set_ylabel(label)
axes[0].set_title('Multi-bubble (fixed): radius, energy, pressure, and camera flash')
if len(events): axes[-1].vlines(t[events]*1e3,0,max(cam_sig),alpha=.2,label='Burst events')
axes[-1].set_xlabel('Time (ms)'); axes[-1].legend(['Camera flash'] + (['Burst events'] if len(events) else []),loc='upper right')
fig.subplots_adjust(left=.12,right=.97,top=.94,bottom=.075,hspace=.35)
fig.savefig(os.path.join(folder,'Figure_4_Multibubble.png'),dpi=220)
fig.savefig(os.path.join(folder,'Figure_4_Multibubble.pdf'))
np.savez_compressed(os.path.join(folder,'Figure_4_model_output.npz'),**out)
with open(os.path.join(folder,'RUN_DIAGNOSTICS.txt'),'w') as f:
    f.write('Unmodified model parameters; export code changes only plotting/saving.\n')
    f.write(f'Total bursts: {bursts}\nNumber of event time steps: {len(events)}\n')
    f.write(f'Mean radius range (um): {radius.mean(axis=1).min():.9g} to {radius.mean(axis=1).max():.9g}\n')
    f.write(f'Camera max (arb): {cam_sig.max():.9g}\n')
    f.write(f'Pressure max (MPa): {pressure.max():.9g}\n')
print('Total bursts:',bursts,'event time steps:',len(events),'camera peak:',cam_sig.max())
