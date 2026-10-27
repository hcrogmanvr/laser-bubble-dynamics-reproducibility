"""Run the unmodified Figure 5 source and export its original six-panel plot.

Run: python export_figure5.py
Requires numpy, matplotlib.
"""
from pathlib import Path
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
source = ROOT / 'multibubble_2T_plasma_original.py'
spec = importlib.util.spec_from_file_location('figure5_original', source)
model = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = model
spec.loader.exec_module(model)

out = model.run_multibubble_2T_plasma(N=8,t_end=1.0e-3,dt=1.0e-7,stride=5,seed=2)
original_show = plt.show

def save_instead_of_show(*args, **kwargs):
    fig = plt.gcf()
    fig.savefig(ROOT / 'Figure_5_Multibubble_2T.png', dpi=300, bbox_inches='tight')
    fig.savefig(ROOT / 'Figure_5_Multibubble_2T.pdf', bbox_inches='tight')

plt.show = save_instead_of_show
try:
    model.plot_results(out, title='Multi-bubble (2T plasma EOS): radius, energy, pressure, Te, camera flash')
finally:
    plt.show = original_show
    plt.close('all')

np.savez_compressed(ROOT / 'Figure_5_model_output.npz', **out)
import csv
fields=['t','R_mean','E_mean','P_mean_MPa','Te_mean','I_cam','laser','burst_flag']
with (ROOT/'Figure_5_numerical_traces.csv').open('w', newline='') as f:
    w=csv.writer(f)
    w.writerow(['time_s','R_mean_m','E_mean_J','P_mean_MPa','Te_mean_K','I_cam_arb','laser_0_1','burst_flag_downsampled'])
    w.writerows(zip(*(out[field] for field in fields)))

lines=[
    'Figure 5: original supplied code reproduced without parameter or equation changes',
    'Inputs: N=8, t_end=0.001 s, dt=1e-7 s, stride=5, seed=2',
    f"Total burst events (all simulation steps): {out['total_bursts']}",
    f"Stored points: {len(out['t'])}",
    f"Recorded/downsampled burst flags: {int(np.count_nonzero(out['burst_flag']))}",
    f"Mean bubble radius range (um): {1e6*np.min(out['R_mean']):.9g} to {1e6*np.max(out['R_mean']):.9g}",
    f"Mean electron temperature range (K): {np.min(out['Te_mean']):.9g} to {np.max(out['Te_mean']):.9g}",
    f"Mean pressure range (MPa): {np.min(out['P_mean_MPa']):.9g} to {np.max(out['P_mean_MPa']):.9g}",
    f"Maximum camera signal (arb): {np.max(out['I_cam']):.9g}",
    f"Fraction of stored camera points at sensor saturation=2.0: {np.mean(np.isclose(out['I_cam'], 2.0)):.4%}",
    '',
    'CAUTIONS:',
    '- Ionization fraction is a logistic proxy, not a calibrated Saha/CR model.',
    '- Pressure is capped; radius and radial speed have numerical clamps.',
    '- Burst triggering and conversion to camera signal are phenomenological.',
    '- Burst flags are downsampled; total burst count uses every simulation step.',
    '- The presented figure reproduces the supplied code; this does not establish physical validation.',
]
(ROOT/'RUN_DIAGNOSTICS.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:9]))
