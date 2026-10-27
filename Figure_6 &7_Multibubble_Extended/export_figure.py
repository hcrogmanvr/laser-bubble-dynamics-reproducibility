"""Run the UNMODIFIED supplied Figure 6 model and export its two figures/data."""
from pathlib import Path
import importlib.util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv

folder = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('original_figure6', folder/'multibubble_full_extended.py')
module = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = module
spec.loader.exec_module(module)

out = module.simulate(N=8, t_end=1e-3, dt=1e-7, stride=5, seed=2, cluster_sigma=2e-6)

# Keep exact model plots; capture their figures when original plot_all() calls plt.show().
plt.show = lambda *args, **kwargs: None
module.plot_all(out)
figs = [plt.figure(num) for num in plt.get_fignums()]
assert len(figs) == 2, f'Expected 2 figures; got {len(figs)}'
for fig, name in zip(figs, ['Figure_6_Multibubble_Extended', 'Figure_6_Membrane_Rupture']):
    fig.savefig(folder/f'{name}.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(folder/f'{name}.pdf', bbox_inches='tight', facecolor='white')
plt.close('all')

arrays = {k:v for k,v in out.items() if isinstance(v,np.ndarray)}
np.savez_compressed(folder/'Figure_6_model_output.npz', **arrays)
keys = ['t','R_mean','E_mean','P_over_MPa','Te_mean','I_cam','laser','burst_flag','xion_mean','M_mean','rupt_frac']
with (folder/'Figure_6_numerical_traces.csv').open('w', newline='') as f:
    w=csv.writer(f);w.writerow(keys)
    w.writerows(zip(*[out[k] for k in keys]))

lines = [
 'Figure 6: Original supplied model, unmodified.',
 'Simulation duration: 1.0 ms',
 'Integration dt: 1.0e-7 s',
 'Number of bubbles: 8',
 f"Total model-triggered bursts: {out['total_bursts']}",
 f"Saved time samples: {len(out['t'])}",
 f"Final ruptured fraction: {np.mean(out['final_ruptured']):.6g}",
 f"Peak mean radius (um): {np.max(out['R_mean'])*1e6:.6g}",
 f"Peak mean Te (K): {np.max(out['Te_mean']):.6g}",
 f"Peak modeled overpressure (MPa): {np.max(out['P_over_MPa']):.6g}",
 f"Peak camera response (arb): {np.max(out['I_cam']):.6g}",
 f"Max modeled mass/accretion state: {np.max(out['M_mean']):.6g}",
 f"Max mean ionization fraction: {np.max(out['xion_mean']):.6g}",
 f"Maximum accumulated pressure impulse (Pa s): {np.max(out['final_Jmem']):.6g}",
 'Caveats: The model uses bounded acceleration, velocity and radius; parameterized EOS, rupture and emission thresholds; camera saturates at 2 arbitrary units.',
 'Event flags are stored only at stride=5, so the saved burst_flag cannot recover all per-step events.',
 'No changes were made to the source model equations or parameters.',
]
(folder/'RUN_DIAGNOSTICS.txt').write_text('\n'.join(lines)+'\n')
(folder/'README.txt').write_text('Figure 6 reproducibility package\n\n1. Install numpy and matplotlib.\n2. Run python export_figure.py to reproduce and save both plots, numerical results and diagnostics.\n3. Or run python multibubble_full_extended.py to display the original two plots.\n\nThe provided original code is copied without changes. The exporter calls the original simulate and plot_all functions with the original main() arguments.\nThe model-generated bursts and ionization are proxies, not direct optical measurements or independently validated plasma outcomes.\n')
print('\n'.join(lines))
