Figure 4: original user-supplied multi-bubble model
===================================================
Files:
- multibubble_fixed_full.py : original supplied source script, including plt.show()
- export_figure.py : runner that imports the original simulation and saves its same six-panel visualization
- Figure_4_Multibubble.png / .pdf : plotted result
- Figure_4_model_output.npz : unmodified simulation output arrays
- RUN_DIAGNOSTICS.txt : outcomes from the original default model

Run: python export_figure.py
Dependencies: numpy, matplotlib

Interpretation cautions:
- The 8 counted bursts occur on just ONE numerical time step; they are simultaneous early crossings, not eight separately timed camera flashes.
- The mean radius settles at the imposed lower bound R_min=0.05 micrometers, so the plotted radius trace is bound-limited.
- E and electron temperature are proxies, not a calibrated plasma equation of state.
- This is a reproduction of the supplied code, NOT a validation against experiment.
