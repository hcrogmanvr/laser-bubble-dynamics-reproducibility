FIGURE 5 - Original multibubble 2T plasma EOS simulation

FILES
multibubble_2T_plasma_original.py - original user-provided code, unchanged.
export_figure5.py - rerun original simulation and export figure, traces and diagnostics.
Figure_5_Multibubble_2T.png - high-resolution (300 dpi) six-panel figure.
Figure_5_Multibubble_2T.pdf - vector PDF figure.
Figure_5_model_output.npz - numerical simulation output.
Figure_5_numerical_traces.csv - stored/downsampled timeseries.
RUN_DIAGNOSTICS.txt - reproducibility information and modeling caveats.

REPRODUCE
Install numpy and matplotlib; from this directory run:
  python export_figure5.py

The original file is unchanged; the export script imports it and uses its plotting method.
The original plot contains a title, as provided in source. No equation or parameter was modified.
The model uses proxy EOS/ionization and event-triggered optical emission.
The calculated light-curve is not independent experimental validation.
