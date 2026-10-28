SECTION 3.9 — SEPARATE PHASE-CHANGE CALCULATIONS

This is the phase-change component formerly distributed inside the Section 3.8 package, now labeled for the reorganized manuscript. It supports Figure 13 and Table 7. The scientific solver phase_change.py is unchanged.

Files:
  phase_change.py: governing equations, parameters, integration, energy ledger
  run_phase_controls.py: eight cases (three refinement runs, phase OFF,
      two accommodation controls and two absorption controls)
  plot_figure13.py: redraw manuscript phase-change figure
  results/: archived ON/OFF arrays and controls/convergence CSV and JSON
  figures/: generated Figure 13 PNG and PDF
  requirements.txt: dependencies

Run in this extracted folder:
  python -m pip install -r requirements.txt
  python run_phase_controls.py
  python plot_figure13.py

To redraw archived results only:
  python plot_figure13.py

The full run regenerates every case, including additional sweep arrays, and
checks the normalized energy residual against 1e-6. Table 7 values come from
results/phase_controls.csv and results/phase_convergence.csv.

Assumptions: an existing air/water-vapor cavity, two-sided Hertz-Knudsen
transfer, IAPWS saturation pressure, lumped thermal nodes, 30 microseconds,
representative pulsed heating. This calculation does not describe nucleation
or establish a quantitative fit to observed growth or optical flashes.
