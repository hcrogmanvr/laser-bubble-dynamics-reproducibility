SECTION 3.7 PYTHON CODE

Completed numerical controls and mechanism tests

Python 3.12; install from this folder:
  python -m pip install -r requirements.txt

Verify archived results and regenerate Figures 9–11:
  python run_section_3_7.py --mode verify

Recompute all Section 3.7 controls, including the matched original-ensemble coupling checks:
  python run_section_3_7.py --mode full

Source mapping:
3.7.1–3.7.4 / Tables 3–4: components/historical_single/code/run_controls.py
3.7.5 / Table 5: components/historical_ensemble/diagnostics/coupling_check.py
Source-output exports: components/historical_ensemble/reproduce_defaults.py
Figures 9–11: scripts/plot_control_figures.py

Numerical summaries and archived arrays are included in the corresponding results folders. Table 5's unrounded values are in diagnostics/coupling_results/matched_coupling_controls.json. Default plotting does not rerun the scientific sweeps; use --mode full to recompute. Full sweeps can take several minutes.

The original scripts are preserved byte-for-byte. The top-level runner is a separate packaging utility. Do not rename the components directories: relative paths are intentional. Quenching controls apply to the explicitly identified reconstruction, and the original coupling controls do not correct the historical model defects.
