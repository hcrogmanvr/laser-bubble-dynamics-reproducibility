UPDATED FIGURE 2 SCRIPT

Run with Python 3, numpy, and matplotlib:
  python figure2_updated.py

Original mechanics, energy, burst events, sampled flash and camera calculation
are preserved. The initial import typo is fixed, the flash energy variable is
renamed to avoid confusing it with absorptance A(t), and internal Step labels
are removed from the figure. Energy quantities remain labelled as proxies.

New manuscript observables:
  T(t) = exp[-(R_b(t)/R_sat)^n]                       [Eq. 3]
  A(t) = (1 - T(t)) * (1 - chi(t))                  [Eq. 4]
A denotes absorbed fraction, not absorbance or geometric area.

Create the separate three-panel companion using established parameters:
  python figure2_updated.py --r-sat-m VALUE --n VALUE --chi-csv ionization.csv
The CSV must have columns time_s,chi and cover 0 to 299.99 s inclusive.
Values between supplied times are linearly interpolated; extrapolation is rejected.
Use an ionization history belonging to this simulation and the same time grid.

The manuscript gives chi=1/(1+exp[-(Ts-Tion)/DeltaT]) and
Ts=T0+E/Cskin, but the pasted script provides no calibrated Cskin, Tion, DeltaT
or energy units. No ionization history is inferred from these missing parameters.
Alternatively --neutral-limit explicitly assumes chi=0; this is labelled as a
limiting case, not the plasma-quenching result. R_sat and n are still required.
No optical parameter values in the bundled checks are used for delivered figures.

The companion is diagnostic postprocessing: it does not feed A(t) into the
original energy source. A self-consistent absorption-feedback simulation would
change the original dynamics and requires separate model revision.

MINIMUM RADIUS
The third panel measures the minimum during each completed laser-OFF phase:
[60,90), [150,180), [240,270) s. The CSV also records the final sampled radius
in each OFF phase. This is not a newly imposed time-dependent radius bound.
The current run gives 0.5, 0.5, 0.5 micrometres: no increasing minima.
Do not describe this output as minimum-radius growth.

CURRENT RUN
577 burst triggers. Radius alternates between numerical lower and upper bounds
on almost every timestep; these events are not verified physical collapses.
The unchanged 10-ms timestep does not resolve a 2-ms flash or 5-ms exposure.
Camera averaging therefore has one sample, as in the original.
This package reproduces the supplied model; it does not establish convergence,
experimental validation, or suitability of the dynamics for publication.

OUTPUTS
Figure_2_original_dynamics: nine original panels, PDF and PNG.
Off_phase_minimum_radius: directly extracted minima, PDF and PNG.
original_states.csv; off_phase_minima.csv; diagnostics.json.
With optical inputs: Companion_optical_and_minimum_radius, PDF/PNG,
and optical_observables.csv. Figure 3 is not modified.

Manuscript source inspected: AS manuscript(2).docx, Eqs. 3, 4, 7, 8
and section 3.1. Optical numerical parameters were not established from it.
