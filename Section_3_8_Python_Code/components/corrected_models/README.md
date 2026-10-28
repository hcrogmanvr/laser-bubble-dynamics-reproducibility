Repository note: This is a component archive. The current paper-wide entry point, figure numbering and source provenance are documented in the root README and docs/FIGURE_INDEX.md. Original Figure 4 is now recovered in components/original_figure4.

# Corrected numerical models and verification results

This package contains executed numerical corrections and controls. It does not claim experimental plasma or flash validation. The original Figure 2–3 spring–damping equations still produce no collapse with implicit integration; changing a solver cannot make those equations reproduce the observed phenomenon. `results/source2_3_recheck.json` records fresh 300 s reruns.

## Run

Python 3.12; numpy, scipy, matplotlib. From any directory:

    python code/run_suite.py
    python code/recheck_sources.py
    python code/make_figures.py

The suite takes about two minutes on the authoring host. Output arrays preserve all saved states and event jumps. No random draws are reselected to obtain events. Seed 2 and the original eight initial positions/radii are retained for the corrected ensemble.

## Ensemble correction scope

`corrected_ensemble.py` uses the dimensionally correct incompressible radial acceleration, partial-pressure work, equal/opposite thermal exchange, ambient heat exchange, and (variant 6) effective Saha ionization energy. Adaptive Radau integrates the equations and restarts at exact envelope-gate edges and threshold events. There is no acceleration, velocity, radius, pressure, temperature or camera saturation. Domain events STOP rather than clip: radius 5 nm and speed Mach 0.1; tests move these guards without changing the solution because no guard is reached.

Original numerical source inputs Ce=3e-12 J/K, Ci=2e-11 J/K, local input scale 2 mW, density inventory 1e25 m^-3 at radius 0.2 um, seed and geometry are retained as effective assumptions. These are not measured cavity heat capacities, water-plasma properties or deposited optical power. Initial states are not mechanical equilibria; their short contraction produces the surviving five threshold crossings. They are not repeated laser-induced plasma flashes.

Variant 5 retains logistic activation as an effective pressure switch, not physical ionization. Variant 6 retains a single-stage effective 12.6 eV Saha mapping, now evaluated with a stable root and included in the energy budget. An equilibrium Saha closure at Te while Ti is evolved separately is an assumption, not general two-temperature chemical equilibrium. Hydrogen (13.5984 eV) is used only as an independent algebra benchmark; hydrogen is not identified as the experimental species.

An event removes 35% of positive excess thermal/ionization energy and 25% of liquid potential-flow kinetic energy. Temperatures are solved for consistently and velocity magnitude is reduced; there is no positive temperature kick. Eight percent of the removed energy is assigned to a radiative-energy proxy and the remainder to an unresolved external loss. This closes the jump ledger but does not derive an emission mechanism, intrinsic pulse duration or spectrum. The unbudgeted delayed pressure pulses, post-event emission multiplication and camera saturation from the earlier source are excluded. Coupling remains the source's instantaneous row-normalized pressure heuristic; its work is explicitly tracked as an external exchange and is not a conservative wave-field calculation.

Variant 6 optionally uses a single-front-interface Drude/slab optical accounting test at 1064 nm, collision frequency 1e14 s^-1, assumed path 2R and a unit-index surrounding medium. Plasma absorption plus transmitted light absorbed by the effective absorber both enter the heat budget. Reflected and unabsorbed transmitted energy leave. This is a new conditional optical closure, not a measured bead cross section or a replacement for scattering in the actual microscope. Variant 5's logistic activation is never used as an electron density in that calculation. Cold repaired runs have negligible Saha ionization, so Drude ON/OFF has no resolved effect; no plasma-quenching claim follows.

## Phase-change model

`phase_change.py` describes a separate pre-existing spherical air/water-vapor cavity of initial radius 5.44 um and an effective heated interface/absorber node with the heat capacity of the assumed 1.57 um bead. It is not a fitted replacement for the minute-scale traces. Vapor mass follows a two-sided Hertz–Knudsen flux with IAPWS saturation pressure. Incoming/outgoing vapor enthalpy uses an upwind temperature and the same latent-energy reference in both node and cavity balances. Gas pressure–volume work, viscous dissipation, liquid kinetic energy, surface energy and ambient displacement work appear in the energy ledger.

Assumed effective material inputs: bead density 2500 kg/m^3 and heat capacity 700 J/kg/K; gas conductivity 0.026 W/m/K; external conductance 4 pi (0.6 W/m/K) a. Identifying the heated absorber and interface as one temperature is an unvalidated lumped approximation; spatial thermal transport is not resolved. The pulse train uses representative 250 mW, 100 kHz and 100 ns values, not recovered trace-specific repetition settings. Absorption fractions 1e-5–1e-3 and accommodation 1e-4–1e-2 are sensitivity inputs, not fitted or measured coefficients. Model-domain guards stop at the bead radius or outside 273.16–600 K, or at nonpositive vapor mass. No guard fires in the reported cases. This model addresses evaporation/condensation in an existing cavity, not nucleation, nonequilibrium water chemistry or plasma emission.

## Verification versus validation

The executed checks establish algebra, energy accounting and numerical refinement for the stated equations. They do not establish parameter accuracy or experimental predictive skill. The Rayleigh comparison is with the analytic finite-radius collapse integral at R/R0=0.2; it is not validation of heated gas/vapor collapse. Saha and Drude checks establish their algebraic implementation, not applicability to the experimental mixture.

Quantitative flash prediction and experimental flash validation are outside the study scope. No calibrated flash trace, impulse response, absolute intensity or spectrum is present in the supplied records. `results/experimental_comparison_status.csv` documents the distinction between supported growth characterization and out-of-scope optical validation. The diameter workbook supplies growth histories, not flash calibration. The phase model runs for 30 microseconds and is not quantitatively compared with seconds-long growth regressions.

## Primary formulation records

- Rayleigh and bubble mechanics: https://brennen.caltech.edu/INTMult/multbook.pdf
- Phase transition and compressible dynamics: https://doi.org/10.1017/jfm.2024.954
- Vapor compression experiments: https://doi.org/10.1063/5.0200361
- IAPWS saturation release: https://www.iapws.org/relguide/Supp-sat.html
- Saha convention: https://goldbook.iupac.org/terms/view/08542
- Hydrogen benchmark ionization energy: https://physics.nist.gov/cgi-bin/Ionization/atom.php?element=H
- Collisional-plasma dielectric relation: https://farside.ph.utexas.edu/teaching/em/lectures/node103.html

See the root RELEASE_SHA256SUMS.json for the packaged code/data manifest. Old supplied source files are retained separately from corrected code.
