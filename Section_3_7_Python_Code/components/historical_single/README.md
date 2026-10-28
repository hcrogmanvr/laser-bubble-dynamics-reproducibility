Repository note: This is a component archive. The current paper-wide entry point, figure numbering and source provenance are documented in the root README and docs/FIGURE_INDEX.md. Original Figure 4 is now recovered in components/original_figure4.

# Completed reviewer controls and red revision

The supplied Figure 2 and Figure 3 scripts reproduce 577 events, but the continuous force equations do not produce collapse when integrated with an implicit solver without radius and velocity clipping. The original event interpretation is not supported by this audit. The accompanying Word report gives the actual controls, figures, quantitative findings, red replacement manuscript text and a reviewer response draft.

## Reproduce

Use Python 3.11 or newer with the package versions in `requirements.txt`:

    python -m pip install -r requirements.txt
    python code/run_controls.py
    python code/verify_results.py
    python code/build_report.py

The runner writes CSV/JSON tables and compressed arrays to `results/`; `build_report.py` regenerates figures and the Word report from those data. No internet connection, input fitting or experimental data generation is required to execute the controls. Preserve result outputs before rerunning if the current recorded run must remain immutable. `requirements.txt` records the environment used, not a blanket compatibility claim.

## Files

- `source_originals/`: untouched Figure 2/3 sources and the separately archived reconstruction, with its original provenance statement.
- `code/run_controls.py`: parameter sweeps, original-source reproduction, implicit references, reconstruction audits and the explicitly new dynamic coupling diagnostic.
- `code/verify_results.py`: source hashes, reproduction checks, an independent static equilibrium check and saved comparison checks.
- `code/build_report.py`: figures and Word report generation.
- `results/`: CSV tables, JSON records and NPZ arrays; metadata record original hashes, runtime versions and result scope.
- `figures/`: scientific plots generated from the completed controls.
- `dynamic_coupling_equations.md`: source force balance, new extension, parameters, solver settings and limitations.
- `Reviewer_Controls_Results_and_Red_Revision.docx`: reviewed report and red manuscript replacement text.

## What was completed

Both figure sources: six timesteps (0.01 to 0.0003125 s), factors 0.5/1/2 for radius floor, radius ceiling, speed limit, pressure ceiling and event cooldown, five forcing factors (0.5/0.75/1/1.25/1.5), and four fixed-peak-source ON/OFF schedules (30/60,45/45,60/30,75/15 s), all over 300 s. Figure 3 also has three optical diagnostic null controls. The unclipped continuous equations were separately checked at three solver settings and run through forcing, pressure ceiling and duty-cycle sweeps. All reported reference cases have zero events.

The archived reconstruction was tested for quenching strengths 0 and 0.65 over 1 s, and coupling strengths 0 and 0.08 for three bubbles with seed 7 over 0.1 s, both at 0.1 ms. Reconstruction controls are not validation of the original figure sources. Its mean-field coupling alters pressure output only; its trajectory difference is exactly zero.

A new, explicitly distinguished pressure-mediated extension was solved for three equal bubbles at 40/100/200 micrometer center spacings over 2 s, with paired coupling on/off and coupling-on refinement. It has no collapse events and no radius spacing effect clearly resolved above refinement differences. This is a diagnostic extension, not recovery of original multi-bubble manuscript code. Fixed peak input and fixed geometry comparisons do not establish fixed deposited-energy equivalence.

## Interpretation

At every tested discrete timestep, speed and radius are clipped on every update. Event counts vary with the prescribed waiting time, and maxima follow the imposed radius ceiling. The implicit references stay near the initial 2 micrometer radius. The baseline fine-reference displacements are approximately 1.510 and 1.678 nm. These calculations solve the source spring–damping equations, not Rayleigh–Plesset dynamics.

The source energy coefficients are treated as proxies, without inferred SI energy calibration. The reconstruction's names containing `J` are nominal labels from that code, not experimental calibration. Neither original Figure 2/3 energy equation contains quenching. Diagnostic absorptance is not deposited energy feedback. The reconstruction suppresses energy when quenching is enabled but remains pressure capped and radius bounded, with zero flashes in the reported 1 s window.

No measured optical flash comparison, spectrum, phase transition validation, validated two-temperature EOS, or demonstrated collective flash enhancement is supplied by these tests. The red text in the report must replace incompatible positive claims throughout a revised manuscript before it is sent to reviewers. It is not an endorsement of submission readiness.

SciPy method reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
