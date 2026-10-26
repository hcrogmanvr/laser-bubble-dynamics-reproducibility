# Laser bubble dynamics reproducibility repository

**Compact upload edition:** see COMPACT_EDITION.md for the files omitted to keep this ZIP below 25 MB. All code and inputs needed to regenerate them are retained.

Code, source data, numerical controls and figure reproduction for **Laser-Driven Bubble Dynamics and Optical-Emission Proxies: Numerical Verification and Experimental Growth Characterization**.

This repository follows the manuscript revision that retains the **original submitted Figure 4**, documented in Section 2.2. It includes the exact manuscript images for Figures 1–13 (with Figure 8A and 8B separately), the experimental diameter workbook, original source algorithms, corrected conditional models, phase-change calculations, and the numerical records behind Tables 1–7. All 108 native Word mathematical expressions are indexed and preserved as Office Math XML.

The contribution is numerical verification and descriptive experimental growth characterization. Prescribed detector signals, thermal/ionization proxies and pressure-threshold flags are not presented as calibrated experimental plasma, flash or biological-damage predictions.

## Quick start

Python 3.12 was used. The tested dependency versions are in requirements-lock.txt.

    python -m venv .venv
    # Linux/macOS:
    source .venv/bin/activate
    # Windows PowerShell:
    # .venv\Scripts\Activate.ps1
    python -m pip install -r requirements-lock.txt
    python reproduce.py --mode verify

The verification mode redraws every manuscript figure from the supplied source/arrays, rebuilds numerical tables, calculates reference quantities, and checks provenance and stated numerical results. Files are written to outputs/; logs identify any failed step.

## Recalculate everything

    python reproduce.py --mode full

This reruns:
- the original Figure 2–3 source reproduction, timestep/bound/event-rule/forcing controls, unclipped continuous references, and declared mechanism diagnostics;
- the original Figure 5–6 defaults, matched coupling controls, restricted event-phase observations, and electron-density observer;
- all 26 corrected-ensemble cases, eight phase-change cases, independent closure benchmarks, and fresh continuous Figure 2–3 rechecks;
- the author-confirmed original Figure 4 source;
- every plot, numerical table and reference calculation, followed by verification.

Run time depends on the computer; the full suite takes several minutes on the preparation environment. No random seed, input or event rule is reselected to obtain a preferred result.

## Reproduce one figure

    python reproduce.py --figure 4
    python reproduce.py --figure 4 --figure4-run

The first command plots the archived original Figure 4 arrays; the second reruns its exact submitted algorithm first.

    python reproduce.py --mode full --figure 13
    python reproduce.py --figure 8B

Figures 9–11 share one control-family renderer; Figures 12–13 share one corrected-model renderer. Selecting one may therefore generate its companion figures. See docs/FIGURE_INDEX.md and docs/FIGURE_INDEX.csv for every source, command, caption and interpretation. The additional reference/optical check plot Figure S1 is generated automatically and indexed in docs/SUPPLEMENTARY_INDEX.json.

## Repository layout

| Path | Contents |
| --- | --- |
| components/historical_single/ | Original Figure 2–3 sources and executed discrete/continuous controls supporting Figures 9–11 |
| components/original_figure4/ | Confirmed submitted Figure 4 source, exact execution/export wrapper, full arrays and parameter/run records |
| components/historical_ensemble/ | Original Figure 5–6 source versions, returned arrays, realized geometry, coupling/density observers |
| components/corrected_models/ | Corrected effective-reservoir solvers, separate finite-rate air/water-vapor model, tests and Figures 12–13 |
| data/experimental/ | Original diameter workbook and archived descriptive fit outputs |
| data/reference_arrays/ | Immutable source-array references used by verification |
| figures/manuscript/ | Exact image bytes extracted from the latest manuscript, with hashes |
| figures/editable/ | Figure 1 SVG and editable PowerPoint with native equations |
| scripts/ | Plotting, table export, scalar calculations, diagram rebuilding and repository verification |
| configs/ | Figure-specific parameters, geometry and case definitions |
| docs/tables/ | Exact manuscript Tables 1–7 as CSV; no production-template table |
| docs/equations/ | Native Office Math XML, indexed in docs/equation_index.json |
| outputs/ | Regenerated figures, plotted coordinates, fitted quantities, calculated tables and logs |
| validation/ | Verification report and environment/source manifests from the prepared release |

## Important source distinctions

**Original Figure 4.** The source is multibubble_fixed_full.py: N=8, seed=4, 1 ms, 100 ns steps, and 60/20 µs gating. It produces eight first-update flags, followed by radii at the imposed 50 nm floor and an initial camera maximum of 0.0028. It is not the corrected variant-6 initialization plot that temporarily replaced it. Its proxy energy, pressure, temperature, clipping and event rules are preserved for reproduction.

**Original versus corrected ensembles.** Figures 5–6 retain their historical source algorithms and documented defects. Figure 12 uses the separately corrected radial/energy/event model. Numerical verification of the corrected models does not validate the historical plasma interpretation.

**Figure 7.** The executable redraw uses the archived extended-source t and rupt_frac arrays. That legacy name denotes pressure-threshold flags, not measured membrane rupture. The original image is preserved, but its exact generation history is not independently established. A source-derived redraw is supplied with explicit provenance rather than asserted to be a pixel-identical reconstruction.

**Figure 1.** The editable diagram represents historical implementation scope. Its dashed extensions should be read with the later corrected and separate phase-change formulations in the manuscript. Rendering its SVG is reproducible; rebuilding its editable geometry is optional with python scripts/rebuild_diagram.py.

**Experimental growth.** Diameter units are µm, as confirmed by the author. The workbook's mm header is treated as a labeling error. Recorded Y-error values of 2 are preserved but their statistical meaning is undocumented: they are not plotted as pointwise error bars or used as fit weights. The first expansion segment ends at the first recorded maximum. Slopes and conventional OLS confidence intervals are descriptive; serial observations and repeated cycles are not independent experimental groups. No synthetic physical-model overlay is introduced.

## Tables and equations

    python scripts/export_calculated_tables.py
    python scripts/reference_calculations.py

Tables 2–7 are recalculated from their identified records. Table 1 is the published/setup/model-input provenance table, not a fitted numerical output. Solver residuals may vary below tolerance across platforms; verification tests their stated acceptance threshold and separately checks physical/control quantities at manuscript precision.

See docs/CALCULATION_COVERAGE.md for the relationship between manuscript equation groups and executable source versions. Some early equations are conceptual mappings with unspecified material parameters; they are indexed rather than assigned invented inputs or presented as additional executed simulations.

## Citation and reuse

Citation metadata are in CITATION.cff. The manuscript has no repository DOI assigned in this package; no identifier or online hosting address is invented. The six manuscript authors are credited there. See NOTICE.md for licensing scope and the retained license on the earlier reconstructed component.

RELEASE_SHA256SUMS.json freezes the tracked bytes of this prepared release. Rerunning calculations can change binary-container metadata or solver residuals below tolerance; numerical/source verification is performed by reproduce.py rather than demanding identical PDF/NPZ container bytes.

The source revision is identified by filename and SHA-256 in docs/manuscript_reference.json. The full manuscript and the motivating journal PDF are not republished in this code/data package. Exact source images, captions, numbered tables and native formulas are included for reproduction and review.

## Upload to GitHub

The repository is prepared locally; no GitHub repository has been created or published by this package. See docs/PUBLISHING.md for import instructions and the concrete repository name/description. A separate Git bundle can restore the prepared Git history, if provided with this release.
