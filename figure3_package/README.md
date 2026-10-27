# Figure 3: reproducible 12-panel extension

This package starts from the Figure 3 script supplied in this conversation. It preserves its original mechanical, energy, burst, impulse, and accretion updates. The only feedback change in the updated run is the evolving lower radius bound. The original scalar absorption variable has been renamed `drive` to distinguish it from the new absorption fraction.

## Run

Install numpy and matplotlib, then run:

    python figure3_updated.py --out results

Adjust the explicitly assumed parameters with `--tau0` and `--volume-gain`. Set `--volume-gain 0` to recover the original dynamics. Results contain PNG/PDF figures, full-resolution CSV/NPZ arrays, and a JSON comparison. CSV labels follow the source variable names; `t` is seconds, `R` and `R_floor` meters, `P` pascals, `v` meters/second, and `pressure_impulse` pascal-seconds. E, Te, and flash power use the nominal SI interpretation of the original coefficients. The source does not establish a physical mass unit for `m_np`.

## Added equations and assumptions

1. Original energy drive, unchanged:
   drive_i = min[(1 + chi_abs*m_np[i-1])*laser_i, 3].
   This multiplier can exceed one and is not an absorption fraction.

2. Added optical depth:
   tau_i = tau0*(1 + chi_abs*m_np[i-1]).
   Transmission T_i = exp(-tau_i); absorption fraction A_i = 1 - exp(-tau_i).
   These are absorption-only sample properties, with scattering and reflection neglected. Optical depth is assumed proportional to the original loading proxy at constant effective cross-sectional area. These diagnostics are not used to replace the original energy drive. Consequently this extension is not yet a self-consistent optical energy-budget model. During laser-off intervals the fractions remain defined as sample properties; actual transmitted and absorbed laser powers would be zero. No incident laser power has been added or inferred.

3. Added lower radius bound:
   R_floor[i] = Rmin0*(1 + volume_gain*m_np[i-1])^(1/3).
   This assumes an equivalent retained-volume contribution proportional to the existing nanoparticle proxy. It is applied in BOTH the volume denominator and the radius clipping step. It can decrease when the proxy decreases; monotonic irreversible growth is not assumed. It is a prescribed bound, not an independently resolved bubble minimum or proof of nanoparticle-shell geometry. Its tiny loading-to-radius lag is the original explicit update convention.

Default tau0=1 and volume_gain=1 are uncalibrated choices introduced here, not measured values or recovered manuscript parameters. volume_gain has inverse units of the original loading proxy. Replace these choices with the manuscript's actual equations and calibrated coefficients before claiming quantitative agreement. There are no experimental observations in this package.

## Layout

Six rows and two columns, ordered across each row: (a) radius and evolving lower bound, (b) laser, (c) energy, (d) temperature, (e) pressure, (f) radial velocity, (g) bursts, (h) flash power, (i) cumulative pressure impulse, (j) rupture indicator, (k) transmission fraction, (l) absorption fraction. A small inset in panel (a) resolves the lower bound. Both baseline and updated figures have identical layouts; the baseline has constant minimum radius but the same added optical diagnostic equations.

## Verification and observed differences

The script checks finite arrays, T+A=1, radius >= lower bound, and exact equality of every original state array in a zero-growth regression. The supplied default run gives 577 bursts in both versions. The evolving bound ends at approximately 0.449740 micrometers, starting at 0.4 micrometers. Pressure and cumulative impulse change; the energy, temperature, velocity, bursts, flash, nanoparticle proxy, and rupture arrays remain exactly identical in this run. See comparison.json for the maximum absolute changes.

## Limitations exposed by running the original code

- At every step after initialization, the radius is at one of its two clipping bounds and velocity is at its +/-250 m/s clipping bound. A speed of 250 m/s integrated over 0.01 s gives a raw displacement of 2.5 m, vastly exceeding the 30-micrometer maximum radius. The plotted cycles therefore do not establish resolved growth and collapse. Dense full-duration curves must not be interpreted as physically resolved oscillations.
- Pressure hits the 2 GPa cap for about 88.4% of updated time steps.
- The flash timescale is 0.002 s but dt=0.01 s; the exponential pulses are undersampled. Their sampled numerical integrals cannot be used as the intended radiated burst energies without a temporal-integration correction or a resolved grid.
- Pressure impulse is cumulative from t=0. The declared impulse_window is unused. This extension deliberately preserves that original calculation.
- Temperature and pressure are calculated before burst energy subtraction, while the saved E array is post-burst. The original intra-step ordering is preserved and is not an equilibrium EOS consistency check.
- The nanomass, energy-input, stiffness, and damping coefficients require physical justification and consistent units before quantitative interpretation. The original code's 'physical EOS' comment alone is not validation.

This package completes the requested figure/code extension. It does not resolve these inherited numerical and physical limitations. Numerical convergence and physical calibration remain necessary before publication claims about bubble collapse, temperature, rupture, or optical emission.

## Suggested descriptive caption

Figure 3. Outputs of the supplied single-bubble model with added transmission, absorption, and an accretion-dependent lower radius bound. Panels (a)-(l) show radius with its lower bound, laser modulation, stored energy, electron temperature, pressure, radial velocity, burst events, flash power, cumulative pressure impulse, rupture indicator, transmission fraction, and absorption fraction. Optical fractions are derived from an assumed absorption-only attenuation law; the lower bound follows an assumed volume-loading relationship. Added coefficients are uncalibrated. The radius, velocity, and pressure retain the original numerical clipping rules; the results do not establish resolved physical collapse dynamics.
