# Equations and scope of the new dynamic coupling diagnostic

This file documents a **new diagnostic extension**, not original manuscript code and not a validated bubble cloud model. The archived reconstruction changes output pressures only after independent integrations. It contains no dynamic coupling or spacing parameter.

## Original continuous mechanical equations

The supplied Figure 2/3 sources update a spring–damping force balance:

    dR_i/dt = v_i
    m dv_i/dt = -k (R_i-R0) + 4 pi R_i^2 (P_i-Pamb) - gamma v_i

This is the equation audited here; it is not the Rayleigh–Plesset equation. Source units and calibration must be reassessed before giving its coefficients a physical interpretation. The implicit reference retains the source pressure ceiling but removes radius and speed clamps. It stops if R crosses the source trigger downward rather than ignoring a required burst jump. No crossing occurs in the reported tests.

Figure 2: R0 = 2 micrometers, m = 1e-12, k = 5e7, gamma = 3e4, Pamb = 1e5 Pa. Source power proxy = 8 times gate times forcing multiplier; energy decay time = 15 s; P = min(2e8 E, 1.5e9). Source trigger R < 1.5 micrometers and inward motion. Original numerical radius bounds = 0.5 and 30 micrometers; speed bound = 200 m/s.

Figure 3: R0 = 2 micrometers, m = (4/3) pi 1000 R0^3, k = 6e7, gamma = 4e4, Pamb = 1e5 Pa. Source proxy = 10 times min((1+4 M) gate, 3) times forcing multiplier. Energy decay time = 20 s; P = min(0.25 E / ((4/3) pi R^3), 2e9). Source trigger R < 1.2 micrometers and inward motion. Original base floor = 0.4 micrometers, ceiling = 30 micrometers, speed limit = 250 m/s. M increases at 0.005 per second during OFF and decreases by a prescribed flash-linked burn term. The no-event branch has no flash burn, so M is constant during ON and increases linearly during OFF. The optical depth diagnostics do not enter the source energy equation.

The usual gate is ON for 60 s and OFF for 30 s. All energies and their integrals from these two sources remain proxies, without inferred joule or watt calibration. All coefficients are unchanged unless a control explicitly varies them.

## New pressure-mediated extension

For fixed bubble centers at separation d_ij, use the instantaneous far-field pressure form derived from the volume flux Q_j = 4 pi R_j^2 v_j:

    p_ij = rho / (4 pi d_ij) dQ_j/dt
         = rho/d_ij (R_j^2 a_j + 2 R_j v_j^2).

Insert sum_j p_ij into the external pressure in the source force balance. The coupled accelerations satisfy the linear system:

    sum_j M_ij a_j = F_i - 4 pi R_i^2 rho sum_(j != i) 2 R_j v_j^2 / d_ij
    M_ii = m
    M_ij = 4 pi R_i^2 rho R_j^2 / d_ij,  i != j
    F_i = -k (R_i-R0) + 4 pi R_i^2 (P_i-Pamb) - gamma v_i.

The sign corresponds to neighboring positive pressure opposing expansion. All accelerations are solved simultaneously inside each right-hand-side evaluation. The coupling-off run has diagonal mass matrix m I and no interaction-pressure term. No output pressure correction is substituted for dynamic feedback.

Use N = 3 equal bubbles, equal forcing and initial conditions, rho = 1000 kg/m^3, centers [0,d,2d], and d = 40, 100 or 200 micrometers. End-to-end spacing is 2d. The maximum R/d ranges approximately 0.01–0.05 for the reported reference solutions. These spacings are diagnostic choices, not fitted experimental geometry. The window is 2 s, entirely inside the first ON gate. Source coefficients and forcing per bubble remain fixed across all paired controls.

The extension assumes spherical, nontranslating bubbles and instantaneous pressure exchange. It excludes finite propagation, near-field deformation, translation, coalescence, thermal transport, phase change and optical shielding. Instantaneous coupling is not appropriate for the earliest numerical relaxation transient when travel time d/c exceeds the transient time. This diagnostic therefore cannot establish quantitative acoustic transients, collective collapse, synchronization or flash enhancement.

## Numerical details

Use SciPy Radau with scaled state deltaR/1e-9 and v/1e-6. Restart local integration time at each laser gate switch to avoid loss of resolution when a stiff initial adjustment is below the spacing of global floating-point times. Energy is evaluated analytically over each gate interval on the no-event branch. A supplied analytic Jacobian is used for the single-bubble reference. The coupled runs use numerical Jacobians.

For the 300 s single-bubble baseline, compare (rtol,max_step) = (1e-5,5 s), (1e-7,1 s), (1e-9,0.25 s). Radius absolute tolerances are rtol times 0.01 nm; velocity absolute tolerances are max(1e-6,10 rtol) times 1e-6 m/s. Maximum step is a ceiling, not the actual adaptive step. Saved evaluation points are independent of integration steps. The shared output grid includes 20 ms points plus logarithmically dense startup and default gate-switch samples.

Coupling pairs use rtol = 1e-8 and max_step = 0.05 s; the coupling-on run is repeated at rtol = 1e-10 and max_step = 0.01 s. Refinement differences are empirical indicators, not rigorous error bounds. Coupling-off refinement is not separately estimated here. The measured radius on/off differences are comparable to these differences and show no reliable spacing trend. Do not interpret them as a resolved physical signal.

The preserved discrete source scripts use np.arange(0,300,dt), so their final sample precedes 300 s. Continuous reference outputs include 300 s. Reported event counts refer to those original sampling conventions. Flash tails in the optimized audit use an algebraically equivalent recurrence; agreement against the untouched source is checked and recorded.
