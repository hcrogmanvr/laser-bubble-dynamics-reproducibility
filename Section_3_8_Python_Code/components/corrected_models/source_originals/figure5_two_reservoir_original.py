import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

# ----------------------------
# Utilities
# ----------------------------
def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))

@dataclass
class CameraModel:
    tau: float = 20e-6          # camera decay time (s)
    gain: float = 1e-3          # scaling from photon flux proxy -> intensity
    saturate: float = 2.0       # saturation level (arb)
    baseline: float = 0.0       # baseline intensity (arb)

    def step(self, I_prev, photon_rate, dt):
        # Simple first-order sensor response: dI/dt = -(I - baseline)/tau + gain*photon_rate
        dI = (-(I_prev - self.baseline) / self.tau + self.gain * photon_rate) * dt
        I = I_prev + dI
        return np.clip(I, 0.0, self.saturate)

# ----------------------------
# Physics blocks (2T + EOS + ionization)
# ----------------------------
@dataclass
class PlasmaEOS2T:
    # Gas/plasma constants (simple ideal mixture with ionization fraction)
    mu: float = 18.0            # mean molar mass (water vapor approx) [g/mol] (proxy)
    Z: float = 1.0              # mean ion charge state
    Te_ion: float = 1.2e4       # ionization threshold temperature [K]
    dTe_ion: float = 2.0e3      # smoothness of ionization turn-on [K]
    P_cap: float = 2.0e9        # pressure cap [Pa] (safety)
    kB: float = 1.380649e-23
    NA: float = 6.02214076e23

    def ion_fraction(self, Te):
        # Smooth logistic “ionization” switch (proxy for real Saha/CR model)
        # x_e ~ 0 below Te_ion, ~1 above
        return sigmoid((Te - self.Te_ion) / self.dTe_ion)

    def pressure(self, n0, Ti, Te):
        # Ideal mixture: P = n_i kB Ti + n_e kB Te
        # n_i ~ n0 ; n_e = x_e * Z * n0
        x = self.ion_fraction(Te)
        P = n0 * self.kB * Ti + (x * self.Z * n0) * self.kB * Te
        return np.minimum(P, self.P_cap), x

@dataclass
class TwoTempEnergyModel:
    Tamb: float = 300.0
    # Effective heat capacities per bubble (lumped, proxy)
    Ce: float = 3.0e-12         # J/K  (electron subsystem effective)
    Ci: float = 2.0e-11         # J/K  (ion/vapor subsystem effective)
    tau_ei: float = 2.0e-6      # e-i equilibration time (s)
    tau_cool_e: float = 20e-6   # electron cooling (s) (radiation+conduction proxy)
    tau_cool_i: float = 80e-6   # ion cooling (s)
    # Absorption -> electron heating fraction
    fe_abs_to_e: float = 0.9

    def step(self, Te, Ti, Pabs, dt):
        # Electron subsystem
        dTe_abs = (self.fe_abs_to_e * Pabs / self.Ce) * dt
        dTe_ei  = (-(Te - Ti) / self.tau_ei) * dt
        dTe_cool = (-(Te - self.Tamb) / self.tau_cool_e) * dt
        Te_new = Te + dTe_abs + dTe_ei + dTe_cool

        # Ion subsystem
        dTi_abs = ((1.0 - self.fe_abs_to_e) * Pabs / self.Ci) * dt
        dTi_ei  = ((Te - Ti) / self.tau_ei) * (self.Ce / self.Ci) * dt
        dTi_cool = (-(Ti - self.Tamb) / self.tau_cool_i) * dt
        Ti_new = Ti + dTi_abs + dTi_ei + dTi_cool

        # Floors
        Te_new = np.maximum(Te_new, 50.0)
        Ti_new = np.maximum(Ti_new, 50.0)
        return Te_new, Ti_new

# ----------------------------
# Bubble dynamics (stable Rayleigh-Plesset proxy + coupling)
# ----------------------------
@dataclass
class BubbleDynamics:
    rho: float = 1000.0         # kg/m^3 (water)
    mu: float = 1.0e-3          # Pa*s (viscosity)
    sigma: float = 0.072        # N/m surface tension
    P0: float = 1.0e5           # Pa ambient
    c: float = 1500.0           # m/s sound speed (water)
    # collapse / numerical safety
    R_min: float = 5e-8         # 0.05 µm
    R_max: float = 3e-5         # 30 µm
    dR_max: float = 250.0       # m/s clamp
    accel_cap: float = 2e10     # m/s^2 clamp (safety)

    def step(self, R, dR, P_drive, P_couple, dt):
        # Total pressure inside minus outside
        # Outside pressure = P0 + viscous + surface terms
        # Use a stable “RP-like” second-order ODE:
        # d2R ~ (P_in - P_out)/rho - (3/2)(dR^2)/R  - (4mu/rho)(dR/R) - (2sigma)/(rho R)
        # This is simplified but bounded.

        R_safe = np.maximum(R, self.R_min)

        P_out = self.P0 + (4.0 * self.mu * dR / R_safe) + (2.0 * self.sigma / R_safe)
        P_in = P_drive + P_couple

        rhs = (P_in - P_out) / self.rho - 1.5 * (dR * dR) / R_safe
        rhs = np.clip(rhs, -self.accel_cap, self.accel_cap)

        dR_new = dR + rhs * dt
        dR_new = np.clip(dR_new, -self.dR_max, self.dR_max)

        R_new = R + dR_new * dt
        R_new = np.clip(R_new, self.R_min, self.R_max)

        return R_new, dR_new

# ----------------------------
# Laser + absorption + burst->photons
# ----------------------------
@dataclass
class LaserModel:
    # Duty-cycle modulation
    period: float = 80e-6
    on_frac: float = 0.75
    # Power absorbed per bubble (proxy); scaled by absorption factor
    Pin: float = 2.0e-3         # W effective per bubble (NOT the macroscopic laser power)
    abs_base: float = 0.25      # baseline absorption fraction
    skin_depth: float = 1.0e-7  # m
    abs_R_scale: float = 2.0e-6 # m (bigger bubbles absorb more in this proxy)

    def laser_on(self, t):
        phase = (t % self.period) / self.period
        return 1.0 if phase < self.on_frac else 0.0

    def absorbed_power(self, laser_flag, R):
        # Bounded absorption model:
        # - increases with R (more intercept)
        # - limited by skin-depth factor (thin absorption at small R)
        # - total fraction clipped to [0, 0.95]
        frac_R = 1.0 - np.exp(-R / self.abs_R_scale)
        frac_skin = 1.0 - np.exp(-R / self.skin_depth)
        frac = self.abs_base * frac_R * (0.3 + 0.7 * frac_skin)
        frac = np.clip(frac, 0.0, 0.95)
        return laser_flag * self.Pin * frac, frac

@dataclass
class BurstModel:
    R_collapse: float = 1.5e-7      # collapse threshold (m)
    quench_frac: float = 0.75       # fraction of stored thermal energy removed on burst
    burst_to_rad: float = 0.08      # fraction of burst energy emitted as light
    wavelength_nm: float = 550.0    # representative visible wavelength for photon energy
    burst_cool_Te_kick: float = 4e4 # extra Te kick during burst (K), proxy

    def photon_energy(self):
        # E_ph = hc/lambda
        h = 6.62607015e-34
        c = 299792458.0
        lam = self.wavelength_nm * 1e-9
        return h * c / lam

# ----------------------------
# Multi-bubble coupling
# ----------------------------
@dataclass
class CouplingModel:
    # Simple pressure/shock coupling: each bubble feels weighted sum of neighbors' pressure fluctuations
    k_couple: float = 0.25
    length: float = 2.5e-6       # characteristic decay length (m)
    max_couple: float = 5e8      # Pa cap

    def build_weights(self, positions):
        # positions: (N,3)
        N = positions.shape[0]
        d = positions[:, None, :] - positions[None, :, :]
        r = np.linalg.norm(d, axis=2)
        # Avoid self
        r = np.maximum(r, 1e-12)
        W = np.exp(-r / self.length) / r
        np.fill_diagonal(W, 0.0)
        # Normalize rows so k_couple acts like a gain
        row_sum = W.sum(axis=1, keepdims=True) + 1e-30
        Wn = W / row_sum
        return Wn

    def coupled_pressure(self, W, P_drive):
        # Use deviations from mean as a “shock-like” coupling
        dP = P_drive - P_drive.mean()
        P_c = self.k_couple * (W @ dP)
        return np.clip(P_c, -self.max_couple, self.max_couple)

# ----------------------------
# Main simulation
# ----------------------------
def run_multibubble_2T_plasma(
    N=8,
    t_end=1.0e-3,
    dt=1.0e-7,
    stride=5,
    seed=2
):
    rng = np.random.default_rng(seed)

    # Models
    eos = PlasmaEOS2T(Te_ion=1.2e4, dTe_ion=2.0e3, P_cap=2.0e9)
    thermal = TwoTempEnergyModel(Ce=3e-12, Ci=2e-11, tau_ei=2.0e-6, tau_cool_e=20e-6, tau_cool_i=80e-6)
    dyn = BubbleDynamics(R_min=5e-8, R_max=3e-5, dR_max=250.0)
    laser = LaserModel(period=80e-6, on_frac=0.75, Pin=2.0e-3, abs_base=0.35)
    burst = BurstModel(R_collapse=1.5e-7, quench_frac=0.75, burst_to_rad=0.08, wavelength_nm=550.0)
    cam = CameraModel(tau=20e-6, gain=1e-3, saturate=2.0)
    couple = CouplingModel(k_couple=0.25, length=2.5e-6, max_couple=5e8)

    # Bubble positions (random cluster)
    pos = rng.normal(scale=2.0e-6, size=(N, 3))
    W = couple.build_weights(pos)

    # Initial states
    R = np.full(N, 1.8e-7) + rng.normal(scale=0.2e-7, size=N)
    R = np.clip(R, dyn.R_min, dyn.R_max)
    dR = np.zeros(N)

    Te = np.full(N, thermal.Tamb)
    Ti = np.full(N, thermal.Tamb)

    # “Mass” in bubble (proxy): choose number of neutral particles fixed
    # n0 depends on volume, so pressure rises when R shrinks (physically consistent trend)
    # Choose N_particles so that at R~0.2 µm, n0 ~ 1e25 m^-3 (dense vapor proxy)
    R_ref = 2.0e-7
    V_ref = (4/3) * np.pi * R_ref**3
    n_ref = 1.0e25
    N_particles = n_ref * V_ref  # per bubble, fixed

    # Energy bookkeeping (store a scalar internal energy proxy per bubble)
    # Track E as "thermal energy above ambient" in the two subsystems
    E = np.zeros(N)

    # Camera intensity (global sum of bubble flashes)
    I_cam = 0.0

    # History (downsampled)
    steps = int(np.ceil(t_end / dt))
    store_n = steps // stride + 1

    t_hist = np.zeros(store_n)
    Rm_hist = np.zeros(store_n)
    Em_hist = np.zeros(store_n)
    Pm_hist = np.zeros(store_n)
    Tem_hist = np.zeros(store_n)
    I_hist = np.zeros(store_n)
    laser_hist = np.zeros(store_n)
    bursts_hist = np.zeros(store_n)

    # Burst detection
    total_bursts = 0
    hv = burst.photon_energy()

    si = 0
    for k in range(steps + 1):
        t = k * dt

        # Laser modulation (scalar flag)
        L = laser.laser_on(t)

        # Absorbed power per bubble
        Pabs, abs_frac = laser.absorbed_power(L, R)

        # Two-temperature update
        Te, Ti = thermal.step(Te, Ti, Pabs, dt)

        # Internal energy proxy (above ambient)
        E = thermal.Ce * (Te - thermal.Tamb) + thermal.Ci * (Ti - thermal.Tamb)
        E = np.maximum(E, 0.0)

        # Number density from fixed particle count and current volume
        V = (4/3) * np.pi * np.maximum(R, dyn.R_min)**3
        n0 = N_particles / V

        # Plasma EOS pressure + ionization fraction
        P_drive, xion = eos.pressure(n0, Ti, Te)

        # Multi-bubble coupling pressure
        P_couple = couple.coupled_pressure(W, P_drive)

        # Bubble dynamics step
        R_prev = R.copy()
        dR_prev = dR.copy()
        R, dR = dyn.step(R, dR, P_drive, P_couple, dt)

        # Collapse-triggered burst: downward crossing of R_collapse
        crossed = (R_prev > burst.R_collapse) & (R <= burst.R_collapse) & (dR_prev < 0)

        burst_energy = np.zeros(N)
        if np.any(crossed):
            # Define burst energy as a fraction of stored internal energy + collapse kinetic proxy
            # Kinetic proxy ~ 0.5 rho (4πR^3/3) dR^2
            Vc = (4/3) * np.pi * np.maximum(R_prev[crossed], dyn.R_min)**3
            K = 0.5 * dyn.rho * Vc * (dR_prev[crossed] ** 2)

            Eb = 0.35 * E[crossed] + 0.25 * K  # stable, bounded
            burst_energy[crossed] = Eb

            # Quench (remove energy)
            E[crossed] *= (1.0 - burst.quench_frac)

            # Apply a Te "kick" to push ionization transiently (proxy for discharge)
            Te[crossed] = np.minimum(Te[crossed] + burst.burst_cool_Te_kick, 2.0e5)

            # Recompute Ti from remaining energy budget (keep consistent-ish)
            # Allocate remaining E mainly to ions
            E_rem = np.maximum(E[crossed], 0.0)
            Ti[crossed] = thermal.Tamb + (E_rem / (thermal.Ci + 1e-30))
            Ti[crossed] = np.clip(Ti[crossed], 50.0, 2.0e5)

            total_bursts += int(crossed.sum())

        # Convert burst energy -> photon rate -> camera
        # Photon count ~ (eta_rad * Eb) / hv ; treat as delta-like injection spread over dt
        photon_rate = (burst.burst_to_rad * burst_energy.sum()) / (hv * dt + 1e-30)
        I_cam = cam.step(I_cam, photon_rate, dt)

        # Store downsampled
        if (k % stride) == 0:
            t_hist[si] = t
            Rm_hist[si] = R.mean()
            Em_hist[si] = E.mean()
            Pm_hist[si] = P_drive.mean() * 1e-6   # MPa
            Tem_hist[si] = Te.mean()
            I_hist[si] = I_cam
            laser_hist[si] = L
            bursts_hist[si] = 1.0 if np.any(crossed) else 0.0
            si += 1

    out = {
        "t": t_hist[:si],
        "R_mean": Rm_hist[:si],
        "E_mean": Em_hist[:si],
        "P_mean_MPa": Pm_hist[:si],
        "Te_mean": Tem_hist[:si],
        "I_cam": I_hist[:si],
        "laser": laser_hist[:si],
        "burst_flag": bursts_hist[:si],
        "total_bursts": total_bursts
    }
    return out

# ----------------------------
# Plotting
# ----------------------------
def plot_results(out, title="Multi-bubble 2T plasma EOS + ionization + camera"):
    t = out["t"]
    # Convert time axis to ms if long enough
    t_ms = t * 1e3

    fig, ax = plt.subplots(6, 1, figsize=(12, 10), sharex=True)

    ax[0].plot(t_ms, out["R_mean"] * 1e6)
    ax[0].set_ylabel("Mean R (µm)")
    ax[0].set_title(title)

    ax[1].plot(t_ms, out["laser"])
    ax[1].set_ylabel("Laser (0/1)")

    ax[2].plot(t_ms, out["E_mean"])
    ax[2].set_ylabel("Mean E (J)")

    ax[3].plot(t_ms, out["P_mean_MPa"])
    ax[3].set_ylabel("Mean P (MPa)")

    ax[4].plot(t_ms, out["Te_mean"])
    ax[4].set_ylabel("Mean Te (K)")

    ax[5].plot(t_ms, out["I_cam"], label="Camera flash")
    # show burst events as faint markers
    burst_mask = out["burst_flag"] > 0
    if np.any(burst_mask):
        ax[5].vlines(t_ms[burst_mask], 0, np.max(out["I_cam"]) + 1e-30, alpha=0.2, label="Burst events")
    ax[5].set_ylabel("I_cam (arb)")
    ax[5].set_xlabel("Time (ms)")
    ax[5].legend(loc="upper right")

    plt.tight_layout()
    plt.show()

# ----------------------------
# Run
# ----------------------------
if __name__ == "__main__":
    out = run_multibubble_2T_plasma(
        N=8,
        t_end=1.0e-3,     # 1 ms
        dt=1.0e-7,        # 0.1 µs step
        stride=5,
        seed=2
    )
    print(f"Done. Total bursts (all bubbles): {out['total_bursts']}")
    print(f"Stored points: {len(out['t'])}")
    plot_results(out, title="Multi-bubble (2T plasma EOS): radius, energy, pressure, Te, camera flash")