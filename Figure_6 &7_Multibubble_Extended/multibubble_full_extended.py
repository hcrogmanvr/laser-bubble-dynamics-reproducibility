"""
multibubble_full_extended.py

FULL model implementing (in order, clearly marked in code):

✅ Step 4.4 — Membrane rupture model (pressure impulse ∫(P-P0)dt + peak pressure)
✅ Step 4.5 — Mass accretion & opacity feedback (M(t) increases absorption + emission → runaway)
✅ Step 4.6 — Replace logistic ionization with Saha EOS (density + Te dependent ionization)
✅ Step 5   — Spatial shock propagation (collapse emits propagating shocklets with time-of-flight)

Core already included:
- Multi-bubble coupling (pressure interaction)
- Two-temperature model (Te, Ti) with e–i coupling and cooling
- Collapse-triggered burst and camera response

Dependencies: numpy, matplotlib (no scipy)

---------------------------------------------------------------------------
How to use
---------------------------------------------------------------------------
Run this file in Jupyter or python:

    %run multibubble_full_extended.py

Adjust the parameter block inside main() for your system (laser duty, N bubbles, etc.).
"""

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
import math
import heapq
import time


# -----------------------------
# Helpers
# -----------------------------
def sigmoid(x):
    # stable sigmoid
    x = np.clip(x, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-x))


# -----------------------------
# Models / Parameters
# -----------------------------
@dataclass
class Medium:
    rho: float = 1000.0          # kg/m^3
    mu: float = 1e-3             # Pa*s
    sigma: float = 0.072         # N/m
    c_sound: float = 1500.0      # m/s
    P0: float = 1e5              # Pa


@dataclass
class Laser:
    period: float = 80e-6        # s
    on_frac: float = 0.75
    Pin_eff: float = 2.0e-3      # W per bubble (effective local absorbed scale)
    abs_base: float = 0.35       # baseline absorption fraction
    abs_R_scale: float = 2e-6    # m
    skin_depth: float = 1e-7     # m

    def on(self, t):
        ph = (t % self.period) / self.period
        return 1.0 if ph < self.on_frac else 0.0

    def absorbed_power(self, t, R, abs_mult):
        """
        Absorption proxy:
        - increases with R (intercept)
        - limited by skin-depth at very small R
        - multiplied by abs_mult (mass accretion feedback)
        """
        L = self.on(t)
        frac_R = 1.0 - np.exp(-R / self.abs_R_scale)
        frac_skin = 1.0 - np.exp(-R / self.skin_depth)
        frac = self.abs_base * frac_R * (0.3 + 0.7 * frac_skin)
        frac = np.clip(frac * abs_mult, 0.0, 0.95)
        return L * self.Pin_eff * frac, frac, L


@dataclass
class TwoTemp:
    Tamb: float = 300.0
    Ce: float = 3e-12            # J/K (effective)
    Ci: float = 2e-11            # J/K
    tau_ei: float = 2.0e-6       # s
    tau_cool_e: float = 20e-6    # s
    tau_cool_i: float = 80e-6    # s
    fe_abs_to_e: float = 0.9

    def step(self, Te, Ti, Pabs, dt):
        # electron
        dTe_abs = (self.fe_abs_to_e * Pabs / self.Ce) * dt
        dTe_ei = (-(Te - Ti) / self.tau_ei) * dt
        dTe_cool = (-(Te - self.Tamb) / self.tau_cool_e) * dt
        Te_new = Te + dTe_abs + dTe_ei + dTe_cool

        # ions/vapor
        dTi_abs = ((1.0 - self.fe_abs_to_e) * Pabs / self.Ci) * dt
        dTi_ei = ((Te - Ti) / self.tau_ei) * (self.Ce / self.Ci) * dt
        dTi_cool = (-(Ti - self.Tamb) / self.tau_cool_i) * dt
        Ti_new = Ti + dTi_abs + dTi_ei + dTi_cool

        Te_new = np.maximum(Te_new, 50.0)
        Ti_new = np.maximum(Ti_new, 50.0)
        return Te_new, Ti_new

    def energy_above_ambient(self, Te, Ti):
        E = self.Ce * (Te - self.Tamb) + self.Ci * (Ti - self.Tamb)
        return np.maximum(E, 0.0)


@dataclass
class BubbleDynamics:
    R_min: float = 5e-8
    R_max: float = 3e-5
    dR_max: float = 250.0        # m/s
    accel_cap: float = 2e10      # m/s^2

    def step(self, R, dR, P_in, medium: Medium, dt):
        """
        Stable RP-like update:
        d2R = (P_in - P_out)/rho - 1.5(dR^2)/R
        with visc + surface in P_out and strong regularization.
        """
        R_safe = np.maximum(R, self.R_min)
        P_out = (medium.P0
                 + 4.0 * medium.mu * dR / R_safe
                 + 2.0 * medium.sigma / R_safe)

        rhs = (P_in - P_out) / medium.rho - 1.5 * (dR * dR) / R_safe
        rhs = np.clip(rhs, -self.accel_cap, self.accel_cap)

        dR_new = np.clip(dR + rhs * dt, -self.dR_max, self.dR_max)
        R_new = np.clip(R + dR_new * dt, self.R_min, self.R_max)
        return R_new, dR_new


# ===== STEP 4.6: SAHA IONIZATION EOS =====
@dataclass
class SahaEOS:
    """
    Saha ionization (single ionization) with density n and electron temperature Te.
    Returns ionization fraction x in [0,1] and pressure P.

    Equation:
      x^2/(1-x) = A
      A = (2/n) * (2π m_e kB Te / h^2)^(3/2) * exp(-Eion/(kB Te))
    Solve:
      x^2 + A x - A = 0
      x = (-A + sqrt(A^2 + 4A)) / 2
    """
    Eion_eV: float = 12.6        # eV (proxy for water vapor fragments)
    P_cap: float = 2.0e9         # Pa safety cap

    def __post_init__(self):
        self.kB = 1.380649e-23
        self.h = 6.62607015e-34
        self.me = 9.10938356e-31
        self.Eion = self.Eion_eV * 1.602176634e-19  # J

    def ion_frac(self, n, Te):
        Te = np.maximum(Te, 50.0)
        n = np.maximum(n, 1e10)

        pref = (2.0 / n) * ((2.0 * np.pi * self.me * self.kB * Te) / (self.h * self.h))**1.5
        A = pref * np.exp(-self.Eion / (self.kB * Te))
        A = np.clip(A, 0.0, 1e12)

        x = (-A + np.sqrt(A*A + 4.0*A)) / 2.0
        return np.clip(x, 0.0, 1.0)

    def pressure(self, n_total, Ti, Te):
        """
        Pressure with neutral+ion at Ti and electrons at Te:
          P = n_total kB Ti + x n_total kB Te
        """
        x = self.ion_frac(n_total, Te)
        P = n_total * self.kB * Ti + (x * n_total) * self.kB * Te
        return np.minimum(P, self.P_cap), x


# ===== STEP 4.5: MASS ACCRETION & OPACITY FEEDBACK =====
@dataclass
class Accretion:
    k_acc: float = 4e5           # 1/s (tune)
    P_thresh: float = 2.0e6      # Pa threshold to start accretion
    dP: float = 5.0e5            # Pa smoothness
    M_cap: float = 10.0          # cap for stability
    alpha_abs: float = 0.25      # absorption multiplier: abs_mult = 1 + alpha_abs*M
    beta_emit: float = 0.35      # emission multiplier: emit_mult = 1 + beta_emit*M

    def step(self, M, P_drive, dt):
        drive = sigmoid((P_drive - self.P_thresh) / self.dP)
        M_new = M + (self.k_acc * drive) * dt
        return np.clip(M_new, 0.0, self.M_cap)

    def abs_mult(self, M):
        return 1.0 + self.alpha_abs * M

    def emit_mult(self, M):
        return 1.0 + self.beta_emit * M


# ===== STEP 4.4: MEMBRANE RUPTURE (PRESSURE IMPULSE) =====
@dataclass
class Membrane:
    Pyield: float = 6.0e6        # Pa (peak-yield threshold)
    Jcrit: float = 2.0e4         # Pa*s impulse threshold
    use_peak: bool = True
    use_impulse: bool = True

    def update(self, J, P_total, P0, dt):
        # accumulate only overpressure (avoid counting ambient)
        J_new = J + np.maximum(P_total - P0, 0.0) * dt
        return J_new

    def ruptured(self, J, P_peak):
        cond = False
        if self.use_impulse:
            cond = cond or (J >= self.Jcrit)
        if self.use_peak:
            cond = cond or (P_peak >= self.Pyield)
        return cond


# ===== STEP 5: SPATIAL SHOCK PROPAGATION =====
@dataclass
class Shock:
    """
    Each collapse emits a shocklet:
      P(r,t) = (A / r) * exp(-(t - t_arr)/tau_s) for t >= t_arr
      t_arr = t_c + r/c

    We implement time-of-flight using a heap of (arrival_time, target_index, amplitude).
    """
    tau_s: float = 3.0e-6        # s decay
    A0: float = 3.0e3            # Pa*m per Joule of burst energy (tune)
    P_cap: float = 8.0e8         # Pa cap
    r_floor: float = 0.8e-6      # m minimum distance in 1/r

    def amplitude(self, Eburst, r):
        r_eff = max(r, self.r_floor)
        return np.clip(self.A0 * Eburst / r_eff, -self.P_cap, self.P_cap)


# Burst & camera (ties to Results narrative)
@dataclass
class BurstAndCamera:
    R_collapse: float = 1.5e-7
    burst_to_rad: float = 0.08
    wavelength_nm: float = 550.0
    quench_frac_E: float = 0.75
    Te_kick: float = 4e4
    cam_tau: float = 20e-6
    cam_gain: float = 1e-3
    cam_sat: float = 2.0

    def photon_energy(self):
        h = 6.62607015e-34
        c = 299792458.0
        lam = self.wavelength_nm * 1e-9
        return h * c / lam

    def cam_step(self, Iprev, photon_rate, dt):
        # dI/dt = -I/tau + gain*photon_rate
        I = Iprev + (-(Iprev)/self.cam_tau + self.cam_gain * photon_rate) * dt
        return float(np.clip(I, 0.0, self.cam_sat))


# -----------------------------
# Coupling model (pressure coupling between bubbles)
# -----------------------------
@dataclass
class Coupling:
    k_couple: float = 0.25
    length: float = 2.5e-6
    max_couple: float = 5e8

    def build_weights(self, pos):
        N = pos.shape[0]
        d = pos[:, None, :] - pos[None, :, :]
        r = np.linalg.norm(d, axis=2)
        r = np.maximum(r, 1e-12)
        W = np.exp(-r / self.length) / r
        np.fill_diagonal(W, 0.0)
        W = W / (W.sum(axis=1, keepdims=True) + 1e-30)
        return W

    def coupled_pressure(self, W, P_drive):
        dP = P_drive - P_drive.mean()
        P_c = self.k_couple * (W @ dP)
        return np.clip(P_c, -self.max_couple, self.max_couple)


# -----------------------------
# Main simulation
# -----------------------------
def simulate(
    N=8,
    t_end=1e-3,
    dt=1e-7,
    stride=5,
    seed=2,
    cluster_sigma=2.0e-6
):
    rng = np.random.default_rng(seed)

    # Instantiate models
    medium = Medium()
    laser = Laser()
    therm = TwoTemp()
    dyn = BubbleDynamics()
    eos = SahaEOS(Eion_eV=12.6, P_cap=2.0e9)
    acc = Accretion()
    mem = Membrane(Pyield=6.0e6, Jcrit=2.0e4, use_peak=True, use_impulse=True)
    shock = Shock(tau_s=3.0e-6, A0=3.0e3, P_cap=8.0e8)
    bc = BurstAndCamera()
    couple = Coupling()

    # Positions and coupling weights
    pos = rng.normal(scale=cluster_sigma, size=(N, 3))
    W = couple.build_weights(pos)

    # Initial state
    R = np.clip(1.8e-7 + rng.normal(scale=0.2e-7, size=N), dyn.R_min, dyn.R_max)
    dR = np.zeros(N)
    Te = np.full(N, therm.Tamb)
    Ti = np.full(N, therm.Tamb)

    # Fixed particle count for density scaling (same as earlier)
    R_ref = 2.0e-7
    V_ref = (4.0/3.0) * np.pi * (R_ref**3)
    n_ref = 1.0e25
    N_particles = n_ref * V_ref  # per bubble

    # STEP 4.5 state
    M = np.zeros(N)              # accreted mass proxy

    # STEP 4.4 state
    Jmem = np.zeros(N)           # pressure impulse
    rupt = np.zeros(N, dtype=bool)
    Ppeak = np.zeros(N)

    # STEP 5 shock state: heap of arrivals + running shock pressure per bubble
    arrivals = []  # heap items: (t_arrival, j_target, amp_at_target)
    Pshock = np.zeros(N)

    # Camera
    hv = bc.photon_energy()
    Icam = 0.0

    # History (downsample)
    steps = int(np.ceil(t_end / dt))
    store_n = steps // stride + 1
    t_hist = np.zeros(store_n)
    Rm_hist = np.zeros(store_n)
    Em_hist = np.zeros(store_n)
    Pm_hist = np.zeros(store_n)
    Tem_hist = np.zeros(store_n)
    I_hist = np.zeros(store_n)
    laser_hist = np.zeros(store_n)
    burst_hist = np.zeros(store_n)
    xion_hist = np.zeros(store_n)
    M_hist = np.zeros(store_n)
    rupt_frac_hist = np.zeros(store_n)

    total_bursts = 0
    si = 0

    # For speed: precompute pairwise distances for shocks
    # (used only when collapse occurs)
    rij = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)

    t0 = time.time()

    for k in range(steps + 1):
        t = k * dt

        # ---- STEP 5: apply shock arrivals + decay existing shock pressure
        # Exponential decay of existing shock field
        if shock.tau_s > 0:
            Pshock *= math.exp(-dt / shock.tau_s)

        # Apply arrivals due now
        while arrivals and arrivals[0][0] <= t:
            _, j, amp = heapq.heappop(arrivals)
            Pshock[j] = np.clip(Pshock[j] + amp, -shock.P_cap, shock.P_cap)

        # ---- STEP 4.5: absorption multiplier from accreted mass
        abs_mult = acc.abs_mult(M)

        # Laser absorbed power
        Pabs, abs_frac, Lflag = laser.absorbed_power(t, R, abs_mult)

        # Two-temperature update
        Te, Ti = therm.step(Te, Ti, Pabs, dt)
        E = therm.energy_above_ambient(Te, Ti)

        # Density from fixed particle count
        V = (4.0/3.0) * np.pi * np.maximum(R, dyn.R_min)**3
        n_total = N_particles / V

        # ---- STEP 4.6: SAHA EOS
        P_drive, xion = eos.pressure(n_total, Ti, Te)

        # ---- Coupling pressure + shock pressure
        P_couple = couple.coupled_pressure(W, P_drive)
        P_total = P_drive + P_couple + Pshock

        # Track peaks for membrane model
        Ppeak = np.maximum(Ppeak, P_total)

        # ---- Bubble dynamics
        R_prev = R.copy()
        dR_prev = dR.copy()
        R, dR = dyn.step(R, dR, P_total, medium, dt)

        # ---- Collapse detection + burst
        crossed = (R_prev > bc.R_collapse) & (R <= bc.R_collapse) & (dR_prev < 0.0)

        burst_energy = np.zeros(N)
        if np.any(crossed):
            total_bursts += int(crossed.sum())

            # Burst energy: fraction of stored internal energy + kinetic proxy
            Vc = (4.0/3.0) * np.pi * np.maximum(R_prev[crossed], dyn.R_min)**3
            K = 0.5 * medium.rho * Vc * (dR_prev[crossed]**2)
            Eb = 0.35 * E[crossed] + 0.25 * K
            burst_energy[crossed] = Eb

            # Quench energy after burst
            E[crossed] *= (1.0 - bc.quench_frac_E)
            # Push Te transiently (plasma spike proxy)
            Te[crossed] = np.minimum(Te[crossed] + bc.Te_kick, 2.0e5)

            # ---- STEP 5: emit shocks to neighbors with time-of-flight
            # For each collapsing bubble i, schedule arrivals to all j != i
            idx_i = np.where(crossed)[0]
            for i in idx_i:
                Ei = float(burst_energy[i])
                if Ei <= 0:
                    continue
                for j in range(N):
                    if j == i:
                        continue
                    r = float(rij[i, j])
                    t_arr = t + r / medium.c_sound
                    amp = shock.amplitude(Ei, r)  # already includes 1/r
                    heapq.heappush(arrivals, (t_arr, j, amp))

            # Reset peaks after rupture check windows if desired (optional); keep here simple.

        # ---- STEP 4.5: mass accretion update (driven by pressure)
        # Use P_drive (local) or P_total (including shocks). P_total makes accretion more explosive.
        M = acc.step(M, P_total, dt)

        # ---- STEP 4.4: membrane impulse + rupture logic
        # Once ruptured, keep it ruptured (absorbing boundary)
        for i in range(N):
            if rupt[i]:
                continue
            Jmem[i] = mem.update(Jmem[i], P_total[i], medium.P0, dt)
            if mem.ruptured(Jmem[i], Ppeak[i]):
                rupt[i] = True

        # ---- Camera response (collapse-triggered photons)
        photon_rate = (bc.burst_to_rad * burst_energy.sum()) / (hv * dt + 1e-30)
        # emission boosted by accreted mass (opacity feedback)
        photon_rate *= float(np.mean(acc.emit_mult(M)))
        Icam = bc.cam_step(Icam, photon_rate, dt)

        # ---- Store
        if (k % stride) == 0:
            t_hist[si] = t
            Rm_hist[si] = R.mean()
            Em_hist[si] = E.mean()
            Pm_hist[si] = (P_total.mean() - medium.P0) * 1e-6  # MPa overpressure
            Tem_hist[si] = Te.mean()
            I_hist[si] = Icam
            laser_hist[si] = Lflag
            burst_hist[si] = 1.0 if np.any(crossed) else 0.0
            xion_hist[si] = xion.mean()
            M_hist[si] = M.mean()
            rupt_frac_hist[si] = rupt.mean()
            si += 1

    runtime = time.time() - t0

    out = {
        "t": t_hist[:si],
        "R_mean": Rm_hist[:si],
        "E_mean": Em_hist[:si],
        "P_over_MPa": Pm_hist[:si],
        "Te_mean": Tem_hist[:si],
        "I_cam": I_hist[:si],
        "laser": laser_hist[:si],
        "burst_flag": burst_hist[:si],
        "xion_mean": xion_hist[:si],
        "M_mean": M_hist[:si],
        "rupt_frac": rupt_frac_hist[:si],
        "total_bursts": total_bursts,
        "runtime_s": runtime,
        "final_ruptured": rupt,
        "final_Jmem": Jmem,
        "final_Ppeak": Ppeak,
    }
    return out


def plot_all(out, title="Full model: rupture + accretion + Saha EOS + shock propagation"):
    t = out["t"]
    t_ms = t * 1e3

    fig, ax = plt.subplots(8, 1, figsize=(12, 13), sharex=True)
    ax[0].plot(t_ms, out["R_mean"] * 1e6)
    ax[0].set_ylabel("⟨R⟩ (µm)")
    ax[0].set_title(title)

    ax[1].plot(t_ms, out["laser"])
    ax[1].set_ylabel("Laser")

    ax[2].plot(t_ms, out["E_mean"])
    ax[2].set_ylabel("⟨E⟩ (J)")

    ax[3].plot(t_ms, out["P_over_MPa"])
    ax[3].set_ylabel("⟨P-P0⟩ (MPa)")

    ax[4].plot(t_ms, out["Te_mean"])
    ax[4].set_ylabel("⟨Te⟩ (K)")

    ax[5].plot(t_ms, out["xion_mean"])
    ax[5].set_ylabel("⟨x_ion⟩")

    ax[6].plot(t_ms, out["M_mean"])
    ax[6].set_ylabel("⟨M⟩")

    ax[7].plot(t_ms, out["I_cam"], label="Camera")
    mask = out["burst_flag"] > 0
    if np.any(mask):
        ax[7].vlines(t_ms[mask], 0, max(out["I_cam"].max(), 1e-12), alpha=0.2, label="Bursts")
    ax[7].set_ylabel("I_cam (arb)")
    ax[7].set_xlabel("Time (ms)")
    ax[7].legend(loc="upper right")

    plt.tight_layout()
    plt.show()

    # Extra plot: rupture fraction
    plt.figure(figsize=(12, 3))
    plt.plot(t_ms, out["rupt_frac"])
    plt.ylabel("Ruptured fraction")
    plt.xlabel("Time (ms)")
    plt.title("Membrane rupture accumulation (Step 4.4)")
    plt.ylim(-0.05, 1.05)
    plt.grid(True)
    plt.show()


def main():
    # -----------------------------
    # USER-TUNABLE BLOCK
    # -----------------------------
    out = simulate(
        N=8,
        t_end=1.0e-3,      # 1 ms
        dt=1.0e-7,         # 0.1 µs
        stride=5,
        seed=2,
        cluster_sigma=2.0e-6
    )

    print(f"Runtime: {out['runtime_s']:.3f} s")
    print(f"Total bursts: {out['total_bursts']}")
    print(f"Final ruptured fraction: {out['final_ruptured'].mean():.2f}")

    plot_all(out)


if __name__ == "__main__":
    main()