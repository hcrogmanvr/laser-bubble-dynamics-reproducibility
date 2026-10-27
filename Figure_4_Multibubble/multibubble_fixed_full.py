# multibubble_fixed_full.py
# Stable, event-based multi-bubble coupling with camera flash output.
# Fixes applied:
#   (1) Regularize 1/R and clip acceleration
#   (2) Enforce subsonic |dR| <= v_frac_c * c_sound
#   (3) Collapse quenches energy AND pressure
#   (4) Camera flash is event-based: burst energy -> intensity, with decay + saturation

import numpy as np
import matplotlib.pyplot as plt
import math
import time

def laser_schedule(t, on=60e-6, off=20e-6):
    period = on + off
    return 1.0 if (t % period) < on else 0.0

class Medium:
    def __init__(self, rho=1000.0, mu=1e-3, sigma=0.072, c_sound=1500.0, P_inf=1e5):
        self.rho = rho
        self.mu = mu
        self.sigma = sigma
        self.c_sound = c_sound
        self.P_inf = P_inf

class Camera:
    def __init__(self, tau=8e-6, gain=2e-3, saturate=1.0):
        self.tau = tau
        self.gain = gain
        self.saturate = saturate

class BubbleState:
    def __init__(self, R0, dR0=0.0):
        self.R = float(R0)
        self.dR = float(dR0)
        self.E = 0.0
        self.P = 0.0
        self.Te = 300.0

def simulate_multibubble(
    n_bubbles=8,
    t_end=1e-3,
    dt=1e-7,
    seed=4,
    medium=None,
    cam=None,
    R0=0.7e-6,
    R_min=5e-8,
    R_max=50e-6,
    Pin=1.0,
    abs_frac=0.25,
    E_tau=80e-6,
    E_to_P=3.5e8,
    P_cap=2.5e9,
    R_collapse=1.8e-7,
    burst_frac=0.7,
    rebound_factor=0.6,
    k_couple=0.12,
    r0=8e-6,
    a_clip=2e14,
    v_frac_c=0.8,
):
    if medium is None:
        medium = Medium()
    if cam is None:
        cam = Camera()

    rng = np.random.default_rng(seed)

    # Positions with min separation (prevents singular coupling)
    pos = []
    min_sep = 2.5e-6
    for _ in range(n_bubbles):
        for _try in range(5000):
            p = rng.normal(scale=r0, size=3)
            if all(np.linalg.norm(p - q) > min_sep for q in pos):
                pos.append(p)
                break
        else:
            pos.append(rng.normal(scale=r0, size=3))
    pos = np.array(pos)

    # Precompute inv(r^3) coupling matrix
    rij = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv_r3 = np.where(rij > 0, 1.0 / (rij**3), 0.0)
    np.fill_diagonal(inv_r3, 0.0)
    inv_r3 = np.minimum(inv_r3, 1.0 / (min_sep**3))

    bubbles = [BubbleState(R0 * (1 + 0.08 * rng.normal()), 0.0) for _ in range(n_bubbles)]

    n_steps = int(np.ceil(t_end / dt)) + 1
    t = np.linspace(0, t_end, n_steps)

    R_hist = np.zeros((n_steps, n_bubbles))
    E_hist = np.zeros((n_steps, n_bubbles))
    P_hist = np.zeros((n_steps, n_bubbles))
    Te_hist = np.zeros((n_steps, n_bubbles))
    laser_hist = np.zeros(n_steps)
    burst_hist = np.zeros((n_steps, n_bubbles))
    cam_hist = np.zeros(n_steps)

    def energy_to_Te(E):
        # proxy temperature (replace with real plasma EOS later)
        return 300.0 + 2500.0 * np.log1p(max(E, 0.0))

    def energy_to_pressure(E):
        # bounded pressure proxy (replace with real plasma EOS later)
        P = E_to_P * max(E, 0.0)
        return P_cap * (P / (P + P_cap))  # saturating cap

    v_max = v_frac_c * medium.c_sound
    last_R = np.array([b.R for b in bubbles], float)
    total_bursts = 0

    for k in range(n_steps):
        tk = t[k]
        laser_on = laser_schedule(tk)
        laser_hist[k] = laser_on

        # Energy update
        for b in bubbles:
            dE = (Pin * abs_frac * laser_on) - (b.E / E_tau)
            b.E = max(0.0, b.E + dE * dt)
            b.Te = energy_to_Te(b.E)
            b.P = energy_to_pressure(b.E)

        # Coupled pressure from neighbors (monopole-like scaling)
        P_vec = np.array([b.P for b in bubbles], float)
        vol_scale = np.array([(b.R**3) / (R0**3) for b in bubbles], float)
        P_couple = k_couple * (inv_r3 @ ((P_vec - medium.P_inf) * vol_scale)) * (r0**3)

        # Rayleigh–Plesset-like radius update (simplified)
        for i, b in enumerate(bubbles):
            R = b.R
            dR = b.dR

            R_eff = max(R, R_min)  # FIX (1): regularize 1/R
            visc = 4.0 * medium.mu * dR / R_eff
            tens = 2.0 * medium.sigma / R_eff

            P_inside = b.P + P_couple[i]
            rhs = (P_inside - medium.P_inf - tens - visc) / medium.rho

            d2R = (rhs - 1.5 * dR * dR) / R_eff
            d2R = float(np.clip(d2R, -a_clip, a_clip))  # FIX (1)

            dR_new = dR + d2R * dt
            dR_new = float(np.clip(dR_new, -v_max, v_max))  # FIX (2)

            R_new = float(np.clip(R + dR_new * dt, R_min, R_max))

            # Collapse detection + burst
            collapse = (last_R[i] > R_collapse) and (R_new <= R_collapse) and (dR_new < 0)
            if collapse:
                total_bursts += 1
                E_before = b.E
                Eburst = burst_frac * E_before
                burst_hist[k, i] = Eburst

                # FIX (3): collapse quenches energy AND pressure
                b.E = 0.2 * E_before
                b.P = 0.1 * b.P
                b.Te = energy_to_Te(b.E)

                # rebound
                dR_new = abs(dR_new) * rebound_factor
                R_new = max(R_new, R_collapse)

            b.R = R_new
            b.dR = dR_new
            last_R[i] = R_new

        # FIX (4): camera response to event bursts
        burst_power = burst_hist[k].sum() / max(dt, 1e-30)
        cam_prev = cam_hist[k-1] if k > 0 else 0.0
        cam_now = cam_prev * math.exp(-dt / cam.tau) + cam.gain * burst_power
        cam_hist[k] = min(cam_now, cam.saturate)

        # Record
        R_hist[k] = [b.R for b in bubbles]
        E_hist[k] = [b.E for b in bubbles]
        P_hist[k] = [b.P for b in bubbles]
        Te_hist[k] = [b.Te for b in bubbles]

    out = dict(t=t, R=R_hist, E=E_hist, P=P_hist, Te=Te_hist,
               laser=laser_hist, burst=burst_hist, cam=cam_hist, pos=pos)
    return out, total_bursts

if __name__ == "__main__":
    medium = Medium()
    cam = Camera(tau=8e-6, gain=2e-3, saturate=1.0)

    t0 = time.time()
    out, total_bursts = simulate_multibubble(medium=medium, cam=cam)
    print(f"Done. Total bursts (all bubbles): {total_bursts}")
    print(f"Runtime: {time.time()-t0:.3f} s")

    t = out["t"]
    R_um = out["R"] * 1e6
    E = out["E"]
    P_MPa = out["P"] / 1e6
    Te = out["Te"]
    laser = out["laser"]
    cam_sig = out["cam"]
    burst_events = (out["burst"].sum(axis=1) > 0)

    plt.figure(figsize=(12, 10))
    gs = plt.GridSpec(6, 1, hspace=0.25)

    ax1 = plt.subplot(gs[0])
    ax1.plot(t*1e3, R_um.mean(axis=1))
    ax1.set_ylabel("Mean R (µm)")
    ax1.set_title("Multi-bubble (fixed): radius, energy, pressure, and camera flash")

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.plot(t*1e3, laser)
    ax2.set_ylabel("Laser (0/1)")

    ax3 = plt.subplot(gs[2], sharex=ax1)
    ax3.plot(t*1e3, E.mean(axis=1))
    ax3.set_ylabel("Mean E (arb)")

    ax4 = plt.subplot(gs[3], sharex=ax1)
    ax4.plot(t*1e3, P_MPa.mean(axis=1))
    ax4.set_ylabel("Mean P (MPa)")

    ax5 = plt.subplot(gs[4], sharex=ax1)
    ax5.plot(t*1e3, Te.mean(axis=1))
    ax5.set_ylabel("Mean Te (K)")

    ax6 = plt.subplot(gs[5], sharex=ax1)
    ax6.plot(t*1e3, cam_sig, label="Camera flash")
    idx = np.where(burst_events)[0]
    if len(idx) > 0:
        ax6.vlines(t[idx]*1e3, 0, cam_sig.max(), alpha=0.2, label="Burst events")
    ax6.set_ylabel("I_cam (arb)")
    ax6.set_xlabel("Time (ms)")
    ax6.legend(loc="upper right")

    plt.show()
