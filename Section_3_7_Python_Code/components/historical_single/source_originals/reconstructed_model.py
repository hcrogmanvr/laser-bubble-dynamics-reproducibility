from dataclasses import dataclass
import numpy as np

@dataclass
class Params:
    wavelength_nm: float = 1064.0
    power_W: float = 0.8
    on_s: float = 60.0
    off_s: float = 30.0
    R0_um: float = 0.5
    Rmax_um: float = 30.0
    rho_liquid: float = 1000.0
    P_ambient_Pa: float = 101325.0
    collapse_radius_um: float = 1.0
    collapse_velocity_m_s: float = -20.0
    heat_capacity_J_K: float = 7.08e-13
    tau_loss_s: float = 0.8
    burst_fraction: float = 0.25
    absorption_base: float = 0.55
    quench_strength: float = 0.65
    saturation_length_um: float = 8.0
    T0_K: float = 300.0
    ionization_threshold_K: float = 5000.0
    ionization_width_K: float = 700.0
    pressure_cap_Pa: float = 5.0e6
    gamma: float = 5.0/3.0
    coupling_strength: float = 0.08
    camera_tau_s: float = 0.01
    impulse_threshold_Pa_s: float = 0.02
    damage_scale_Pa_s: float = 0.05
    dt_s: float = 1e-4
    duration_s: float = 300.0


def laser_gate(t, on_s, off_s):
    period = on_s + off_s
    return 1.0 if (t % period) < on_s else 0.0


def transmission(R_um, L_um):
    # Phenomenological bounded transmission law used to represent optical darkening.
    return np.exp(-np.clip(R_um / max(L_um, 1e-12), 0, 20))


def ionization_fraction(T, threshold, width):
    z = np.clip((T - threshold) / max(width, 1e-12), -60, 60)
    return 1.0 / (1.0 + np.exp(-z))


def absorbed_fraction(R_um, T, p: Params):
    geometric = p.absorption_base * (1.0 - transmission(R_um, p.saturation_length_um))
    x = ionization_fraction(T, p.ionization_threshold_K, p.ionization_width_K)
    return np.clip(geometric * (1.0 - p.quench_strength * x), 0.0, 1.0)


def pressure_from_energy(E, R_m, p: Params):
    V = 4.0 * np.pi * max(R_m, 1e-9)**3 / 3.0
    raw = (p.gamma - 1.0) * E / V + p.P_ambient_Pa
    return min(raw, p.pressure_cap_Pa)


def temperature_proxy(E, p: Params):
    return p.T0_K + E / max(p.heat_capacity_J_K, 1e-30)


def simulate_single(p: Params):
    n = int(p.duration_s / p.dt_s) + 1
    t = np.linspace(0.0, p.duration_s, n)
    R = np.zeros(n); V = np.zeros(n); E = np.zeros(n)
    T = np.zeros(n); P = np.zeros(n); gate = np.zeros(n); trans = np.zeros(n)
    burst = np.zeros(n); flash = np.zeros(n); damage = np.zeros(n)
    R[0] = p.R0_um * 1e-6
    T[0] = p.T0_K
    P[0] = p.P_ambient_Pa
    camera_state = 0.0
    cumulative_impulse = 0.0

    for i in range(n-1):
        dt = p.dt_s
        gate[i] = laser_gate(t[i], p.on_s, p.off_s)
        T[i] = temperature_proxy(E[i], p)
        af = absorbed_fraction(R[i]*1e6, T[i], p)
        Pin = p.power_W * gate[i]
        dE = Pin * af - E[i] / p.tau_loss_s
        P[i] = pressure_from_energy(E[i], R[i], p)

        # Reduced Rayleigh-Plesset inertial balance: R*Rddot + 3/2 Rdot^2 = (Pb-Pinf)/rho
        R_safe = max(R[i], 1e-8)
        accel = ((P[i] - p.P_ambient_Pa) / p.rho_liquid - 1.5 * V[i]**2) / R_safe
        # Weak restoring term limits expansion to the intended reduced-order regime.
        accel -= 20.0 * max(R[i] - p.Rmax_um*1e-6, 0.0)
        accel = float(np.clip(accel, -5.0e5, 5.0e5))

        V[i+1] = float(np.clip(V[i] + accel * dt, -100.0, 100.0))
        R[i+1] = float(np.clip(R[i] + V[i+1] * dt, 0.2e-6, p.Rmax_um*1e-6))
        E[i+1] = max(0.0, E[i] + dE * dt)

        collapse = (R[i+1] <= p.collapse_radius_um*1e-6 and V[i+1] < p.collapse_velocity_m_s)
        if collapse:
            released = p.burst_fraction * E[i+1]
            burst[i+1] = released
            E[i+1] *= (1.0 - p.burst_fraction)
            V[i+1] *= -0.25
            R[i+1] = max(R[i+1], p.R0_um*1e-6)
            cumulative_impulse += max(P[i] - p.P_ambient_Pa, 0.0) * dt

        camera_state += dt * (-camera_state / p.camera_tau_s) + burst[i+1]
        flash[i+1] = camera_state
        trans[i] = transmission(R[i]*1e6, p.saturation_length_um)
        damage[i] = np.clip((cumulative_impulse - p.impulse_threshold_Pa_s) / p.damage_scale_Pa_s, 0, 1)

    gate[-1] = laser_gate(t[-1], p.on_s, p.off_s)
    T[-1] = temperature_proxy(E[-1], p)
    P[-1] = pressure_from_energy(E[-1], R[-1], p)
    trans[-1] = transmission(R[-1]*1e6, p.saturation_length_um)
    damage[-1] = damage[-2]
    return dict(t=t, R_um=R*1e6, velocity_m_s=V, energy_J=E, temperature_K=T,
                pressure_Pa=P, gate=gate, transmission=trans, burst_J=burst,
                camera_flash=flash, membrane_damage=damage)


def simulate_multi(p: Params, n_bubbles=5, seed=7):
    rng = np.random.default_rng(seed)
    runs = []
    for j in range(n_bubbles):
        pj = Params(**p.__dict__)
        pj.R0_um *= float(rng.uniform(0.9, 1.1))
        pj.absorption_base *= float(rng.uniform(0.9, 1.1))
        runs.append(simulate_single(pj))
    # Mean-field pressure coupling proxy.
    pressures = np.stack([r['pressure_Pa'] for r in runs])
    mean_excess = np.mean(np.maximum(pressures - p.P_ambient_Pa, 0.0), axis=0)
    coupled_pressures = pressures + p.coupling_strength * mean_excess[None, :]
    for j, r in enumerate(runs):
        r['pressure_coupled_Pa'] = coupled_pressures[j]
    return runs
