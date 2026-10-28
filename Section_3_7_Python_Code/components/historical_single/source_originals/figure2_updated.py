"""Original Figure 2 dynamics plus manuscript optical diagnostics.

Run: python figure2_updated.py
Optics: python figure2_updated.py --r-sat-m VALUE --n VALUE --chi-csv FILE
chi CSV columns: time_s,chi (must cover the simulation; no extrapolation).
An explicit --neutral-limit uses chi=0, labelled as a limiting assumption.
No optical parameters or ionization calibration are invented by default.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def simulate():
    # Original parameters and update ordering preserved.
    dt = 0.01
    t = np.arange(0, 300.0, dt)
    N = len(t)
    laser_on, laser_off = 60.0, 30.0
    laser = ((t % (laser_on + laser_off)) < laser_on).astype(float)
    R0, Rmin, Rmax, Rc = 2e-6, 0.5e-6, 30e-6, 1.5e-6
    k_restore, gamma_damp, m_eff, v_max = 5e7, 3e4, 1e-12, 200.0
    E_tau, E_abs_scale, E_to_P, P_max, P_ext = 15.0, 8.0, 2e8, 1.5e9, 1e5
    burst_fraction, burst_cooldown = 0.4, 0.5
    eta_rad, tau_flash, camera_exposure = 1e-3, 2e-3, 5e-3
    R, v, E, P, Mach, burst_flag, E_burst = [np.zeros(N) for _ in range(7)]
    R[0] = R0
    last_burst_time, armed = -np.inf, True
    for i in range(1, N):
        dE = (laser[i] * E_abs_scale - E[i-1] / E_tau) * dt
        E[i] = max(E[i-1] + dE, 0.0)
        P[i] = min(E_to_P * E[i], P_max)
        F_restore = -k_restore * (R[i-1] - R0)
        F_pressure = 4 * np.pi * R[i-1]**2 * (P[i] - P_ext)
        F_damp = -gamma_damp * v[i-1]
        a = (F_restore + F_pressure + F_damp) / m_eff
        v[i] = np.clip(v[i-1] + a * dt, -v_max, v_max)
        R[i] = np.clip(R[i-1] + v[i] * dt, Rmin, Rmax)
        Mach[i] = abs(v[i]) / 1500.0
        collapsing = (R[i] < Rc) and (v[i] < 0)
        if collapsing and armed and (t[i] - last_burst_time > burst_cooldown):
            burst_flag[i] = 1.0
            Eb = burst_fraction * E[i]
            E_burst[i] = Eb
            E[i] -= Eb
            last_burst_time, armed = t[i], False
        if R[i] > 1.2 * Rc:
            armed = True
    # Original sampled flash and camera calculation retained for reproducibility.
    I_flash = np.zeros(N)
    for k in np.where(E_burst > 0)[0]:
        E_radiated = eta_rad * E_burst[k]  # renamed to avoid confusion with A(t)
        mask = t >= t[k]
        I_flash[mask] += (E_radiated / tau_flash) * np.exp(-(t[mask]-t[k])/tau_flash)
    w = max(1, int(camera_exposure / dt))
    I_cam = np.convolve(I_flash, np.ones(w)/w, mode='same')
    return dict(t=t, R=R, v=v, E=E, P=P, Mach=Mach, laser=laser,
                burst_flag=burst_flag, E_burst=E_burst, I_cam=I_cam)


def off_phase_minima(s):
    """Minimum in each COMPLETE OFF window; not a changing radius bound.

    Discrete OFF windows: [60,90), [150,180), [240,270).
    Record both actual minimum and final sampled OFF-phase radius.
    """
    rows = []
    for start in np.arange(60.0, s['t'][-1], 90.0):
        end = start + 30.0
        if end > s['t'][-1] + 0.01 + 1e-9:
            continue
        ix = np.flatnonzero((s['t'] >= start) & (s['t'] < end))
        k = ix[np.argmin(s['R'][ix])]
        rows.append([start, end, s['t'][k], s['R'][k], s['R'][ix[-1]]])
    return np.asarray(rows, dtype=float).reshape(-1, 5)


def optical_observables(s, r_sat, n, chi):
    """Manuscript Eqs. (3)-(4); diagnostic only, no feedback into dynamics."""
    if not (np.isfinite(r_sat) and r_sat > 0 and np.isfinite(n) and n > 0):
        raise ValueError('R_sat and n must be finite and positive.')
    chi = np.asarray(chi, dtype=float)
    if chi.shape != s['t'].shape or not np.all(np.isfinite(chi)) or np.any((chi < 0) | (chi > 1)):
        raise ValueError('chi must match the time grid and lie in [0,1].')
    with np.errstate(over='ignore'):
        T = np.exp(-(s['R']/r_sat)**n)
    return T, (1-T)*(1-chi)


def save_figure(fig, out, name):
    fig.savefig(out / (name+'.pdf'), bbox_inches='tight')
    fig.savefig(out / (name+'.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=Path('figure2_outputs'))
    p.add_argument('--r-sat-m', type=float)
    p.add_argument('--n', type=float)
    source = p.add_mutually_exclusive_group()
    source.add_argument('--chi-csv', type=Path)
    source.add_argument('--neutral-limit', action='store_true')
    args = p.parse_args()
    supplied = [args.r_sat_m is not None, args.n is not None, bool(args.chi_csv or args.neutral_limit)]
    if any(supplied) and not all(supplied):
        p.error('Optics requires --r-sat-m, --n, and either --chi-csv or --neutral-limit.')
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    s = simulate()
    t = s['t']
    minima = off_phase_minima(s)
    fig, axs = plt.subplots(9, 1, figsize=(10, 18), sharex=True)
    series = [('R',1e6,'Bubble radius (µm)'),('laser',1,'Laser gate'),
              ('E',1,'Internal E (proxy)'),('P',1e-6,'Pressure (MPa)'),
              ('v',1,'Velocity (m/s)'),('Mach',1,'Mach number'),
              ('burst_flag',1,'Burst flag'),('E_burst',1,'Burst E (proxy)'),
              ('I_cam',1,'Optical signal (proxy)')]
    for j,(key,scale,label) in enumerate(series):
        axs[j].plot(t,s[key]*scale,linewidth=0.7)
        axs[j].set_ylabel(label)
        axs[j].text(0.01,0.8,f'({chr(97+j)})',transform=axs[j].transAxes)
    axs[5].axhline(1,ls='--',color='gray')
    events = t[s['E_burst']>0]
    axs[8].vlines(events,0,s['I_cam'].max(),alpha=0.4,label='Collapse triggers')
    axs[8].legend(loc='upper right')
    axs[8].set_xlabel('Time (s)')
    fig.tight_layout()
    save_figure(fig,out,'Figure_2_original_dynamics')
    np.savetxt(out/'off_phase_minima.csv',minima,delimiter=',',comments='',
               header='off_start_s,off_end_s,min_time_s,min_radius_m,end_radius_m')
    np.savetxt(out/'original_states.csv',np.column_stack([s[k] for k in s]),
               delimiter=',',header=','.join(s),comments='')
    # Always provide the computable minimum-radius diagnostic.
    fig,ax=plt.subplots(figsize=(8,3))
    ax.plot(minima[:,1],minima[:,3]*1e6,'o-',label='Minimum during each OFF phase')
    ax.axhline(0.5,ls='--',color='gray',label='Fixed numerical lower bound')
    ax.set(xlabel='OFF-phase end time (s)',ylabel='Minimum radius (µm)',ylim=(0,1))
    ax.legend()
    fig.tight_layout()
    save_figure(fig,out,'Off_phase_minimum_radius')
    if all(supplied):
        if args.chi_csv:
            data=np.genfromtxt(args.chi_csv,delimiter=',',names=True)
            tt=np.atleast_1d(data['time_s']); cc=np.atleast_1d(data['chi'])
            if len(tt)<2 or not np.all(np.isfinite(tt)) or np.any(np.diff(tt)<=0):
                raise ValueError('chi CSV needs strictly increasing finite times.')
            if not np.all(np.isfinite(cc)) or np.any((cc<0)|(cc>1)):
                raise ValueError('chi CSV values must be finite and in [0,1].')
            if tt[0]>t[0] or tt[-1]<t[-1]:
                raise ValueError('chi CSV must cover the complete simulation time grid.')
            chi=np.interp(t,tt,cc)
            label='Supplied ionization history; diagnostic postprocessing'
        else:
            chi=np.zeros_like(t)
            label='Neutral limit: chi = 0 (assumed); diagnostic postprocessing'
        T,A=optical_observables(s,args.r_sat_m,args.n,chi)
        fig,axs=plt.subplots(3,1,figsize=(10,7),sharex=True)
        axs[0].plot(t,T,lw=0.7); axs[0].set_ylabel('Transmission T(t)')
        axs[1].plot(t,A,lw=0.7); axs[1].set_ylabel('Absorbed fraction A(t)')
        for ax in axs[:2]: ax.set_ylim(-0.02,1.02)
        axs[2].plot(minima[:,1],minima[:,3]*1e6,'o-')
        axs[2].set(ylabel='OFF-phase minimum (µm)',xlabel='Time (s)')
        fig.suptitle(label,fontsize=10)
        fig.tight_layout()
        save_figure(fig,out,'Companion_optical_and_minimum_radius')
        np.savetxt(out/'optical_observables.csv',np.column_stack((t,T,A,chi)),
                   delimiter=',',header='time_s,T,A,chi',comments='')
    diagnostics=dict(burst_count=int(len(events)),first_burst_times_s=events[:10].tolist(),
        emitted_energy_proxy=float((1e-3*s['E_burst']).sum()),
        lower_bound_fraction=float(np.mean(s['R']==0.5e-6)),
        upper_bound_fraction=float(np.mean(s['R']==30e-6)),
        off_phase_minima_um=(minima[:,3]*1e6).tolist(),
        camera_boxcar_samples=1,optics_generated=bool(all(supplied)),
        r_sat_m=args.r_sat_m,n=args.n,chi_source=str(args.chi_csv) if args.chi_csv else ('assumed neutral limit' if args.neutral_limit else None),
        limitations=['Original explicit mechanics and clipping preserved; numerical convergence unverified.',
        'dt=10 ms does not resolve 2-ms flash decay or 5-ms camera exposure.',
        'Energy and optical signal remain proxies; no joule or watt calibration supplied.',
        'Optics is postprocessing, not a coupled absorption-feedback model.',
        'Off-phase minima are measured from R(t), not imposed as a growth law.'])
    (out/'diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    print(json.dumps(diagnostics,indent=2))
    if not all(supplied):
        print('Companion optical figure pending: supply R_sat, n and chi history. See --help.')


if __name__ == '__main__':
    main()
