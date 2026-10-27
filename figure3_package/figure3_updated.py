"""User's Figure 3 model plus explicitly UNCALIBRATED diagnostic extensions.
Run: python figure3_updated.py --out results
Dependencies: numpy, matplotlib. No experimental validation is implied.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def simulate(growing_floor=False, tau0=1.0, volume_gain=1.0):
    # Original constants, parameters and operation order are retained.
    kB, c_sound, rho, P_amb = 1.380649e-23, 1500.0, 1000.0, 1e5
    dt, t_end = 0.01, 300.0
    t = np.arange(0, t_end, dt); N = len(t)
    laser = ((t % 90.0) < 60.0).astype(float)
    R0, Rmin, Rmax, Rc = 2e-6, 0.4e-6, 30e-6, 1.2e-6
    k_restore, gamma_damp = 6e7, 4e4
    m_eff = 4/3*np.pi*rho*R0**3; v_max = 250.0
    E_tau, E_abs_scale = 20.0, 10.0
    gamma_e, P_cap, Ce = 1.25, 2e9, 2e-9
    burst_fraction, burst_cooldown = 0.35, 0.5
    eta_rad, tau_flash = 1e-3, 2e-3
    impulse_crit, impulse_window = 3e3, 0.01
    beta_acc, beta_burn, chi_abs = 5e-3, 2e-3, 4.0
    keys = ['R','v','E','Te','P','burst_flag','E_burst','I_flash',
            'pressure_impulse','m_np','rupture','drive','tau','T','A','R_floor']
    s = {key:np.zeros(N) for key in keys}
    R,v,E,Te,P,burst_flag,E_burst,I_flash,pressure_impulse,m_np,rupture = [s[k] for k in keys[:11]]
    R[0] = R0
    s['R_floor'][0] = Rmin
    s['tau'][0] = tau0; s['T'][0] = np.exp(-tau0)
    s['A'][0] = -np.expm1(-tau0); s['drive'][0] = laser[0]
    last_burst = -np.inf; armed = True
    for i in range(1, N):
        # ADDED: volume-equivalent radius bound. m_np has no specified SI units.
        # volume_gain therefore has inverse units of the original m_np proxy.
        # Previous-step mass is used consistently in optics and mechanics.
        floor = Rmin*(1 + volume_gain*m_np[i-1])**(1/3) if growing_floor else Rmin
        s['R_floor'][i] = floor
        # ADDED: absorption-only Beer-Lambert diagnostic, assuming constant
        # cross-sectional area and optical depth proportional to loading.
        # tau0=1 and volume_gain=1 are USER-ADJUSTABLE, UNCALIBRATED assumptions.
        s['tau'][i] = tau0*(1 + chi_abs*m_np[i-1])
        s['T'][i] = np.exp(-s['tau'][i])
        s['A'][i] = -np.expm1(-s['tau'][i])
        # Original A renamed drive: it is NOT the fractional absorption A(t).
        drive = min((1 + chi_abs*m_np[i-1])*laser[i], 3.0)
        s['drive'][i] = drive
        dE = (E_abs_scale*drive - E[i-1]/E_tau)*dt
        E[i] = max(E[i-1] + dE, 0.0)
        Te[i] = E[i]/max(Ce, 1e-30)
        V = 4/3*np.pi*max(R[i-1], floor)**3
        P[i] = min((gamma_e-1)*E[i]/V, P_cap)
        Fp = 4*np.pi*R[i-1]**2*(P[i]-P_amb)
        Fr = -k_restore*(R[i-1]-R0)
        Fd = -gamma_damp*v[i-1]
        a = (Fp + Fr + Fd)/m_eff
        v[i] = np.clip(v[i-1]+a*dt, -v_max, v_max)
        R[i] = np.clip(R[i-1]+v[i]*dt, floor, Rmax)
        collapsing = (R[i] < Rc) and (v[i] < 0)
        if collapsing and armed and (t[i]-last_burst > burst_cooldown):
            burst_flag[i] = 1
            E_burst[i] = burst_fraction*E[i]
            E[i] -= E_burst[i]
            last_burst = t[i]; armed = False
        if R[i] > 1.2*Rc:
            armed = True
        if burst_flag[i]:
            Aflash = eta_rad*E_burst[i]
            mask = t >= t[i]
            I_flash[mask] += (Aflash/tau_flash)*np.exp(-(t[mask]-t[i])/tau_flash)
        pressure_impulse[i] = pressure_impulse[i-1]+max(P[i]-P_amb, 0)*dt
        if pressure_impulse[i] > impulse_crit:
            rupture[i] = 1
        dm = beta_acc*(1-laser[i])-beta_burn*I_flash[i]
        m_np[i] = max(m_np[i-1]+dm*dt, 0)
    s.update(t=t,laser=laser)
    return s


def plot(s, path, title):
    fig, axes = plt.subplots(6, 2, figsize=(13, 16), sharex=True, layout='constrained')
    quantities = [('R',1e6,'Radius (µm)'),('laser',1,'Laser (on/off)'),
        ('E',1,'Stored energy (J)*'),('Te',1,'Electron temperature (K)*'),
        ('P',1e-6,'Pressure (MPa)'),('v',1,'Radial velocity (m/s)'),
        ('burst_flag',1,'Burst event'),('I_flash',1,'Flash power (W)*'),
        ('pressure_impulse',1,'Cumulative impulse (Pa s)'),('rupture',1,'Rupture indicator'),
        ('T',1,'Transmission fraction T(t)'),('A',1,'Absorption fraction A(t)')]
    for j,(ax,(key,scale,label)) in enumerate(zip(axes.flat,quantities)):
        ax.plot(s['t'],s[key]*scale, color='#195b9b',lw=0.75)
        ax.set_ylabel(label,fontsize=9)
        ax.set_title(f'({chr(97+j)})',loc='left',fontsize=10)
        ax.grid(alpha=0.22)
        if key=='R':
            ax.plot(s['t'],s['R_floor']*1e6,'--',color='#d14932',lw=1.3,label='Minimum radius bound')
            ax.legend(fontsize=8,loc='upper right')
            inset=ax.inset_axes([0.55,0.16,0.41,0.40])
            inset.plot(s['t'],s['R_floor']*1e6,color='#d14932',lw=1.1)
            inset.set_title('Minimum bound (µm)',fontsize=7)
            inset.tick_params(labelsize=6)
            inset.grid(alpha=0.2)
        if key in ('T','A'): ax.set_ylim(-0.03,1.03)
    for ax in axes[-1,:]: ax.set_xlabel('Time (s)')
    fig.suptitle(title+'\nAdded optics and radius-growth parameters are uncalibrated',fontsize=12)
    fig.supxlabel('*Nominal SI interpretation of original coefficients; not experimentally calibrated.',fontsize=9)
    fig.savefig(path.with_suffix('.png'),dpi=200)
    fig.savefig(path.with_suffix('.pdf'))
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('results'))
    parser.add_argument('--tau0',type=float,default=1.0)
    parser.add_argument('--volume-gain',type=float,default=1.0)
    args=parser.parse_args()
    if args.tau0 < 0 or args.volume_gain < 0: parser.error('Coefficients must be nonnegative')
    args.out.mkdir(parents=True,exist_ok=True)
    baseline=simulate(False,args.tau0,args.volume_gain)
    updated=simulate(True,args.tau0,args.volume_gain)
    assert np.allclose(updated['T']+updated['A'],1)
    assert np.all(updated['R'] >= updated['R_floor'])
    assert all(np.all(np.isfinite(v)) for v in updated.values())
    # Exact zero-growth regression: every original state must match.
    zero=simulate(True,args.tau0,0.0)
    original_keys=['R','v','E','Te','P','burst_flag','E_burst','I_flash','pressure_impulse','m_np','rupture']
    assert all(np.array_equal(zero[k],baseline[k]) for k in original_keys)
    report={'assumed_tau0':args.tau0,'assumed_volume_gain':args.volume_gain,
        'zero_growth_regression':'All 11 original state arrays match exactly',
        'max_absolute_changes':{k:float(np.max(np.abs(updated[k]-baseline[k]))) for k in original_keys},
        'baseline_bursts':int(baseline['burst_flag'].sum()),
        'updated_bursts':int(updated['burst_flag'].sum()),
        'R_floor_initial_um':float(updated['R_floor'][0]*1e6),
        'R_floor_final_um':float(updated['R_floor'][-1]*1e6),
        'radius_bound_fraction':float(np.mean((updated['R'][1:]==updated['R_floor'][1:]) | (updated['R'][1:]==30e-6))),
        'velocity_clip_fraction':float(np.mean(np.abs(updated['v'][1:])==250.0)),
        'pressure_cap_fraction':float(np.mean(updated['P'][1:]==2e9))}
    (args.out/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    for name,s in [('baseline',baseline),('updated',updated)]:
        keys=list(s)
        np.savetxt(args.out/f'{name}_data.csv',np.column_stack([s[k] for k in keys]),delimiter=',',header=','.join(keys),comments='')
        np.savez_compressed(args.out/f'{name}_data.npz',**s)
        plot(s,args.out/f'Figure_3_{name}',f'Figure 3 — {name.capitalize()} dynamics with added diagnostics')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
