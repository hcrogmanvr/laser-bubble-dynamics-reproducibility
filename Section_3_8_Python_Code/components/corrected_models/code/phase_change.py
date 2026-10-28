"""Conditional air/water-vapor cavity with finite-rate mass transfer and closed energy ledger.
An existing cavity is assumed; this does not model nucleation or plasma chemistry.
"""
from dataclasses import dataclass,asdict
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.integrate import quad

def saturation_pressure(T):
    # IAPWS SR1-86(1992), saturation curve, 273.16 K <= T <= 647.096 K.
    nu=1-np.asarray(T)/647.096
    terms=-7.85951783*nu+1.84408259*nu**1.5-11.7866497*nu**3+22.6807411*nu**3.5-15.9618719*nu**4+1.80122502*nu**7.5
    return 22.064e6*np.exp(647.096/np.asarray(T)*terms)

@dataclass
class PhaseParameters:
    duration_us:float=30.
    max_step_us:float=.04
    rtol:float=1e-7
    phase_on:bool=True
    accommodation:float=1e-3
    absorption_fraction:float=1e-4
    radius_um:float=5.44
    bead_radius_um:float=1.57
    ambient:float=101325.
    T0:float=300.
    rho:float=997.
    viscosity:float=.00089
    sigma:float=.07197
    Rair:float=287.
    Rv:float=461.5
    cvair:float=718.
    cvv:float=1410.
    latent:float=2.438e6
    bead_density:float=2500.
    bead_cp:float=700.
    conductivity:float=.6
    gas_conductivity:float=.026
    average_laser_W:float=.25
    frequency_Hz:float=1e5
    pulse_width_us:float=.1

def run_phase(p):
    R0=p.radius_um*1e-6;V0=4*np.pi/3*R0**3
    mv0=saturation_pressure(p.T0)*V0/(p.Rv*p.T0)
    mg=(p.ambient+2*p.sigma/R0-saturation_pressure(p.T0))*V0/(p.Rair*p.T0)
    Cw=4*np.pi/3*(p.bead_radius_um*1e-6)**3*p.bead_density*p.bead_cp
    Gb=4*np.pi*p.conductivity*p.bead_radius_um*1e-6
    def u_v(T):return p.latent-p.Rv*p.T0+p.cvv*(T-p.T0)
    def eval_state(y,power):
        r,v,Tg,Tw,m=y[:5];R=r*1e-6;mv=m*1e-15;V=4*np.pi/3*R**3;vd=4*np.pi*R*R*v
        pv=mv*p.Rv*Tg/V;pg=mg*p.Rair*Tg/V;pb=pg+pv
        # Two-sided Hertz-Knudsen flux; incoming and outgoing temperatures differ.
        md=4*np.pi*R*R*p.accommodation*(saturation_pressure(Tw)/np.sqrt(2*np.pi*p.Rv*Tw)-pv/np.sqrt(2*np.pi*p.Rv*Tg)) if p.phase_on else 0.
        Gg=4*np.pi*p.gas_conductivity*R
        Qwall=Gg*(Tw-Tg)
        Tflux=Tw if md>=0 else Tg
        hflux=p.latent+(p.cvv+p.Rv)*(Tflux-p.T0)
        Tgdot=(Qwall-pb*vd+md*(hflux-u_v(Tg)))/(mg*p.cvair+mv*p.cvv)
        Twdot=(power-Gb*(Tw-p.T0)-Qwall-md*hflux)/Cw
        acc=((pb-p.ambient-2*p.sigma/R-4*p.viscosity*v/R)/p.rho-1.5*v*v)/R
        f=np.zeros_like(y);f[:5]=[v,acc*1e-6,Tgdot*1e-6,Twdot*1e-6,md*1e9]
        f[5:]=np.array([power,Gb*(Tw-p.T0),16*np.pi*p.viscosity*R*v*v])*1e6 # pJ/us
        U=mg*p.cvair*(Tg-p.T0)+mv*u_v(Tg)+Cw*(Tw-p.T0)
        B=U+2*np.pi*p.rho*R**3*v*v+4*np.pi*p.sigma*R*R+p.ambient*V
        return f,dict(R_m=R,velocity_m_s=v,Tgas_K=Tg,Twall_K=Tw,vapor_kg=mv,partial_vapor_Pa=pv,pressure_Pa=pb,mass_flux_kg_s=md,total_energy_J=B)
    y=np.r_[p.radius_um,0.,p.T0,p.T0,mv0/1e-15,np.zeros(3)]
    ts=[0.];ys=[y.copy()];B0=eval_state(y,0)[1]['total_energy_J'];termination=None
    period=1e6/p.frequency_Hz
    edges=sorted(set([0.,p.duration_us]+[x for k in range(int(p.duration_us/period)+1) for x in [k*period,k*period+p.pulse_width_us] if 0<x<p.duration_us]))
    for a,b in zip(edges[:-1],edges[1:]):
        on=(.5*(a+b)%period)<p.pulse_width_us
        power=(p.average_laser_W/p.frequency_Hz)/(p.pulse_width_us*1e-6)*p.absorption_fraction*on
        def radius_domain(t,z):return z[0]-p.bead_radius_um
        def temperature_domain(t,z):return min(z[2]-273.16,z[3]-273.16,600-z[2],600-z[3])
        def mass_domain(t,z):return z[4]
        for f in [radius_domain,temperature_domain,mass_domain]:f.terminal=True;f.direction=-1
        sol=solve_ivp(lambda t,z:eval_state(z,power)[0],(a,b),y,method='Radau',rtol=p.rtol,atol=[1e-10,1e-8,1e-7,1e-7,1e-10,1e-12,1e-12,1e-12],max_step=p.max_step_us,events=[radius_domain,temperature_domain,mass_domain])
        if not sol.success:raise RuntimeError(sol.message)
        ts.extend(sol.t[1:]);ys.extend(sol.y.T[1:]);y=sol.y[:,-1]
        if sol.status==1:termination='domain_guard';break
    ts=np.array(ts);ys=np.array(ys)
    obs=[eval_state(z,0)[1] for z in ys];o={k:np.array([z[k] for z in obs]) for k in obs[0]}
    ledger=ys[:,5:]*1e-12 # Stored nJ ledger state actually pJ: derivative W * 1e6, dt us -> J*1e12.
    residual=o['total_energy_J']-B0-ledger[:,0]+ledger[:,1]+ledger[:,2]
    scale=max(ledger[:,0].max(),1e-18)
    summary=dict(parameters=asdict(p),completion_us=float(ts[-1]),termination=termination,min_radius_um=float(o['R_m'].min()*1e6),max_radius_um=float(o['R_m'].max()*1e6),max_Twall_K=float(o['Twall_K'].max()),max_Tgas_K=float(o['Tgas_K'].max()),vapor_mass_ratio_final=float(o['vapor_kg'][-1]/mv0),max_speed_m_s=float(abs(o['velocity_m_s']).max()),energy_residual_relative=float(abs(residual).max()/scale),absorbed_J=float(ledger[-1,0]),initial_vapor_kg=float(mv0),Cwall_J_K=float(Cw),Gwall_W_K=float(Gb))
    return dict(t_us=ts,states=ys,energy_residual_J=residual,**o),summary

def rayleigh_benchmark(rtol=1e-8):
    rho=1000.;dp=101325.;R0=30e-6;cut=.2
    # Scale R by R0 and time by R0 sqrt(rho/dp).
    fun=lambda t,y:[y[1],(-1-1.5*y[1]**2)/y[0]]
    def hit(t,y):return y[0]-cut
    hit.terminal=True;hit.direction=-1
    s=solve_ivp(fun,(0,2),[1.,0.],rtol=rtol,atol=rtol*1e-3,max_step=.01,events=hit,dense_output=True)
    exact=quad(lambda r:np.sqrt(1.5/(r**-3-1)),cut,1,epsabs=1e-12)[0]
    energy=s.y[0]**3*(1+1.5*s.y[1]**2)
    return dict(rtol=rtol,time_numeric_us=float(s.t[-1]*R0*np.sqrt(rho/dp)*1e6),time_exact_us=float(exact*R0*np.sqrt(rho/dp)*1e6),relative_time_error=float(abs(s.t[-1]/exact-1)),max_energy_invariant_error=float(abs(energy-1).max()),cutoff_R_over_R0=cut)
