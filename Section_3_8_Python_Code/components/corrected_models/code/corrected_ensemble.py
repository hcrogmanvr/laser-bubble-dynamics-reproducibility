"""Numerical repair of the supplied effective-reservoir sources, not water-plasma validation.
SI physics internally; solve_ivp independent variable is microseconds, radius state micrometers.
No state clipping, artificial temperature kicks, pressure cap, or photon saturation.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from scipy.optimize import brentq

KB=1.380649e-23; H=6.62607015e-34; ME=9.1093837015e-31; EV=1.602176634e-19
NPART=1e25*(4*np.pi/3)*(2e-7)**3
EUNIT=1e-9

def saha(n,T,ion_eV=12.6):
    """Stable single-stage effective root and derivatives; no subtractive cancellation.
    The 12.6 eV source value is an effective input, not a species identification.
    """
    logA=np.log(2/n)+1.5*np.log(2*np.pi*ME*KB*T/H**2)-ion_eV*EV/(KB*T)
    # Equivalent root 2/(1+sqrt(1+4/A)); logaddexp prevents overflow.
    x=2*np.exp(-np.logaddexp(0, .5*np.logaddexp(0,np.log(4)-logA)))
    q=x*(1-x)/(2-x)
    return x,q*(1.5/T+ion_eV*EV/(KB*T*T)),q

def drude(n_e,R,nu=1e14,wavelength=1064e-9):
    """Single-interface, incoherent slab attenuation; no multiple internal reflections.
    nu and path 2R are assumptions. Returned power fractions sum to one.
    """
    w=2*np.pi*299792458/wavelength
    eps=1-n_e*EV**2/(8.8541878128e-12*ME*w*(w+1j*nu))
    z=np.sqrt(eps); reflect=np.abs((1-z)/(1+z))**2
    depth=4*w*R*z.imag/299792458
    transmit=(1-reflect)*np.exp(-depth)
    absorb=(1-reflect)*(-np.expm1(-depth))
    return reflect,transmit,absorb

@dataclass
class Parameters:
    variant:int=6
    N:int=8
    seed:int=2
    duration_us:float=1000.
    max_step_us:float=.1
    rtol:float=1e-7
    coupling:bool=True
    drude_on:bool=True
    one_temperature:bool=False
    power_multiplier:float=1.
    spacing_multiplier:float=1.
    rho:float=1000.
    viscosity:float=.001
    sigma:float=.072
    ambient:float=1e5
    Ce:float=3e-12
    Ci:float=2e-11
    tau_ei:float=2e-6
    tau_e:float=20e-6
    tau_i:float=80e-6
    ion_eV:float=12.6
    cutoff_um:float=.005
    mach_limit:float=.1

def geometry(p):
    rng=np.random.default_rng(p.seed)
    pos=rng.normal(scale=2e-6,size=(p.N,3))*p.spacing_multiplier
    rad=1.8e-7+rng.normal(scale=.2e-7,size=p.N)
    if np.any(rad<=0):raise ValueError('Nonpositive initial radius; seed/geometry invalid.')
    d=np.linalg.norm(pos[:,None]-pos[None,:],axis=-1)
    w=np.zeros_like(d); mask=~np.eye(p.N,dtype=bool)
    w[mask]=np.exp(-d[mask]/2.5e-6)/d[mask]
    w/=np.where(w.sum(axis=1,keepdims=True)>0,w.sum(axis=1,keepdims=True),1)
    return pos,rad,d,w

def evaluate(y,p,w,gate):
    a=y[:5*p.N].reshape(5,p.N); R=a[0]*1e-6; v=a[1]; Te=a[2]; Ti=a[3]; M=a[4]
    V=4*np.pi/3*R**3; vd=4*np.pi*R**2*v; n=NPART/V
    if p.variant==6:
        x,xt,q=saha(n,Te,p.ion_eV); xr=3*q/R; ion=NPART*p.ion_eV*EV
    else:
        x=expit((Te-12000)/2000);xt=x*(1-x)/2000;xr=np.zeros_like(R);ion=0.
    pe=n*KB*x*Te; ph=n*KB*Ti; pb=pe+ph
    pc=.25*(w@(pb-pb.mean())) if p.coupling and p.N>1 else np.zeros(p.N)
    accel=((pb+pc-p.ambient-2*p.sigma/R-4*p.viscosity*v/R)/p.rho-1.5*v*v)/R
    af=.35*(-np.expm1(-R/2e-6))*(.3+.7*(-np.expm1(-R/1e-7)))
    # Bounded optical fraction, implemented smoothly rather than clipping states.
    if p.variant==6: af=-np.expm1(-af*(1+.25*M))
    refl,tr,ap=drude(x*n,R) if p.drude_on and p.variant==6 else (np.zeros(p.N),np.ones(p.N),np.zeros(p.N))
    inc=gate*2e-3*p.power_multiplier
    # Plasma absorption and transmitted absorber heating both enter the thermal budget.
    Qabs=inc*(ap+tr*af)
    qe=p.Ce*(Te-300)/p.tau_e; qi=p.Ci*(Ti-300)/p.tau_i
    exchange=p.Ce*(Te-Ti)/p.tau_ei
    dTe=(.9*Qabs-exchange-qe-pe*vd-ion*xr*v)/(p.Ce+ion*xt)
    dTi=(.1*Qabs+exchange-qi-ph*vd)/p.Ci
    if p.one_temperature:
        common=(Qabs-qe-qi-pb*vd-ion*xr*v)/(p.Ce+p.Ci+ion*xt)
        dTe=common;dTi=common
    dm=4e5*expit((pb-2e6)/5e5)*(1-M/10) if p.variant==6 else np.zeros(p.N)
    f=np.zeros_like(y);f[:5*p.N]=np.array([v,accel*1e-6,dTe*1e-6,dTi*1e-6,dm*1e-6]).ravel()
    # Global integrated ledgers: absorbed, bath loss, viscous loss, coupling work, removed burst energy.
    rates=np.array([Qabs.sum(),(qe+qi).sum(),(16*np.pi*p.viscosity*R*v*v).sum(),(pc*vd).sum(),0])
    f[5*p.N:]=rates*1e-6/EUNIT
    internal=p.Ce*Te+p.Ci*Ti+ion*x
    mech=2*np.pi*p.rho*R**3*v*v+4*np.pi*p.sigma*R**2+p.ambient*V
    return f,dict(R=R,v=v,Te=Te,Ti=Ti,M=M,x=x,pb=pb,pabs=Qabs,internal=internal,mechanical=mech,refl=refl,tr=tr,ap=ap)

def apply_event(y,i,p,w):
    """Phenomenological source threshold retained only as a conservative diagnostic jump.
    Remove 35% of excess thermal+ionization energy and 25% of liquid kinetic energy.
    No extra energy is introduced. 8% is assigned to the optical proxy, 92% unresolved loss.
    """
    z=y.copy();a=z[:5*p.N].reshape(5,p.N); before=evaluate(z,p,w,0)[1]
    e=before['internal'][i]
    b=z.copy();b[:5*p.N].reshape(5,p.N)[2:4,i]=300
    e0=evaluate(b,p,w,0)[1]['internal'][i]
    loss=.35*max(e-e0,0)
    old=a[2:4,i].copy()
    if loss>0:
        def objective(s):
            a[2:4,i]=300+s*(old-300)
            return evaluate(z,p,w,0)[1]['internal'][i]-(e-loss)
        scale=brentq(objective,0,1,xtol=1e-12);a[2:4,i]=300+scale*(old-300)
    kinetic=2*np.pi*p.rho*before['R'][i]**3*a[1,i]**2
    a[1,i]*=np.sqrt(.75)
    removal=loss+.25*kinetic;z[-1]+=removal/EUNIT
    after=evaluate(z,p,w,0)[1]
    residual=(before['internal'].sum()+before['mechanical'].sum()-after['internal'].sum()-after['mechanical'].sum())-removal
    return z,removal,residual

def run(p):
    pos,r0,dist,w=geometry(p)
    y=np.r_[np.array([r0*1e6,np.zeros(p.N),np.full(p.N,300.),np.full(p.N,300.),np.zeros(p.N)]).ravel(),np.zeros(5)]
    initial=evaluate(y,p,w,0)[1];B0=(initial['internal']+initial['mechanical']).sum()
    ts=[0.];ys=[y.copy()];events=[];reset_residual=[]; current=0.; evaluations=0;terminated=None
    armed=r0>1.5e-7
    edges=sorted(set([0.,p.duration_us]+[x for k in range(int(p.duration_us/80)+1) for x in (80*k,80*k+60) if 0<x<p.duration_us]))
    for lo,hi in zip(edges[:-1],edges[1:]):
        gate=float((.5*(lo+hi)%80)<60)
        while current<hi-1e-10:
            ev=[];which=[]
            for i in range(p.N):
                if armed[i]:
                    def cross(t,z,i=i):return z[i]-.15
                    cross.direction=-1;cross.terminal=True;ev.append(cross);which.append(('burst',i))
                else:
                    def rearm(t,z,i=i):return z[i]-.18
                    rearm.direction=1;rearm.terminal=True;ev.append(rearm);which.append(('arm',i))
            def small(t,z):return np.min(z[:p.N])-p.cutoff_um
            small.direction=-1;small.terminal=True;ev.append(small);which.append(('radius_domain',-1))
            def mach(t,z):return p.mach_limit*1500-np.max(np.abs(z[p.N:2*p.N]))
            mach.direction=-1;mach.terminal=True;ev.append(mach);which.append(('mach_domain',-1))
            atol=np.r_[np.full(p.N,1e-10),np.full(p.N,1e-7),np.full(2*p.N,1e-7),np.full(p.N,1e-10),np.full(5,1e-12)]
            sol=solve_ivp(lambda t,z:evaluate(z,p,w,gate)[0],(current,hi),y,method='Radau',rtol=p.rtol,atol=atol,max_step=p.max_step_us,events=ev)
            if not sol.success:raise RuntimeError(sol.message)
            evaluations+=sol.nfev;ts.extend(sol.t[1:]);ys.extend(sol.y.T[1:]);current=sol.t[-1];y=sol.y[:,-1].copy()
            hit=next((j for j,x in enumerate(sol.t_events) if len(x)),None)
            if hit is None:break
            kind,i=which[hit]
            if kind.endswith('domain'):terminated=kind;break
            if kind=='arm':armed[i]=True
            else:
                y,loss,res=apply_event(y,i,p,w);events.append((current,i,loss));reset_residual.append(res);armed[i]=False
                ts.append(current);ys.append(y.copy())
        if terminated:break
    ts=np.array(ts);ys=np.array(ys)
    obs={k:[] for k in initial}
    for t,z in zip(ts,ys):
        o=evaluate(z,p,w,float((t%80)<60))[1]
        for k,v in o.items():obs[k].append(v)
    obs={k:np.array(v) for k,v in obs.items()}
    B=(obs['internal']+obs['mechanical']).sum(axis=1)
    led=ys[:,5*p.N:]*EUNIT
    residual=B-B0-led[:,0]+led[:,1]+led[:,2]-led[:,3]+led[:,4]
    scale=max(np.max(np.abs(led[:,0])),np.max(np.abs(B-B0)),1e-18)
    mask=~np.eye(p.N,dtype=bool)
    ratio=np.max((obs['R'][:,:,None]+obs['R'][:,None,:])[:,mask]/dist[mask]) if p.N>1 else 0.
    out=dict(t_us=ts,states=ys,positions_m=pos,distances_m=dist,events=np.array(events).reshape(-1,3),energy_residual_J=residual,**obs)
    summary=dict(parameters=asdict(p),completion_us=float(current),terminated=terminated,events=len(events),min_radius_um=float(obs['R'].min()*1e6),max_radius_um=float(obs['R'].max()*1e6),peak_Te_K=float(obs['Te'].max()),peak_Ti_K=float(obs['Ti'].max()),peak_speed_m_s=float(abs(obs['v']).max()),max_ionization=float(obs['x'].max()),absorbed_J=float(led[-1,0]),removed_J=float(led[-1,4]),energy_residual_relative=float(abs(residual).max()/scale),event_reset_residual_J=float(max(abs(np.array(reset_residual)),default=0)),max_pair_overlap_ratio=float(ratio),nfev=evaluations)
    return out,summary
