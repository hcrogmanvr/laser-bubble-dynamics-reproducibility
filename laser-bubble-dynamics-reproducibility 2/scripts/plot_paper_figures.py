"""Regenerate Figures 1–8 from the identified source versions and real workbook.

The original submitted Figure 4 is used here. The corrected variant-6
initialization plot that temporarily replaced it is deliberately excluded.
"""
from pathlib import Path
import argparse,importlib.util,sys,json,shutil,csv
import numpy as np
import pandas as pd
from scipy.stats import linregress,t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'outputs/figures';DATA=ROOT/'outputs/plotted_data'
for p in (FIG,DATA):p.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.fonttype':'none','pdf.fonttype':42})
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m
def finish(fig,name):
    import fitz
    for ext in ('pdf','svg'):fig.savefig(FIG/f'Figure{name}.{ext}')
    with fitz.open(FIG/f'Figure{name}.pdf') as pdf:
        pdf[0].get_pixmap(dpi=210).save(str(FIG/f'Figure{name}.png'))
    plt.close(fig)
def coordinates(name,t,series):
    points=np.column_stack([t]+[np.asarray(x) for x,_ in series])
    np.savetxt(DATA/f'Figure{name}.csv.gz',points,delimiter=',',
               header='time,'+','.join(label.replace('\n',' ').replace(',',';') for _,label in series),comments='')
def grid(name,t,series,shape,size,xlabel):
    fig,axs=plt.subplots(*shape,figsize=size,layout='constrained',sharex=True)
    for i,(ax,(values,label)) in enumerate(zip(np.asarray(axs).flat,series)):
        ax.plot(t,values,lw=.65,color='#195b9b');ax.set_ylabel(label,fontsize=9);ax.grid(alpha=.2)
        ax.text(.01,.96,f'({chr(97+i)})',transform=ax.transAxes,va='top',weight='bold',
                bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.5))
    for ax in np.asarray(axs).reshape(shape)[-1]:ax.set_xlabel(xlabel)
    coordinates(name,t,series)
    return fig,axs
def figure1():
    import fitz
    src=ROOT/'figures/editable/Figure1_editable.svg'
    shutil.copyfile(src,FIG/'Figure1.svg')
    with fitz.open(stream=src.read_bytes(),filetype='svg') as svg:
        with fitz.open('pdf',svg.convert_to_pdf()) as pdf:
            pdf.save(FIG/'Figure1.pdf');pdf[0].get_pixmap(dpi=250).save(str(FIG/'Figure1.png'))
def figure2():
    m=module(ROOT/'components/historical_single/source_originals/figure2_updated.py','paper_source2')
    a=m.simulate();np.savez_compressed(DATA/'Figure2_full_arrays.npz',**a)
    series=[(a['R']*1e6,'Radius (µm)'),(a['laser'],'Gate (0/1)'),(a['E'],'Stored energy\n(source units)'),
            (a['P']*1e-6,'Drive pressure\n(nominal MPa)'),(a['v'],'Wall speed (m/s)'),
            (abs(a['v'])/1500,'Mach number'),(a['burst_flag'],'Event flag (0/1)'),
            (a['E_burst'],'Event energy\n(source units)'),(a['I_cam'],'Detector proxy\n(source units/s)')]
    fig,axs=grid('2',a['t'],series,(9,1),(8,10),'Time (s)')
    minima=[float(a['R'][(a['t']>=lo)&(a['t']<hi)].min()*1e6) for lo,hi in [(60,90),(150,180),(240,270)]]
    ins=axs[0].inset_axes([.59,.23,.37,.51]);ins.plot([90,180,270],minima,'o-');ins.set(ylim=(0,1));ins.set_title('OFF minima (µm)',fontsize=7);ins.tick_params(labelsize=6)
    finish(fig,'2')
def figure3():
    m=module(ROOT/'components/historical_single/source_originals/figure3_updated.py','paper_source3')
    a=m.simulate(True);np.savez_compressed(DATA/'Figure3_full_arrays.npz',**a)
    series=[(a['R']*1e6,'Radius (µm)'),(a['laser'],'Gate (0/1)'),(a['E'],'Stored energy\n(source units)'),
            (a['Te'],'Temperature proxy\n(uncalibrated)'),(a['P']*1e-6,'Drive pressure\n(nominal MPa)'),
            (a['v'],'Wall speed (m/s)'),(a['burst_flag'],'Event flag (0/1)'),(a['I_flash'],'Flash proxy\n(source units/s)'),
            (a['pressure_impulse'],'Pressure impulse\n(nominal Pa s)'),(a['rupture'],'Pressure flag (0/1)'),
            (a['T'],'Transmission fraction'),(a['A'],'Absorption fraction')]
    fig,axs=grid('3',a['t'],series,(6,2),(10,11.5),'Time (s)')
    ins=axs[0,0].inset_axes([.57,.17,.40,.43]);ins.plot(a['t'],a['R_floor']*1e6,color='#d14932');ins.set_title('Prescribed bound (µm)',fontsize=7);ins.tick_params(labelsize=6)
    finish(fig,'3')
def figure4():
    # The archived states are outputs of the author-confirmed original source.
    # Rerun that exact source with `python reproduce.py --mode full` or --figure4-run.
    p=ROOT/'components/original_figure4/results/arrays.npz'
    with np.load(p) as archive:a={k:archive[k] for k in archive.files}
    series=[(a['R'].mean(axis=1)*1e6,'Mean radius (µm)'),(a['laser'],'Gate (0/1)'),
            (a['E'].mean(axis=1),'Mean energy\n(source units)'),(a['P'].mean(axis=1)*1e-6,'Mean drive\n(nominal MPa)'),
            (a['Te'].mean(axis=1),'Mean temperature\nproxy (nominal K)'),(a['cam'],'Camera proxy\n(arbitrary units)')]
    fig,axs=grid('4',a['t']*1e3,series,(6,1),(8,9),'Time (ms)')
    finish(fig,'4')
def ensemble(v):
    with np.load(ROOT/f'components/historical_ensemble/returned_arrays/figure{v}_returned_arrays.npz') as z:a={k:z[k] for k in z.files}
    series=[(a['R_mean']*1e6,'Mean radius (µm)'),(a['laser'],'Gate (0/1)'),(a['E_mean'],'Mean energy\n(nominal J)'),
            (a['P_mean_MPa'] if v==5 else a['P_over_MPa'],'Mean drive (MPa)' if v==5 else 'Mean overpressure\n(MPa)'),
            (a['Te_mean'],'Mean effective Tₑ (K)')]
    if v==6:series.extend([(a['xion_mean'],'Effective Saha fraction'),(a['M_mean'],'Loading index M')])
    series.append((a['I_cam'],'Detector response\n(arbitrary units)'))
    fig,_=grid(str(v),a['t']*1e3,series,(len(series),1),(8,9 if v==5 else 11),'Time (ms)');finish(fig,str(v))
def figure7():
    with np.load(ROOT/'components/historical_ensemble/returned_arrays/figure6_returned_arrays.npz') as z:t=z['t'];flags=z['rupt_frac']
    fig,ax=plt.subplots(figsize=(8,3.2),layout='constrained');ax.step(t*1e3,flags,where='post')
    ax.set(xlabel='Time (ms)',ylabel='Activated pressure flags / 8',ylim=(-.02,1.02));ax.grid(alpha=.2)
    coordinates('7',t,[(flags,'activated_flag_fraction')]);finish(fig,'7')
def growth_analysis():
    book=ROOT/'data/experimental/diameter_growth.xlsx';data={};summaries=[];points=[]
    for sheet,color in zip(pd.ExcelFile(book).sheet_names,['#c2185b','#1976d2','#278c50','#d42e2e']):
        df=pd.read_excel(book,sheet_name=sheet,header=None).iloc[3:,:2].apply(pd.to_numeric,errors='coerce').dropna()
        df.columns=['time_min','diameter_um'];df=df.sort_values('time_min').reset_index(drop=True)
        df['time_s']=df.time_min*60;df['D2_um2']=df.diameter_um**2;peak=int(df.diameter_um.idxmax());grow=df.iloc[:peak+1]
        fit=linregress(grow.time_s,grow.D2_um2);ci=student_t.ppf(.975,len(grow)-2)*fit.stderr
        row=dict(trace=sheet,n=len(grow),start_s=float(grow.time_s.iloc[0]),end_s=float(grow.time_s.iloc[-1]),
                 slope_um2_s=float(fit.slope),ci_low=float(fit.slope-ci),ci_high=float(fit.slope+ci),
                 intercept_um2=float(fit.intercept),r_squared=float(fit.rvalue**2))
        summaries.append(row);data[sheet]=(grow,row,color)
        for i,r in df.iterrows():points.append(dict(trace=sheet,**r.to_dict(),in_growth_window=i<=peak))
    pd.DataFrame(summaries).to_csv(DATA/'growth_fit_summary.csv',index=False)
    pd.DataFrame(points).to_csv(DATA/'growth_points.csv',index=False)
    return data,summaries
def figure8(which):
    data,summaries=growth_analysis()
    if which=='8A':
        red=next(k for k in data if 'Red' in k);df,r,c=data[red]
        fig,ax=plt.subplots(2,1,figsize=(8,7),layout='constrained');tt=np.linspace(0,60,200)
        ax[0].scatter(df.time_s,df.D2_um2,color=c,s=22,label='Measured red trace')
        ax[0].plot(tt,r['intercept_um2']+r['slope_um2_s']*tt,label=f'Empirical fit: {r["slope_um2_s"]:.2f} µm²/s')
        ax[1].plot(df.time_s,df.diameter_um,'o-',color=c)
        for a,label in zip(ax,['D² (µm²)','Diameter D (µm)']):a.set(xlim=(0,60),xlabel='Time (s)',ylabel=label)
        ax[0].legend(fontsize=8)
    else:
        fig,ax=plt.subplots(2,1,figsize=(8.5,7),layout='constrained')
        for i,(sheet,(df,r,c)) in enumerate(data.items()):
            tt=np.linspace(r['start_s'],r['end_s'],100);ax[0].scatter(df.time_s,df.D2_um2,color=c,s=16)
            ax[0].plot(tt,r['intercept_um2']+r['slope_um2_s']*tt,color=c,label=f'{sheet.split()[0]}: {r["slope_um2_s"]:.2f} µm²/s')
            ax[1].errorbar(i,r['slope_um2_s'],yerr=[[r['slope_um2_s']-r['ci_low']],[r['ci_high']-r['slope_um2_s']]],fmt='o',color=c,capsize=5)
        ax[0].set(xlabel='Time (s)',ylabel='D² (µm²)');ax[0].legend(fontsize=8,ncol=2)
        ax[1].set(xticks=range(4),xticklabels=[x['trace'].split()[0] for x in summaries],ylabel='Fitted slope (µm²/s)',
                  xlabel='Recorded trace; intervals are conditional OLS 95% CIs')
    for i,a in enumerate(ax):a.grid(alpha=.2);a.text(.01,.96,f'({chr(97+i)})',transform=a.transAxes,va='top',weight='bold')
    finish(fig,which)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--figures',nargs='+',default=['1','2','3','4','5','6','7','8A','8B']);args=parser.parse_args()
    for f in args.figures:
        if f in ['1','2','3','4','7']:globals()['figure'+f]()
        elif f in ['5','6']:ensemble(int(f))
        elif f in ['8A','8B']:figure8(f)
        else:raise ValueError(f)
        print('Rendered Figure',f,flush=True)
if __name__=='__main__':main()
