#!/usr/bin/env python3
"""Entry point for the manuscript's figures, tables and numerical verification."""
from pathlib import Path
import argparse,subprocess,sys,time,json,shutil,platform
ROOT=Path(__file__).resolve().parent
LOGS=ROOT/'outputs/logs'
def step(label,script,*args):
    LOGS.mkdir(parents=True,exist_ok=True);dest=LOGS/(label+'.log')
    print(label,flush=True);started=time.time()
    with dest.open('w') as f:
        result=subprocess.run([sys.executable,str(ROOT/script),*map(str,args)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    if result.returncode:
        print(dest.read_text()[-8000:],file=sys.stderr)
        raise SystemExit(f'{label} failed; see {dest.relative_to(ROOT)}')
    print(f'  PASS ({time.time()-started:.1f} s)',flush=True)
def corrected_plots():
    step('corrected_figures','components/corrected_models/code/make_figures.py')
    source=ROOT/'components/corrected_models/figures';dest=ROOT/'outputs/figures';dest.mkdir(parents=True,exist_ok=True)
    for p in source.iterdir():
        if p.name.startswith(('Figure12_','Figure13_','FigureS1_')):shutil.copyfile(p,dest/p.name)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['figures','verify','full'],default='figures',
                        help='figures: redraw archived outputs; verify: redraw and check; full: rerun every numerical suite first')
    parser.add_argument('--figure',choices=['1','2','3','4','5','6','7','8A','8B','9','10','11','12','13'],
                        help='Redraw one figure; --mode full reruns only the relevant numerical family')
    parser.add_argument('--figure4-run',action='store_true',help='Recalculate the author-confirmed original Figure 4 source')
    args=parser.parse_args();started=time.time()
    full=args.mode=='full';fig=args.figure
    if full:
        if fig is None or fig in ['2','3','9','10','11']:
            step('historical_single_controls','components/historical_single/code/run_controls.py')
        if fig is None or fig in ['5','6','7']:
            step('historical_ensemble_defaults','components/historical_ensemble/reproduce_defaults.py')
            step('historical_ensemble_coupling','components/historical_ensemble/diagnostics/coupling_check.py')
            step('historical_density_observer','components/historical_ensemble/diagnostics/density_check.py')
        if fig is None or fig in ['12','13']:
            step('corrected_models_all_cases','components/corrected_models/code/run_suite.py')
            step('continuous_sources_recheck','components/corrected_models/code/recheck_sources.py')
    if args.figure4_run or (full and (fig is None or fig=='4')):
        step('original_figure4','components/original_figure4/run_and_export.py')
    if fig is None:
        step('paper_figures_1_to_8','scripts/plot_paper_figures.py')
        step('control_figures_9_to_11','scripts/plot_control_figures.py')
        corrected_plots()
    elif fig in ['9','10','11']:step('control_figures_9_to_11','scripts/plot_control_figures.py')
    elif fig in ['12','13']:corrected_plots()
    else:step('paper_figure_'+fig,'scripts/plot_paper_figures.py','--figures',fig)
    if fig is None:
        step('calculated_tables','scripts/export_calculated_tables.py')
        step('reference_calculations','scripts/reference_calculations.py')
        if args.mode in ['verify','full']:step('verification','scripts/verify_repository.py')
    manifest=dict(mode=args.mode,figure=fig,figure4_rerun=bool(args.figure4_run or (full and (fig is None or fig=='4'))),elapsed_s=time.time()-started,
                  python=platform.python_version(),successful=True,
                  scope='Numerical reproduction and verification; no experimental plasma/flash validation claimed.')
    (ROOT/'outputs/run_manifest.json').write_text(json.dumps(manifest,indent=2))
    print('Completed. Figures, plotted coordinates, tables and logs are in outputs/.',flush=True)
if __name__=='__main__':main()
