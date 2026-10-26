"""Verify provenance, source reproduction, printed numerical tables and ledgers."""
from pathlib import Path
import json,hashlib,csv,re,subprocess,sys
import numpy as np
import pandas as pd
from PIL import Image
from xml.etree import ElementTree
ROOT=Path(__file__).resolve().parents[1];checks=[]
def passed(name,detail):checks.append(dict(check=name,status='PASS',detail=detail))
def read(name):return json.loads((ROOT/name).read_text())
def printed_equal(actual,printed):
    s=printed.strip().lower().replace('−','-');value=float(s)
    mantissa,exponent=(s.split('e')+[0])[:2] if 'e' in s else (s,0)
    decimals=len(mantissa.split('.')[1]) if '.' in mantissa else 0
    tolerance=.501*10**(int(exponent)-decimals)+1e-14*max(abs(value),1e-12)
    assert abs(float(actual)-value)<=tolerance,(actual,printed,tolerance)
def compare_archived(source,actual,label):
    with np.load(source,allow_pickle=False) as a,np.load(actual,allow_pickle=False) as b:
        for key in a.files:
            assert key in b.files,(label,key)
            assert a[key].shape==b[key].shape,(label,key,'shape')
            if a[key].dtype.kind=='b':
                assert np.array_equal(a[key],b[key]),(label,key)
                continue
            scale=max(float(np.max(abs(a[key]))) if a[key].size else 0.,1e-300)
            difference=float(np.max(abs(a[key]-b[key]))) if a[key].size else 0.
            assert difference/scale<1e-10,(label,key,difference/scale)
    passed(label,'Every archived array reproduced within 1e-10 of its peak scale.')
def main():
    records=read('docs/figure_registry.json');assert len(records)==14
    for r in records:assert hashlib.sha256((ROOT/r['manuscript_image']).read_bytes()).hexdigest()==r['sha256'],r['figure']
    passed('Manuscript image provenance','All 14 exact manuscript image files match their recorded hashes.')
    formulas=read('docs/equation_index.json');assert len(formulas)==108
    for r in formulas:ElementTree.parse(ROOT/r['ooxml'])
    passed('Equation inventory','All 108 native Office Math XML expressions are preserved and parse.')
    assert len(read('docs/table_index.json'))==7
    for f in ROOT.joinpath('outputs/figures').glob('*.png'):
        with Image.open(f) as im:im.verify()
    for n in ['1','2','3','4','5','6','7','8A','8B','9','10','11']:
        assert (ROOT/f'outputs/figures/Figure{n}.png').is_file(),n
    for n in ['12','13']:assert list(ROOT.joinpath('outputs/figures').glob(f'Figure{n}_*.png')),n
    passed('Figure rendering','Every numbered manuscript figure has a readable regenerated PNG.')
    for n in [2,3]:compare_archived(ROOT/f'data/reference_arrays/historical_figure{n}.npz',ROOT/f'outputs/plotted_data/Figure{n}_full_arrays.npz',f'Original Figure {n} states')
    compare_archived(ROOT/'data/reference_arrays/figure4_submitted_archived.npz',ROOT/'components/original_figure4/results/arrays.npz','Original Figure 4 states')
    for n in [5,6]:compare_archived(ROOT/f'data/reference_arrays/historical_figure{n}.npz',ROOT/f'components/historical_ensemble/returned_arrays/figure{n}_returned_arrays.npz',f'Original ensemble {n} returned arrays')
    summary=read('components/original_figure4/results/run_summary.json')
    assert summary['total_burst_flags']==8 and summary['event_rows']==[0]
    assert np.isclose(summary['max_camera_signal'],.0028,rtol=1e-12)
    passed('Submitted Figure 4 interpretation','Eight first-update flags and camera maximum 0.0028; source version unchanged.')
    subprocess.run([sys.executable,str(ROOT/'components/historical_single/code/verify_results.py')],check=True)
    passed('Historical source controls','Source hashes, equilibrium check, refinement records and saved controls pass their original verifier.')
    # Compare calculated quantities to the manuscript at its printed precision.
    def manuscript(n):return list(csv.reader((ROOT/f'docs/tables/Table{n}_manuscript.csv').open()))[1:]
    tables={n:pd.read_csv(ROOT/f'outputs/tables/Table{n}_calculated.csv') for n in range(2,8)}
    for raw,r in zip(manuscript(2),tables[2].to_dict('records')):
        assert raw[0].lower()==r['trace'].split()[0].lower()
        printed_equal(r['n'],raw[1]);printed_equal(r['slope_um2_s'],raw[3]);printed_equal(r['r_squared'],raw[5])
        lo,hi=re.split('[–—]',raw[4]);printed_equal(r['ci_low'],lo);printed_equal(r['ci_high'],hi)
    for raw,r in zip(manuscript(3),tables[3].to_dict('records')):
        for value,s in zip([r['timestep_ms'],r['events_source2'],r['flash_ratio_source2'],r['flash_ratio_source3']],raw):printed_equal(value,s)
        assert r['events_source2']==r['events_source3']
    for raw,r in zip(manuscript(4),tables[4].to_dict('records')):
        for key,s in zip(['source','relative_tolerance','max_step_s','max_displacement_nm','difference_from_finest_nm'],raw):printed_equal(r[key],s)
        assert r['events']==0
    for raw,r in zip(manuscript(5),tables[5].to_dict('records')):
        for key,s in zip(['source','triggers','peak_mean_Te_K','peak_pressure_MPa','mean_active_phase_coherence'],[raw[0].split()[-1],*raw[2:]]):printed_equal(r[key],s)
    for raw,r in zip(manuscript(6),tables[6].to_dict('records')):
        for key,s in zip(['variant','events','peak_temperature_K','minimum_radius_nm'],[raw[0],*raw[2:5]]):printed_equal(r[key],s)
        assert r['energy_residual']<1e-6
    for raw,r in zip(manuscript(7),tables[7].to_dict('records')):
        for key,s in zip(['peak_radius_um','peak_wall_temperature_K','final_vapor_mass_ratio'],raw[1:4]):printed_equal(r[key],s)
        assert r['energy_residual']<1e-6
    passed('Tables 2–7','Calculated counts, growth slopes/intervals, extrema and control quantities agree at manuscript precision; energy residuals pass the stated tolerance.')
    # Residual magnitudes can vary below solver tolerance; retain both archived
    # printed values and freshly calculated values without forcing agreement.
    ensemble=read('components/corrected_models/results/corrected_ensemble_controls.json')
    phase=read('components/corrected_models/results/phase_controls.json')
    assert len(ensemble)==26 and len(phase)==8
    assert all(r['terminated'] is None and r['energy_residual_relative']<1e-6 for r in ensemble)
    assert all(r['energy_residual_relative']<1e-6 for r in phase)
    passed('Corrected numerical ledgers','26 ensemble and eight phase cases complete within the declared energy-ledger tolerance.')
    benchmark=read('components/corrected_models/results/closure_checks.json')
    ray=[r for r in benchmark if r['test']=='Rayleigh finite-radius benchmark'][-1]
    assert ray['relative_time_error']<1e-8
    passed('Independent radial benchmark',f'Finite-radius Rayleigh relative time error {ray["relative_time_error"]:.3g}.')
    out=ROOT/'outputs/verification.json';out.write_text(json.dumps(dict(successful=True,checks=checks,scope='Numerical reproduction and verification, not experimental plasma validation.'),indent=2))
    for c in checks:print(c['status'],c['check'])
if __name__=='__main__':main()
