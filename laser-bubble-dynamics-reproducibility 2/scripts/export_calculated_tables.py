"""Build numerical Tables 2–7 from their identified calculation records."""
from pathlib import Path
import csv,json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/tables';OUT.mkdir(parents=True,exist_ok=True)
def read(p):return json.loads((ROOT/p).read_text())
def save(n,rows):pd.DataFrame(rows).to_csv(OUT/f'Table{n}_calculated.csv',index=False)
def main():
    fit=ROOT/'outputs/plotted_data/growth_fit_summary.csv'
    if not fit.exists():
        from plot_paper_figures import growth_analysis
        growth_analysis()
    save(2,pd.read_csv(fit).to_dict('records'))
    legacy=read('components/historical_single/results/legacy_controls.json')
    rows=[]
    for a in [r for r in legacy if r['model'].startswith('Figure 2') and r['test']=='timestep']:
        b=next(r for r in legacy if r['model'].startswith('Figure 3') and r['test']=='timestep' and r['dt_s']==a['dt_s'])
        rows.append(dict(timestep_ms=a['dt_s']*1e3,events_source2=a['events'],events_source3=b['events'],
                         flash_ratio_source2=a['sampled_flash_integral_proxy']/a['radiated_event_energy_proxy'],
                         flash_ratio_source3=b['sampled_flash_integral_proxy']/b['radiated_event_energy_proxy']))
    save(3,rows)
    rows=read('components/historical_single/results/continuous_convergence.json')
    save(4,[dict(source=int(r['model'].split()[1]),relative_tolerance=r['rtol'],max_step_s=r['max_step_s'],
                 max_displacement_nm=r['R_max_excursion_nm'],difference_from_finest_nm=r['max_radius_difference_from_finest_nm'],events=r['events']) for r in rows])
    controls=read('components/historical_ensemble/diagnostics/coupling_results/matched_coupling_controls.json')['cases']
    save(5,[dict(source=r['figure'],coupling=r['coupling'],triggers=r['events'],
                 peak_mean_Te_K=r['peak_stored_mean_Te_K'],peak_pressure_MPa=r['peak_stored_mean_pressure_MPa'],
                 mean_active_phase_coherence=r['active_phase_coherence_mean']) for r in controls])
    controls=read('components/corrected_models/results/corrected_ensemble_controls.json');rows=[]
    for v in [5,6]:
        for label,test in [('Two T; coupling ON','refinement'),('Two T; coupling OFF','coupling_off'),('One T; coupling ON','one_temperature')]:
            r=next(x for x in controls if x['parameters']['variant']==v and x['test']==test and (test!='refinement' or x['parameters']['rtol']==1e-8))
            rows.append(dict(variant=v,control=label,events=r['events'],peak_temperature_K=r['peak_Te_K'],
                             minimum_radius_nm=r['min_radius_um']*1e3,energy_residual=r['energy_residual_relative']))
    save(6,rows)
    controls=read('components/corrected_models/results/phase_controls.json');rows=[]
    for test in ['refinement','phase_off','accommodation_low','accommodation_high','absorption_low','absorption_high']:
        r=next(x for x in controls if x['test']==test and (test!='refinement' or x['parameters']['rtol']==1e-8))
        rows.append(dict(case='Reference ON' if test=='refinement' else test.replace('_',' '),peak_radius_um=r['max_radius_um'],
                         peak_wall_temperature_K=r['max_Twall_K'],final_vapor_mass_ratio=r['vapor_mass_ratio_final'],
                         energy_residual=r['energy_residual_relative']))
    save(7,rows)
    print('Exported calculated Tables 2–7. Table 1 is the documented input/provenance table, not a fitted output.')
if __name__=='__main__':main()
