"""Section 3.8 corrected ensembles and shared independent/phase checks."""
from pathlib import Path
import argparse,subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
def call(path):
    subprocess.run([sys.executable,str(ROOT/path)],cwd=ROOT,check=True)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['verify','full'],default='verify')
    args=parser.parse_args()
    if args.mode=='full':
        call('components/corrected_models/code/run_suite.py')
        call('components/corrected_models/code/recheck_sources.py')
    manifest=json.loads((ROOT/'ORIGINAL_SOURCE_SHA256.json').read_text())
    for rel,expected in manifest.items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==expected,rel
    res=ROOT/'components/corrected_models/results'
    rows=json.loads((res/'corrected_ensemble_controls.json').read_text())
    assert len(rows)==26
    assert all(row['terminated'] is None for row in rows)
    assert max(row['energy_residual_relative'] for row in rows)<1e-6
    fine=[row for row in rows if row['name'] in ('figure5_refine2','figure6_refine2')]
    assert len(fine)==2 and all(row['events']==5 for row in fine)
    checks=json.loads((res/'closure_checks.json').read_text())
    saha=next(row for row in checks if row['test'].startswith('Hydrogen Saha'))
    assert saha['max_absolute_fraction_error']<1e-12
    phase=json.loads((res/'phase_controls.json').read_text())
    assert len(phase)==8 and max(row['energy_residual_relative'] for row in phase)<1e-6
    call('components/corrected_models/code/make_figures.py')
    print('PASS Section 3.8: source identity, 26 archived ensemble controls and Figure 12.')
    print('Shared suite also includes 8 separate phase cases and Figure 13 (Section 3.9 in the reorganized manuscript).')
if __name__=='__main__':main()
