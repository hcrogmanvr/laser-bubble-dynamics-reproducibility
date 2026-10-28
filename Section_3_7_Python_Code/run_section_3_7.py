"""Section 3.7 reproduction. Default verifies archived results and redraws Figures 9–11."""
from pathlib import Path
import argparse,subprocess,sys,json,csv
ROOT=Path(__file__).resolve().parent
def call(path):
    subprocess.run([sys.executable,str(ROOT/path)],cwd=ROOT,check=True)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['verify','full'],default='verify')
    args=parser.parse_args()
    if args.mode=='full':
        call('components/historical_single/code/run_controls.py')
        call('components/historical_ensemble/reproduce_defaults.py')
        call('components/historical_ensemble/diagnostics/coupling_check.py')
    call('components/historical_single/code/verify_results.py')
    # Check source identity against the supplied scripts, not only saved outputs.
    import hashlib
    manifest=json.loads((ROOT/'ORIGINAL_SOURCE_SHA256.json').read_text())
    for rel,expected in manifest.items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==expected,rel
    summary=ROOT/'components/historical_ensemble/diagnostics/coupling_results/matched_coupling_controls.json'
    rows=json.loads(summary.read_text())
    destination=ROOT/'outputs/tables';destination.mkdir(parents=True,exist_ok=True)
    # Export all source-level columns; manuscript Table 5 selects and rounds these values.
    if isinstance(rows,dict):
        (destination/'Table5_source_summary.json').write_text(json.dumps(rows,indent=2))
    else:
        keys=list(dict.fromkeys(k for r in rows for k in r))
        with (destination/'Table5_source_summary.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader()
            writer.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
    call('scripts/plot_control_figures.py')
    print('PASS Section 3.7: source identity, archived controls and Figures 9–11.')
if __name__=='__main__':main()
