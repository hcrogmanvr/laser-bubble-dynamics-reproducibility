from pathlib import Path
import json,numpy as np
from historical_controls import continuous
out=Path(__file__).resolve().parents[1]/'results';out.mkdir(exist_ok=True);rows=[]
for variant in [2,3]:
    for tol,step in [(1e-7,1.),(1e-9,.25)]:
        o,s=continuous(variant,rtol=tol,max_step=step);rows.append(s)
        if tol==1e-9:np.savez_compressed(out/f'unclipped_source{variant}.npz',**o)
(out/'source2_3_recheck.json').write_text(json.dumps(rows,indent=2))
