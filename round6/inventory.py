"""Round 6 / TOI-4342: which TESS sectors see the star, and which MAST products exist."""
import sys, json
sys.path.insert(0, '.')
from planet_hunt import config  # noqa: cache redirect
import numpy as np
from tess_stars2px import tess_stars2px_function_entry as t2p
from astroquery.mast import Observations, conf
conf.timeout = 600
TIC, RA, DEC = 354944123, 324.3894166666666, -77.97914722222222
out = {}
res = t2p(TIC, RA, DEC)
rows = []
for i in range(len(res[3])):
    rows.append(dict(sector=int(res[3][i]), camera=int(res[4][i]), ccd=int(res[5][i]),
                     col=float(res[6][i]), row=float(res[7][i])))
out['tess_point'] = rows
print('tess-point sectors:', [r['sector'] for r in rows])
obs = Observations.query_criteria(target_name=str(TIC), obs_collection=['TESS', 'HLSP'])
df = obs.to_pandas()
cols = ['obs_collection', 'provenance_name', 'obs_id', 'sequence_number', 't_exptime', 'dataproduct_type', 't_min', 't_max']
df = df[[c for c in cols if c in df.columns]].sort_values(['sequence_number', 'provenance_name'])
print(df.to_string())
out['mast'] = df.to_dict(orient='records')
json.dump(out, open('runs/round6_toi4342/data/inventory.json', 'w'), indent=1, default=str)
