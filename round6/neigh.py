from common import *
import json
d = json.load(open(OUT / 'contam' / 'contamination.json'))
rows = [dict(source_id=r['source_id'], sep=r['sep_arcsec'], pa=r['pa_deg'], tmag=r['Tmag_est'], G=r['G'],
             need1150=r['depth_needed_1150ppm'], need284=r['depth_needed_284ppm']) for r in d['neighbours'] if r['sep_arcsec'] < 110]
json.dump(rows, open(OUT / 'data' / 'neighbours_local.json', 'w'), indent=1); print(len(rows))
