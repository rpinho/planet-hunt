import csv,sys,math
ra0,de0=324.38941957,-77.97914675
def sex2deg(r,d):
    try: return float(r),float(d)
    except: pass
    h,m,s=[float(x) for x in r.split(':')]; ra=15*(h+m/60+s/3600)
    sg=-1 if d.strip().startswith('-') else 1; dd,dm,ds=[abs(float(x)) for x in d.split(':')]
    return ra, sg*(dd+dm/60+ds/3600)
for fn in sys.argv[1:]:
    with open(fn, newline='') as fh:
        for row in csv.DictReader(fh):
            try: ra,de=sex2deg(row['RA'],row['Dec'])
            except Exception as e: continue
            sep=math.degrees(math.acos(min(1,math.sin(math.radians(de))*math.sin(math.radians(de0))+math.cos(math.radians(de))*math.cos(math.radians(de0))*math.cos(math.radians(ra-ra0)))))*3600
            if sep<300: print(fn.split('/')[-1], round(sep,1), row.get('TOI') or row.get('CTOI'), row['TIC ID'], row.get('Period (days)'))
