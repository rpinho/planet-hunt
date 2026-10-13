import sys
from astropy.io import fits
for f in sys.argv[1:]:
    try:
        h=fits.getheader(f)
        print(f.split('/')[-1], h.get('NAXIS1'),h.get('NAXIS2'), 'SURVEY',h.get('SURVEY'),'PLATEID',h.get('PLATEID'),'EMUL',h.get('EMULSION'),'FILT',h.get('FILTER'),'DATE',h.get('DATE-OBS'),'EXP',h.get('EXPOSURE'),'EPOCH',h.get('EPOCH'), 'PLTSCALE',h.get('PLTSCALE'),'XPIX',h.get('XPIXELSZ'), 'CDELT', h.get('CDELT1'), h.get('CD1_1'))
    except Exception as e: print(f, 'ERR', e)
