import csv,sys
for f in sys.argv[1:]:
    for r in csv.DictReader(open(f)):
        print(f.split('/')[-1], r['tceid'], 'sectors=',r.get('sectors'), 'tce_sectors=',r.get('tce_sectors'), 'lastUpdate=',r.get('lastUpdate'))
