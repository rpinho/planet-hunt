import csv,sys,glob,os
keys=['tceid','tce_plnt_num','tce_period','tce_period_err','tce_time0bt','tce_time0bt_err','tce_duration','tce_depth','tce_prad','tce_model_snr','tce_max_mult_ev','tce_num_transits','tce_ror','tce_impact','tce_dikco_msky','tce_dicco_msky','tce_ws_mesmed','tce_bin_oedp_stat','tce_hap_stat','tce_ptemp_stat']
for f in sorted(glob.glob(os.path.join(sys.argv[1],'tce_*.csv'))):
    rows=list(csv.DictReader(open(f)))
    print('##',os.path.basename(f), 'rows',len(rows))
    for r in rows:
        print('   ', {k:r.get(k) for k in keys if k in r})
