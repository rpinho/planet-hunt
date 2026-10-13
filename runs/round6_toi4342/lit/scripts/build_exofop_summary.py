import json, csv, os, sys, glob, datetime

BASE = sys.argv[1]  # .../exofop
T = json.load(open(os.path.join(BASE, 'raw', 'target.json')))
W = json.load(open(os.path.join(BASE, 'sg1_windows.json')))

# ---------------- imaging ----------------
imaging = [
    {
        "telescope": "Gemini South (8 m)", "instrument": "Zorro", "type": "Speckle",
        "filter": "832 nm (40 nm width); only the 832 nm channel curve was uploaded",
        "date_ut": "2021-07-23", "user": "furlan (PI Howell)", "exofop_tag": "99434",
        "exofop_contrast_summary": "delta 5.87 mag @ 0.5\"",
        "pixel_scale_arcsec": 0.01, "psf_arcsec": 0.02,
        "contrast_dmag": {"0.1": 4.90, "0.5": 5.87, "1.0": 6.08, "1.17_max_sep": 6.15},
        "contrast_source": "TOI4342I-ef20210723-832_sensitivity.dat, smooth-fit section (annulus section gives 4.74 @0.10, 5.88 @0.525, 6.09 @1.025)",
        "files": ["TOI4342I-ef20210723-832_sensitivity.dat", "TOI4342I-ef20210723-832_sensitivity.pdf",
                  "TOI4342I-ef20210723-20210723_832_plot.pdf", "TOI4342I-ef20210723-832.fits"],
        "file_format": "sensitivity.dat: header lines start with '#'. Section 1 '# annulus range&center vs. measured relative mag. limit': 4 cols = sep_inner_arcsec, sep_outer_arcsec, sep_center_arcsec, delta_mag(832nm, 5-sigma limit). Section 2 '# fit of angular separation vs. relative mag. limit': 2 cols = sep_arcsec (0.00-1.17 step 0.01), delta_mag(832nm).",
        "companions_detected": "none reported (no stellar companions on ExoFOP; Tey+2023 report no companions)"
    },
    {
        "telescope": "SOAR (4.1 m)", "instrument": "HRCam", "type": "Speckle",
        "filter": "Cousins I (879 nm, 289 nm width)",
        "date_ut": "2021-10-01", "user": "ziegler", "exofop_tag": "100581",
        "exofop_contrast_summary": "delta 6.5 mag @ 1\" (ExoFOP table entry; the uploaded plot shows ~5.1 mag at 1.0\" and ~6.5 mag at ~3.15\" - the ExoFOP summary looks inconsistent with the plot; Tey+2023 text says sensitive to a 5-mag-fainter star at 1\")",
        "pixel_scale_arcsec": 0.01575, "psf_arcsec": 0.06364,
        "contrast_dmag": {"0.05": 1.0, "0.1": 2.05, "0.15": 2.7, "0.2": 4.4, "0.5": 5.3, "1.0": 5.07, "3.15": 6.55},
        "contrast_source": "Read by eye from TIC354944123I-cz20211001.pdf (points: 0.05,0.10,0.15,0.20,1.0,3.15 arcsec; 0.5\" value is the plotted broken-line fit). NO machine-readable SOAR contrast table is on ExoFOP.",
        "files": ["TIC354944123I-cz20211001.pdf"],
        "file_format": "PDF plot only (Delta magnitude I-band vs Delta arcsec, with ACF inset)",
        "companions_detected": "none reported"
    },
]

# ---------------- time series ----------------
notes_by_tag = {
    "97637": "KarenC (collins): defocused intended-full on .02 cut short by weather; hint of on-time ~1.2 ppt ingress, data clouded, unreliable at 1.4 ppt level (ExoFOP obsnotes). Excluded by Tey+2023/Parc+2026.",
    "97875": "Schwarz: 1.59 ppt transit of .01 (b) detected 39 min early (6.5 sigma vs then-ephemeris); Tc=2459385.5679; 22 px (8.58\") aperture; detrended vs Y(FITS). Parc label B1.",
    "98323": "Schwarz: 1.10 ppt transit of .01 (b) detected 23 min late (3.8 sigma vs then-ephemeris); Tc=2459407.7376; 16 px (6.24\"). Parc label B2.",
    "98329": "Schwarz: intended full on .02 (c) cut short (weather), effectively ingress-only; possible ingress (1.08 ppt, questionable, strong systematics). Excluded by Tey+2023/Parc+2026. Data window 2407.432-2407.605 BTJD (from measurements.tbl).",
    "98878": "Schwarz: .02 (c) full, 1.29 ppt, ~14.5 min late (1.1 sigma); Tc=2459439.6264 (notes file says 9349.6264 - a typo, data are at 2459439.5x); NEB check cleared field to 2.5'. Simultaneous with tag 98886 (two CTIO 1m). Parc label C1/C2.",
    "98886": "Schwarz: .02 (c) full on second CTIO 1m, 1.29 ppt, ~11.5 min late (0.9 sigma); Tc=2459439.6244 (notes typo 9349.6244). Parc label C1/C2.",
    "98887": "Schwarz: .01 (b) full from LCO-SSO, 1.26 ppt, ~13 min early (1.6 sigma); Tc=2459440.9662. Parc label B3.",
    "99618": "Gan: .01 (b) full, ~5 min early ~1.5 ppt on-target (obsnotes). Parc label B4.",
    "99695": "Collins: supplemental reduction (20 px aperture) of tag 99618, same night.",
    "99904": "Schwarz: .02 (c) full, 2.38 ppt, ~2.5 min early (0.19 sigma); Tc=2459471.6804; tracking imprecise. Parc label C3 (detrended vs sky/pixel).",
}
file_by_tag = {
    "97637": "TIC354944123-02_20210530_LCO-CTIO-1m0_ip_20px_bjd-flux-err-am.dat",
    "97875": "TIC354944123-01_20210619_LCO-SAAO-1m0_ip_bjd-flux-err-yfit_subset.dat",
    "98323": "TIC354944123-01_20210712_LCO-CTIO-1m0_ip_bjd-flux-err-am_subset.dat",
    "98329": "TIC354944123-02_20210711_LCO-SAAO-1m0_ip_measurements.tbl",
    "98878": "TIC354944123-02_20210813-a_LCO-CTIO-1m0_ip_bjd-flux-err-am.dat",
    "98886": "TIC354944123-02_20210813-b_LCO-CTIO-1m0_ip_bjd-flux-err-am.dat",
    "98887": "TIC354944123-01_20210814_LCO-SSO-1m0_ip_bjd-flux-err-air_subset.dat",
    "99618": "TIC354944123-01_20210831_LCO-CTIO-1m0_ip_20px_bjd-flux-err-xcen.dat",
    "99695": "TIC354944123-01_20210831_LCO-CTIO-1m0_ip_20px_bjd-flux-err-xcen.dat",
    "99904": "TIC354944123-02_20210914_LCO-CTIO-1m0_ip_bjd-flux-err-sky_subset.dat",
}
wmap = {w['file']: w for w in W}
toi_from_name = {"97637": "TOI 4342.02 (then TOI-18.02)", "97875": "TOI 4342.01 (then TOI-18.01)", "98323": "TOI 4342.01 (then TOI-18.01)", "98329": "TOI 4342.02 (then TOI-18.02)"}
tseries = []
for ts in T['time_series']:
    tag = ts['tstag']
    f = file_by_tag.get(tag)
    w = wmap.get(f, {})
    rec = {
        "date_ut": ts['tsdate'], "telescope": ts['tstel'], "camera": ts['tscam'], "filter": ts['tsfilt'],
        "toi": ts['tstoi'] or toi_from_name.get(tag, ''), "coverage_exofop": ts['tstcov'],
        "duration_min_exofop": ts['tsdur'], "n_images": ts['tsnum'], "aperture_px": ts['tspar'], "psf_arcsec": ts['tspsf'],
        "user": ts['tsuser'], "exofop_tag": tag, "exofop_note": ts['tsnotes'],
        "result_notes": notes_by_tag.get(tag, ''),
        "data_file_used_for_window": f,
        "window_btjd_first": w.get('btjd_first'), "window_btjd_last": w.get('btjd_last'), "window_span_h": w.get('span_h'),
        "n_points_in_file": w.get('n'), "rms_ppt_raw_normalized": w.get('rms_ppt'),
        "window_source": "actual BJD_TDB time stamps in the downloaded data file, BTJD = BJD_TDB - 2457000",
        "predicted_events_near_window": w.get('events', {}),
    }
    tseries.append(rec)

cand_offsets = {
    "description": "Nominal candidate ephemerides propagated to each SG1 window (no period uncertainty applied). A: P=18.9083, T0=1675.535 BTJD. B: P=1.95769, T0=1658.952 BTJD.",
    "A_18.9083d": "No SG1 window lies within 2 days of a nominal A transit (closest: 2021-09-14 window, nominal A transit at BTJD 2469.684, 48 h before the window; epoch 42 -> would need a period error of ~68 min/epoch). SG1 data do not test candidate A.",
    "A_half_period_9.454d": "No overlap either (closest 6.4 h after the 2021-08-14 SSO window).",
    "B_1.95769d": "Only the 2021-08-31 CTIO window (tag 99618/99695) reaches a nominal B transit: Tc=2457.6895 BTJD (epoch 408), 0.19 h after the window ends -> ~35% of the 1.26 h transit falls inside the last minutes of the night (after planet b egress). A crude linear-baseline in-vs-out test on those last ~26 points gives -0.93 +/- 0.30 ppt (i.e., brighter, not a dip; end-of-night systematics likely). Other windows are 2.4-22 h from nominal B events; with epochs 361-415 the B ephemeris uncertainty (sigma_P x ~400) could easily be hours, so this is not a constraint without a proper period error.",
}

# ---------------- spectroscopy ----------------
spectroscopy = []
for s in T['spectroscopy']:
    spectroscopy.append({"telescope": s['stel'], "instrument": s['sinst'], "resolution": s['ssres'], "wavelength": s['swave'],
                         "snr": s['ssnr1'], "snr_wavelength": s['ssnr2'], "date_ut": s['sdate'].strip(), "user": s['specuser'],
                         "group": s['sgroup'], "exofop_tag": s['stag'], "notes": s['snotes']})
spectroscopy_extra = {
    "chiron_stellar_params_brennom_2023-12-14": [
        {k: sp.get(k) for k in ('teff', 'logg', 'vsini', 'met', 'snr', 'fitq', 'stag')}
        for sp in T['stellar_parameters'] if sp.get('inst') == 'CHIRON'],
    "chiron_note": "5 CHIRON stellar-parameter fits (Teff 4027-4327 K, logg 4.65-4.71, vsini ~5.9-6.2 km/s [likely at/below CHIRON resolution floor], [M/H] -0.36 to -0.55). Files TIC354944123O-wb2021*.pdf (not downloaded). CHIRON 1D spectra FITS uploaded 2025-04-04 (not downloaded).",
    "tey2023_rv": "Tey+2023 fit the CHIRON RVs: K_b=19.4+5.1-5.0, K_c=14.4+-3.6 m/s -> 3-sigma upper limits 0.321 and 0.289 MJup (rules out stellar-mass companions only).",
    "precise_rv_elsewhere": "ESPRESSO (58 spectra, Apr 2022-Mar 2023) and NIRPS (Aug-Oct 2023) in Parc+2026; ESPRESSO data on DACE (doi 10.82180/dace-pwgeam18). Not on ExoFOP.",
}

# ---------------- TOI / CTOI ----------------
toi_rows = []
with open(os.path.join(BASE, 'raw', 'toi_rows_with_header.csv'), newline='') as fh:
    for r in csv.DictReader(fh):
        toi_rows.append({k: r[k] for k in ('TIC ID', 'TOI', 'TESS Disposition', 'TFOPWG Disposition', 'Epoch (BJD)', 'Epoch (BJD) err',
                                             'Period (days)', 'Period (days) err', 'Duration (hours)', 'Duration (hours) err',
                                             'Depth (ppm)', 'Depth (ppm) err', 'Planet Radius (R_Earth)', 'Planet SNR', 'Sectors',
                                             'Date TOI Alerted (UTC)', 'Date Modified', 'Comments')})

summary = {
    "target": {"tic": 354944123, "toi": 4342, "names": T['basic_info']['star_names'], "confirmed_planets": T['basic_info']['confirmed_planets'],
               "ra_deg": T['coordinates']['ra'], "dec_deg": T['coordinates']['dec'], "tic_contamination_ratio": T['basic_info']['tic_contamination_ratio'],
               "former_internal_name": "TOI-18.01 / TOI-18.02 (Evan Tey & Chelsea Huang custom QLP targets; renamed TOI-4342.01/.02 on 2021-07-28)"},
    "imaging": imaging,
    "tseries": tseries,
    "tseries_candidate_overlap": cand_offsets,
    "tseries_obsnotes_master": "ExoFOP observing notes (mlund, 2025-11-17): .01 and .02 both Master Disp VP, Phot Disp VPC, Spec Disp RVP; TFOPWG disposition CP (retired as VP 2023-03-03). 'No more SG1 observations needed.' Field cleared of NEBs to 2.5' (2021-08-13).",
    "spectroscopy": spectroscopy,
    "spectroscopy_extra": spectroscopy_extra,
    "toi_rows": toi_rows,
    "ctoi_rows": [],
    "toi_ctoi_check": {
        "toi_list_rows_total": 8151, "ctoi_list_rows_total": 5138,
        "entries_on_tic_354944123": "Only TOI 4342.01 (P=5.538329 d) and TOI 4342.02 (P=10.688608 d). Zero CTOIs on this TIC (ExoFOP 'CTOIs 0'). No TOI or CTOI within 300 arcsec of the target.",
        "periods_near_18p9_9p45_37p8_28p4_12p6": "none on this star",
        "periods_near_1p958_0p979_3p915": "none on this star",
        "spoc_tce_catalogs_checked": "MAST SPOC DV TCE stats for s0001-s0036, s0001-s0069, s0001-s0096 (multi-sector; only TCE-1 P=5.538 d and TCE-2 P=10.689 d), single-sector s13 (10.689 d), s27/s66/s93/s94 (no TCE), s67 (10.711 d), s101/s102/s103 (5.53-5.54 d). NO SPOC TCE at ~18.9 d or ~1.96 d (or aliases) in any run. Rows saved in lit/mast/tce_*.csv."
    },
    "planet_table_exofop": {
        "TOI 4342.01": {"epoch_bjd": "2458654.532426 +/- 0.0015552", "period_d": "5.538329 +/- 0.000038", "duration_h": "2.226 +/- 0.419", "depth_ppm": "1407 +/- 130", "rp_re": "2.27 +/- 1.12", "snr": 8.50, "source": "SPOC s13+s27 (TOI release 2021-07-28, params 2023-03-03)"},
        "TOI 4342.02": {"epoch_bjd": "2458659.348023 +/- 0.0024498", "period_d": "10.688608 +/- 0.000115", "duration_h": "2.667 +/- 0.346", "depth_ppm": "1501 +/- 162", "rp_re": "2.46 +/- 0.50", "snr": 8.19},
        "confirmed_b_exoplanet_archive_default(Parc+2026)": {"epoch_bjd": "2458654.53479 (+0.00084/-0.00093)", "period_d": "5.5382592 +/- 0.0000034", "rp_re": "2.329 +/- 0.086", "mass_me": "7.3 +/- 1.3", "inc_deg": "88.83 +/- 0.22", "b": "0.382 +/- 0.062", "rp_rs": "0.03571 +/- 0.00073", "dur_h": "2.182 +/- 0.024", "a_au": 0.0519},
        "confirmed_c_exoplanet_archive_default(Parc+2026)": {"epoch_bjd": "2458659.3486 +/- 0.0017", "period_d": "10.688662 +/- 0.00015 (Table 2; Table 8 of the paper gives +/- 0.000015)", "rp_re": "2.349 +/- 0.094", "mass_me": "4.8 +/- 1.4", "inc_deg": "89.61 +/- 0.26", "b": "0.20 +/- 0.12 (Table 8 gives 0.320)", "rp_rs": "0.03602 +/- 0.00092", "dur_h": "2.854 +/- 0.053", "a_au": 0.0804},
    },
    "stellar": {
        "tic_v8.2": {"teff": "3883 +/- 157", "logg": "4.653 +/- 0.011", "radius_rsun": "0.598 +/- 0.018", "mass_msun": "0.587 +/- 0.020", "density_gcc": "3.86 +/- 0.21", "lum_lsun": 0.0733, "dist_pc": "61.52 +/- 0.10", "tmag": 11.0318, "V": 12.669, "J": 9.832, "H": 9.179, "K": 9.018, "Gaia_G": 11.9822},
        "chiron_brennom": "Teff 4027-4327 K, logg 4.65-4.71, [M/H] -0.36..-0.55, vsini 5.9-6.2 km/s (5 epochs)",
        "pm_mas_yr": {"ra": 120.399, "dec": -91.4254},
        "stellar_companions": "none listed",
    },
    "files_downloaded": {
        "imaging_dir": os.path.join(BASE, 'imaging'),
        "imaging": sorted(os.listdir(os.path.join(BASE, 'imaging'))),
        "tseries_dir": os.path.join(BASE, 'tseries'),
        "tseries": sorted(os.listdir(os.path.join(BASE, 'tseries'))),
        "tseries_note": "Downloaded all .dat (BJD_TDB, rel_flux_T1_n, rel_flux_err_T1_n, detrend param), AIJ measurements .tbl, notes.txt, lightcurve.png and .plotcfg files. Skipped plate-solved FITS images, NEBcheck zips and other PNG diagnostics. The 0-byte NEB-table_5px.txt is empty on ExoFOP itself. The gan 20210831 measurements.tbl has BJD_TDB = 0 (not computed in that reduction); use the collins 20px .dat/.tbl (same night) or JD_UTC.",
        "raw_dir": os.path.join(BASE, 'raw'),
        "raw": sorted(os.listdir(os.path.join(BASE, 'raw'))),
        "sg1_windows_json": os.path.join(BASE, 'sg1_windows.json'),
    },
    "urls_checked": [
        "https://exofop.ipac.caltech.edu/tess/target.php?id=354944123",
        "https://exofop.ipac.caltech.edu/tess/target.php?id=354944123&json",
        "https://exofop.ipac.caltech.edu/tess/edit_obsnotes.php?id=354944123",
        "https://exofop.ipac.caltech.edu/tess/get_file.php?id=<fid> (fids 226014-226017, 232918 imaging; 51 light-curve fids listed in raw/ts_filelist.txt)",
        "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
        "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv",
        "https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/354944123/tces/",
        "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html",
        "https://archive.stsci.edu/missions/tess/catalogs/tce/<run>_dvr-tcestats.csv (12 runs, see lit/mast/)",
    ],
    "not_accessed_or_not_available": [
        "No SOAR machine-readable contrast curve exists on ExoFOP (PDF only).",
        "Gemini/Zorro 562 nm channel curve not uploaded to ExoFOP (only 832 nm).",
        "CHIRON spectra FITS and CHIRON stellar-parameter PDFs not downloaded (not needed).",
        "Plate-solved FITS images and NEBcheck zips skipped.",
    ],
    "fetched_utc": datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
}
json.dump(summary, open(os.path.join(BASE, 'exofop_summary.json'), 'w'), indent=1, default=str)
print('wrote', os.path.join(BASE, 'exofop_summary.json'), len(tseries), 'tseries')
