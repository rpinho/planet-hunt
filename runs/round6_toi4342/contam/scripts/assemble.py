import sys, json
base = sys.argv[1]
g = json.load(open(f'{base}/gaia_raw.json')); nb = json.load(open(f'{base}/neighbours.json')); ar = json.load(open(f'{base}/archival_raw.json'))
gt = g['gaia_target']; gt['source_id'] = 6355915466181029376
gt['notes'] = ('RUWE 0.98 (<1.4 single-star-like). AEN 0.16 mas with sig 39.8 is formally significant but typical of the G~12 '
  'calibration floor in DR3 (RUWE is the better binarity indicator at this brightness); ipd_frac_multi_peak=0 and '
  'ipd_gof_harmonic_amplitude=0.011 show no resolved-double signature; non_single_star=0, duplicated_source=False. '
  'RV -4.43 +/- 0.36 km/s from 25 transits; rv_amplitude_robust 5.5 km/s with rv_chisq_pvalue 0.23 (>0.01) -> no significant RV variability '
  '(Gaia RV epoch precision for an M0 at G_RVS~10.5 is a few km/s, so this only excludes large-K SB1s). Corrected BP/RP excess C*=0.004 (no blend). '
  'phot_variable_flag=VARIABLE: G std 6.0 mmag over 51 epochs (range 25 mmag) - consistent with spot modulation on an M dwarf. '
  'TIC v8.2 lists Gaia DR2 6355915466180482944 (DR2 id; DR3 id differs) and an epoch-2000 position.')
mov = [n for n in nb['neighbours'] if n['pmra'] is not None and abs(n['pmra']-gt['pmra']) < 20 and abs(n['pmdec']-gt['pmdec']) < 20]
for n in mov: n['note'] = ('similar proper motion to target but parallax %.2f mas vs target 16.25 mas -> ~%.0f pc vs 61.5 pc; not a bound companion '
                           '(tangential velocities ~%.0f vs 44 km/s); chance pm alignment' % (n['parallax'], 1000/n['parallax'], 4.74*((n['pmra']**2+n['pmdec']**2)**0.5)/n['parallax']))
m1150 = [dict(source_id=n['source_id'], sep_arcsec=n['sep_arcsec'], pa_deg=n['pa_deg'], G=n['G'], Tmag_est=n['Tmag_est'], dTmag=n['dTmag'], depth_needed=round(n['depth_needed_1150ppm'], 4)) for n in nb['neighbours'] if n['mimic_1150']]
m284 = [dict(source_id=n['source_id'], sep_arcsec=n['sep_arcsec'], pa_deg=n['pa_deg'], G=n['G'], Tmag_est=n['Tmag_est'], dTmag=n['dTmag'], depth_needed=round(n['depth_needed_284ppm'], 4)) for n in nb['neighbours'] if n['mimic_284']]
arch = []
for a in ar:
    if 'error' in a: arch.append(a); continue
    keep = {k: a[k] for k in ['survey', 'band', 'epoch', 'pixel_scale_arcsec', 'fwhm_arcsec', 'offset_of_target_from_now_arcsec', 'offset_over_fwhm',
            'target_measured_minus_predicted_arcsec', 'astrom_rms_arcsec', 'n_astrom_stars', 'current_position_brightness_frac_of_target_peak',
            'resid_flux_at_current_pos_sigma', 'target_mag_in_proxy', 'limit_mag_at_current_position', 'm5_uncapped', 'delta_mag_limit',
            'delta_T_limit_same_colour', 'delta_T_limit_typical_field_colour', 'band_minus_T_target', 'band_minus_T_field_median',
            'usable_for_background_test', 'note', 'injected_stamps_mag_snr']}
    if not a['usable_for_background_test']:
        keep['limit_mag_at_current_position_formal'] = keep.pop('limit_mag_at_current_position'); keep['limit_mag_at_current_position'] = None
        keep['delta_mag_limit_formal'] = keep.pop('delta_mag_limit'); keep['delta_mag_limit'] = None
        for k in ['delta_T_limit_same_colour', 'delta_T_limit_typical_field_colour']: keep[k] = None
    arch.append(keep)
sj = next(a for a in arch if a['survey'].startswith('DSS1'))
faint_rec = max(m for m, snr in sj['injected_stamps_mag_snr'] if snr > 5)
dBP_verified = round(faint_rec - sj['target_mag_in_proxy'], 2)
colour_shift = round(sj['band_minus_T_target'] - sj['band_minus_T_field_median'], 2)
out = dict(
  target=dict(toi='TOI-4342', tic=354944123, gaia_dr3=6355915466181029376, Tmag_tic=nb['tic']['Tmag']),
  gaia_target=gt, nss=dict(nss_two_body_orbit=g['nss']['nss_two_body_orbit'], nss_acceleration_astro=g['nss']['nss_acceleration_astro'],
      nss_non_linear_spectro=g['nss']['nss_non_linear_spectro'], nss_vim_fl=g['nss']['nss_vim_fl'],
      summary='No entries in any Gaia DR3 NSS table (two-body orbit, astrometric acceleration, non-linear spectro, VIM). Target is in vari_summary (in_vari_classification_result=True).'),
  tic=nb['tic'], tic_contratio=nb['tic']['contratio'], tic_numcont=int(nb['tic']['numcont']),
  tmag_target_stassun_from_gaia=round(nb['T0'], 3),
  neighbours_note='All Gaia DR3 sources within 150" (57). dTmag relative to Stassun-estimated target T=%.3f. depth_needed = depth*(1+f)/f with f=10^(-0.4 dT); assumes the full neighbour flux is in the aperture (upper bound on dilution, so a conservative mimic flag). Nearest Gaia source is 16.8" away; none within 16".' % nb['T0'],
  neighbours=nb['neighbours'],
  mimic_capable={'1150ppm': m1150, '284ppm': m284,
     'summary': '%d of 57 neighbours could in principle produce 1150 ppm (needed eclipse <100%%); %d could produce 284 ppm. Brightest/closest: 60.9" G=14.96 (needs 2.6%%/0.7%%), 94.5" G=14.97 (2.7%%/0.7%%), 140.5" G=14.04 (1.2%%/0.3%%), 37.6" G=17.65 (22%%/5.4%%), 53.0" G=17.17 (17%%/4.2%%), 18.9" G=18.44 (59%%/15%%), 24.1" G=19.23 (74%%/18%%). All are resolved by Gaia and TESS centroid/difference-image or ground-based photometry can test them individually.' % (len(m1150), len(m284))},
  archival=arch,
  archival_summary=dict(
    surveys_used=['DSS1 SERC-J 1976.66 (oldest; DSS2-blue not available for this field)', 'DSS2 AAO-SES R 1992.64', 'DSS2 SERC-I 1995.67', '2MASS JHKs 2000.54', 'SkyMapper DR4 r/i 2014.74-2014.76'],
    not_available=['Pan-STARRS (no coverage south of Dec -30)', 'DSS1 red / POSS-I (northern; no southern DSS1 red plate for this field: STScI and ESO DSS servers return errors)', 'SDSS/DECaLS/DECaPS (no coverage)'],
    not_attempted=['VISTA VHS / unWISE (VHS ~2010-2015 and WISE ~2010 are too recent: displacement 1.7-2.5 arcsec < FWHM; WISE FWHM 6 arcsec)'],
    epoch_now=2026.78, total_pm_mas_yr=round((gt['pmra']**2+gt['pmdec']**2)**0.5, 2),
    displacement_oldest_to_now_arcsec=sj['offset_of_target_from_now_arcsec'],
    oldest_plate_fwhm_arcsec=sj['fwhm_arcsec'],
    dBP_limit_verified_by_real_star_injection=dBP_verified,
    dBP_limit_adopted=sj['delta_mag_limit'], dBP_limit_extrapolated=round(sj['m5_uncapped']-sj['target_mag_in_proxy'], 2),
    dT_limit_same_colour=[dBP_verified, sj['delta_mag_limit']],
    dT_limit_typical_background_colour=[round(dBP_verified+colour_shift, 2), round(sj['delta_mag_limit']+colour_shift, 2)],
    colour_note='Background stars are typically bluer than the M0 target (field median BP-T=%.2f vs target %.2f), so a given BP contrast corresponds to a ~%.2f mag larger T contrast. Same-colour value is the conservative case.' % (sj['band_minus_T_field_median'], sj['band_minus_T_target'], colour_shift),
    caveats=['Photographic density is non-linear; injection adds background-subtracted real-star density stamps linearly. At the current position the plate density is only 1.6% of the target peak (near sky), so the plate is in its usable regime there.',
             'Stamps are calibrated with Gaia BP (proxy for B_J); BP of G~19-20 stars is biased bright, which makes the limit conservative.',
             'The limit applies within ~1 FWHM (3.4") of the 2026 position; Gaia DR3 independently shows no source within 16.8".',
             'All later epochs (1992, 1995, 2000, 2014) have shift/FWHM < 1.5 and/or the current position inside the saturated halo, so they give no independent background constraint.']),
)
needed_1150 = 7.4; needed_284 = 8.9
out['verdict'] = dict(
  bound_companions=('No Gaia evidence of a bound companion: RUWE 0.98, ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.011, non_single_star 0, no NSS solution, '
    'duplicated_source False, C* 0.004, RV scatter consistent with noise (rv_chisq_pvalue 0.23). This excludes companions that perturb the astrometry '
    '(roughly ~0.1-few arcsec, moderate contrast; RUWE<1.4 at 61 pc limits massive companions at ~1-10 au) but not close (<~0.1") faint or near-equal-period-blind configurations; high-resolution imaging (ExoFOP) and RV remain the tests for bound EBs. '
    'No co-moving Gaia source within 150" (the one similar-pm star at 95" has parallax 5.7 mas, not 16.2).'),
  background_hidden_now=('On the 1976.66 SERC-J plate the target sat %.1f" (%.1f x FWHM) from its 2026 position; the 2026 position is at sky level there and shows no source. '
    'Real field stars injected at that spot are recovered at >5 sigma down to BP=%.2f (dBP>=%.2f; adopted 5-sigma dBP~%.2f). In TESS band this excludes background stars with dT < ~%.1f (same colour as target, conservative) '
    'or dT < ~%.1f (typical bluer background colour). 1150 ppm needs dT<=7.4 (100%% eclipse): a same-colour hidden star would need a >=%d%% eclipse - largely ruled out, fully ruled out for typical colours. '
    '284 ppm needs dT<=8.9: background stars with dT~7.3-8.9 (needed eclipse 24-100%%) are NOT excluded by archival imaging; a bright-enough hidden background star for 284 ppm is constrained only partially.' % (
      sj['offset_of_target_from_now_arcsec'], sj['offset_over_fwhm'], faint_rec, dBP_verified, sj['delta_mag_limit'], dBP_verified, dBP_verified+colour_shift,
      int(round(100*1150e-6*10**(0.4*dBP_verified))))))
json.dump(out, open(f'{base}/contamination.json', 'w'), indent=1)
print(json.dumps(out['verdict'], indent=1)); print(json.dumps(out['archival_summary'], indent=1)); print(json.dumps(mov, indent=1))
