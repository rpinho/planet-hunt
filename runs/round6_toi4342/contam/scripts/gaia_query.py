import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from astroquery.utils.tap.core import TapPlus
out = sys.argv[1]
tap = TapPlus(url="https://gea.esac.esa.int/tap-server/tap", verbose=False)
sid = 6355915466181029376
def q(s):
    return tap.launch_job(s).get_results()
res = {}
t = q(f"SELECT * FROM gaiadr3.gaia_source WHERE source_id={sid}")
row = t[0]
keys = ['source_id','ra','dec','ref_epoch','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error','ruwe',
 'astrometric_excess_noise','astrometric_excess_noise_sig','ipd_frac_multi_peak','ipd_gof_harmonic_amplitude','ipd_gof_harmonic_phase',
 'ipd_frac_odd_win','non_single_star','duplicated_source','radial_velocity','radial_velocity_error','rv_amplitude_robust','rv_chisq_pvalue',
 'rv_renormalised_gof','rv_nb_transits','rv_method_used','phot_bp_rp_excess_factor','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','bp_rp',
 'astrometric_params_solved','visibility_periods_used','astrometric_n_good_obs_al','astrometric_chi2_al','phot_variable_flag','classprob_dsc_combmod_star','teff_gspphot']
g = {}
for k in keys:
    if k in t.colnames:
        v = row[k]
        try:
            if np.ma.is_masked(v): v=None
            elif isinstance(v,(bytes,)): v=v.decode()
            elif isinstance(v,(np.integer,)): v=int(v)
            elif isinstance(v,(np.floating,float)): v=float(v)
            elif isinstance(v,(np.bool_,bool)): v=bool(v)
            else: v=str(v)
        except Exception: v=str(v)
        g[k]=v
# corrected bp_rp excess C*
try:
    C=g['phot_bp_rp_excess_factor']; x=g['bp_rp']
    if x<0.5: f=1.154360+0.033772*x+0.032277*x**2
    elif x<4: f=1.162004+0.011464*x+0.049255*x**2-0.005879*x**3
    else: f=1.057572+0.140537*x
    g['C_star_corrected_excess']=C-f
except Exception as e: g['C_star_err']=str(e)
res['gaia_target']=g
nss={}
for tab in ['nss_two_body_orbit','nss_acceleration_astro','nss_non_linear_spectro','nss_vim_fl']:
    try:
        r=q(f"SELECT * FROM gaiadr3.{tab} WHERE source_id={sid}")
        nss[tab]=[{c:(None if np.ma.is_masked(r[0][c]) else str(r[0][c])) for c in r.colnames}] if len(r) else []
    except Exception as e: nss[tab]=f"error: {e}"
for tab in ['vari_summary']:
    try:
        r=q(f"SELECT * FROM gaiadr3.{tab} WHERE source_id={sid}")
        nss[tab]=[{c:str(r[0][c]) for c in r.colnames}] if len(r) else []
    except Exception as e: nss[tab]=f"error: {e}"
res['nss']=nss
# neighbours
ra0,de0=g['ra'],g['dec']
nb=q(f"""SELECT source_id,ra,dec,phot_g_mean_mag,bp_rp,parallax,parallax_error,pmra,pmdec,ruwe,
 DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{ra0},{de0}))*3600 AS sep
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra0},{de0},{150/3600.}))""")
nb.write(out.replace('.json','_neighbours.ecsv'),overwrite=True)
json.dump(res,open(out,'w'),indent=1,default=str)
print(json.dumps(res['gaia_target'],indent=1,default=str))
print({k:(len(v) if isinstance(v,list) else v) for k,v in nss.items()})
print(len(nb))
