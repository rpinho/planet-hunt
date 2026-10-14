"""Round 8A (P1 list): candidate definitions, fixed BEFORE any new data was looked at.

Ephemerides are the discovery values (round 5 / round 4 / round 6) and are the
ones the hold-out test is run against. They are not re-fitted on the new data
before the hold-out depth is measured.

disc_sectors: sectors (and product) the discovery search used. Everything else
is hold-out ("new") data:
  * round-5 hosts search (TOI-1846, -4342, -6008, -249, -4405): SPOC 2-min only,
    S1-S107. Hold-out = FFI-only sectors (no 2-min that sector) and any sector > 107.
    FFI products of sectors that DO have 2-min share photons and are used only as
    product checks, never as hold-out.
  * round-5 faint M dwarf (TIC 350522259): TESS-SPOC FFI sectors searched in round 5.
  * round-4 newest-sector (TIC 11528365): 2-min S107 only.
"""
from __future__ import annotations

CANDS = {
    "toi1846": dict(
        name="TOI-1846 0.944 d", tic=198385543, ra=257.9549583333333, dec=54.52784722222222,
        tmag=11.8467, teff=3512.0, rstar=0.411452, mstar=0.404348,
        P=0.944059523509414, t0=1765.128039741291, dur_h=0.588616713186432, depth_ppm=257.07,
        disc="r5hosts", disc_max_sector=107, disc_product="SPOC-2min",
        known=[dict(name="TOI-1846 b", P=3.93067462792567, t0=2884.349803000223, dur_h=1.43051415609098)],
        eumemic_claim=None),
    "toi4342_d2": dict(
        name="TOI-4342 1.958 d", tic=354944123, ra=324.3894166666666, dec=-77.97914722222222,
        tmag=11.0318, teff=3883.0, rstar=0.598175, mstar=0.586597,
        # round-6 refined ephemeris (2-min fit with rho prior); round-5 discovery P/T0 differ by < 1 min over the baseline
        P=1.9576890755606433, t0=1658.9524010321534, dur_h=1.256898305349431, depth_ppm=284.2,
        disc="r5hosts", disc_max_sector=107, disc_product="SPOC-2min",
        known=[dict(name="TOI-4342 b", P=5.538329, t0=1654.532426, dur_h=2.23),
               dict(name="TOI-4342 c", P=10.688608, t0=1659.348023, dur_h=2.67),
               dict(name="d18 candidate", P=18.908226, t0=2904.5666, dur_h=2.44)],
        eumemic_claim="in Eumemic's public signal table (not novel; credit Eumemic)"),
    "toi6008": dict(
        name="TOI-6008 7.220 d", tic=286201103, ra=310.44874999999996, dec=49.6469,
        tmag=12.2936, teff=3098.0, rstar=0.235048, mstar=0.204861,
        P=7.220152258409202, t0=1717.7952704824131, dur_h=1.107227782969045, depth_ppm=1297.5,
        disc="r5hosts", disc_max_sector=107, disc_product="SPOC-2min",
        known=[dict(name="TOI-6008 b", P=0.857372697894355, t0=1711.3729039998725, dur_h=0.547859807639078)],
        eumemic_claim=None),
    "tic350522259": dict(
        name="TIC 350522259.01 0.739 d", tic=350522259, ra=None, dec=None,
        tmag=14.4605, teff=3230.0, rstar=0.296684, mstar=0.272838,
        P=0.7387734710264852, t0=1354.462780213149, dur_h=0.834, depth_ppm=1127.8,
        disc="r5faint", disc_max_sector=95, disc_product="TESS-SPOC",
        known=[], eumemic_claim=None),
    "tic11528365": dict(
        name="TIC 11528365 8.17 d", tic=11528365, ra=18.6674679621019, dec=-23.9610814839411,
        tmag=10.7034, teff=5746.0, rstar=0.841172, mstar=1.027,
        P=8.171535, t0=4267.422972, dur_h=1.944, depth_ppm=814.0,
        disc="r4new", disc_max_sector=107, disc_product="SPOC-2min", disc_sectors=[107],
        known=[], eumemic_claim=None),
    "toi249": dict(
        name="TOI-249 19.09 d", tic=179985715, ra=14.080041666666665, dec=-38.94923333333333,
        tmag=11.2738, teff=4157.0, rstar=0.603418, mstar=0.63,
        P=19.09284888946109, t0=1357.2946081212626, dur_h=1.1621379988086376, depth_ppm=1120.6,
        disc="r5hosts", disc_max_sector=107, disc_product="SPOC-2min",
        known=[dict(name="TOI-249 b", P=6.61531875224234, t0=1359.5415150001645, dur_h=2.33458443381212)],
        eumemic_claim="on Eumemic's counted list (do not claim)"),
    "toi4405": dict(
        name="TOI-4405 3.601 d", tic=177254988, ra=101.59533333333331, dec=-75.44412777777778,
        tmag=12.3813, teff=4361.0, rstar=0.687025, mstar=0.68,
        P=3.6008190429212275, t0=2038.7766323619169, dur_h=1.7619336436846291, depth_ppm=345.5,
        disc="r5hosts", disc_max_sector=107, disc_product="SPOC-2min",
        known=[dict(name="TOI-4405 b", P=5.92142305135093, t0=2039.6418280000798, dur_h=2.26805215626293)],
        eumemic_claim=None),
}
