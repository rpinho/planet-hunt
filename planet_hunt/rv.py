"""Method B: periodic signals in archival HARPS radial velocities.

Data: the HARPS RVBank (Trifonov et al. 2020, A&A 636, A74; CDS catalogue
J/A+A/636/A74): SERVAL re-reduction of every public HARPS spectrum up to
2018, ~212,000 RVs of ~2,900 stars, corrected for nightly zero-points, plus
activity indicators (CRX, dLW, H-alpha, Na D1/D2 from SERVAL; FWHM, contrast,
BIS from the DRS cross-correlation function).

The recipe (standard "pre-whitening"):
  1. per star, average the RVs of each night; treat data before and after the
     2015 fibre upgrade as two instruments (independent zero-points, jitter);
  2. generalised Lomb-Scargle (GLS) periodogram of the residuals; if the best
     peak has an analytic false-alarm probability (Baluev 2008) below the
     threshold, add a Keplerian at that period and refit everything jointly
     (offsets, a linear trend, all signals, per-instrument jitter);
  3. repeat until the next peak is not significant;
  4. vet every signal: activity indicators with the same period (or 2x,
     1/2x), RV-indicator correlations, 1-day / 1-year / lunar aliases of
     stronger signals, periods at 1 yr or the lunar month, and whether the
     signal is present with the same amplitude in both halves of the data
     (planets are coherent for years; spots and plage are not).
"""
from __future__ import annotations

import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from . import config

RVB_DIR = config.CACHE / "rvbank"
CDS = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/636/A74/"
FIBRE_UPGRADE_BJD = 2457174.5           # HARPS fibre exchange, end of May 2015
LA_SILLA_NOON = 0.196                   # JD fraction of local noon (night boundary)
INDICATORS = ["RHK", "FWHMDRS", "BIS", "Contrast", "dLW", "CRX", "Halpha", "NaD1", "NaD2"]
SIDEREAL_DAY, YEAR, LUNAR = 0.99727, 365.25, 29.53
G, MSUN, MEARTH, DAY = 6.674e-11, 1.989e30, 5.972e24, 86400.0


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
def _download(name: str):
    RVB_DIR.mkdir(parents=True, exist_ok=True)
    path = RVB_DIR / name
    if not path.exists():
        print(f"[rvbank] downloading {name} from CDS ...", flush=True)
        urllib.request.urlretrieve(CDS + name, path)
    return path


def load_rvbank() -> pd.DataFrame:
    """All RVBank rows (cached as a pickle after the first parse, ~50 MB gz)."""
    pkl = RVB_DIR / "rvbank.pkl"
    if pkl.exists():
        return pd.read_pickle(pkl)
    import warnings
    from astropy.table import Table
    readme, data = _download("ReadMe"), _download("rvbank.dat.gz")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = Table.read(data, format="ascii.cds", readme=str(readme)).to_pandas()
    df["Name"] = df["Name"].str.strip()
    df.to_pickle(pkl)
    return df


RHK_DIR = config.CACHE / "rvbank_rhk"


def load_rhk() -> pd.DataFrame:
    """Per-spectrum Ca II H&K emission R'HK for the RVBank spectra
    (Perdelwitz et al. 2024, A&A 683, A125; CDS J/A+A/683/A125), the
    classic chromospheric activity index. Only the three needed columns are
    fetched, through the VizieR TAP service. Columns: Name, BJD, RHK."""
    pkl = RHK_DIR / "rhk.pkl"
    if pkl.exists():
        return pd.read_pickle(pkl)
    RHK_DIR.mkdir(parents=True, exist_ok=True)
    csv = RHK_DIR / "rhk_tap.csv"
    if not csv.exists():
        print("[rhk] downloading R'HK (Target, BJD, RHKp) from VizieR TAP ...", flush=True)
        q = 'SELECT Target, BJD, RHKp FROM "J/A+A/683/A125/table4"'
        url = ("https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync?request=doQuery&lang=adql"
               f"&format=csv&query={urllib.parse.quote(q)}")
        urllib.request.urlretrieve(url, csv.with_suffix(".part"))
        csv.with_suffix(".part").rename(csv)
    df = pd.read_csv(csv).rename(columns={"Target": "Name", "RHKp": "RHK"})
    df["Name"] = df["Name"].astype(str).str.strip()
    df.to_pickle(pkl)
    return df


def add_rhk(df: pd.DataFrame) -> pd.DataFrame:
    """Attach R'HK to RVBank rows: same star, nearest time within 0.005 d
    (the two catalogues' time stamps differ by ~1 minute)."""
    try:
        rhk = load_rhk()
    except Exception as exc:
        print(f"[rhk] not available: {exc!r}", flush=True)
        return df.assign(RHK=np.nan)
    left = df.reset_index().sort_values("BJD")
    right = rhk.dropna(subset=["BJD"]).sort_values("BJD")
    m = pd.merge_asof(left, right, on="BJD", by="Name", tolerance=0.005, direction="nearest")
    return m.set_index("index").sort_index().rename_axis(None)


def load_star_list() -> pd.DataFrame:
    """RVBank target list with decimal coordinates."""
    path = _download("list.dat")
    rows = []
    for line in path.read_text().splitlines():
        if len(line) < 100:
            continue

        def num(a, b):
            s = line[a:b].strip()
            try:
                return float(s)
            except ValueError:
                return np.nan
        ra = 15 * (num(16, 18) + num(19, 21) / 60 + num(22, 29) / 3600)
        dec = num(31, 33) + num(34, 36) / 60 + num(37, 43) / 3600
        if line[30] == "-":
            dec = -dec
        rows.append({"Name": line[0:15].strip(), "ra": ra, "dec": dec,
                     "SpType": line[44:64].strip(), "Teff": num(65, 71), "Bmag": num(72, 78),
                     "Vmag": num(79, 85), "Gmag": num(86, 92), "Dist": num(93, 100)})
    return pd.DataFrame(rows)


def nightly(star: pd.DataFrame, rv_col: str = "DRVmlcnzp") -> pd.DataFrame:
    """Weighted nightly means of the RVs (and plain means of the indicators).
    Instrument 0 = before the 2015 fibre upgrade, 1 = after."""
    d = star[(star["Flag"] == 0) & np.isfinite(star[rv_col]) & np.isfinite(star["e_" + rv_col])
             & (star["e_" + rv_col] > 0) & (star["e_" + rv_col] < 50)].copy()
    if d.empty:
        return pd.DataFrame()
    d["inst"] = (d["BJD"] > FIBRE_UPGRADE_BJD).astype(int)
    d["night"] = np.floor(d["BJD"] - LA_SILLA_NOON).astype(int)
    for c in ("FWHMDRS", "Contrast", "BIS"):
        d.loc[d[c] < -9e6, c] = np.nan
    d["w"] = 1 / d["e_" + rv_col] ** 2
    out = []
    for (night, inst), g in d.groupby(["night", "inst"]):
        w = g["w"].to_numpy()
        row = {"bjd": float(np.sum(w * g["BJD"]) / w.sum()), "inst": int(inst),
               "rv": float(np.sum(w * g[rv_col]) / w.sum()), "err": float(1 / np.sqrt(w.sum())),
               "n_exp": len(g)}
        for c in INDICATORS:
            v = g[c].to_numpy(float)
            row[c] = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
        out.append(row)
    n = pd.DataFrame(out).sort_values("bjd").reset_index(drop=True)
    # 5-sigma clip of gross outliers around the per-instrument median (bad nights)
    keep = np.ones(len(n), bool)
    for i in (0, 1):
        m = n["inst"] == i
        if m.sum() > 5:
            r = n.loc[m, "rv"] - n.loc[m, "rv"].median()
            mad = 1.4826 * np.median(np.abs(r))
            keep[m.to_numpy()] = np.abs(r.to_numpy()) < max(10 * mad, 30)
    return n[keep].reset_index(drop=True)


# --------------------------------------------------------------------------
# Keplerian model
# --------------------------------------------------------------------------
def kepler_E(M, e, n_iter: int = 30):
    E = M + e * np.sin(M) * (1 + e * np.cos(M))
    for _ in range(n_iter):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-10:
            break
    return E


def keplerian(t, P, K, sec, ses, lam, tref):
    """RV of one Keplerian. sec/ses = sqrt(e) cos/sin(omega); lam = mean
    longitude (M + omega) at tref, radians."""
    e = sec ** 2 + ses ** 2
    w = np.arctan2(ses, sec) if e > 0 else 0.0
    M = 2 * np.pi * (t - tref) / P + lam - w
    if e < 1e-6:
        return K * np.cos(M + w)
    E = kepler_E(np.mod(M, 2 * np.pi), e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    return K * (np.cos(nu + w) + e * np.cos(w))


class RVModel:
    """offsets (per instrument) + linear trend + N Keplerians, fit by weighted
    least squares with per-instrument jitter added in quadrature."""

    def __init__(self, t, rv, err, inst):
        self.t, self.rv, self.err, self.inst = map(np.asarray, (t, rv, err, inst))
        self.tref = float(np.round(np.median(self.t)))
        self.insts = sorted(set(self.inst.tolist()))
        self.jit = {i: 0.0 for i in self.insts}
        self.signals: list[dict] = []        # each: P, K, sec, ses, lam
        self.eccentric = False
        self.lin = np.zeros(len(self.insts) + 1)   # offsets..., slope

    # -- pieces --------------------------------------------------------
    def sigma(self):
        return np.sqrt(self.err ** 2 + np.array([self.jit[i] for i in self.inst]) ** 2)

    def kep_sum(self, signals=None, skip: int | None = None):
        out = np.zeros_like(self.t, dtype=float)
        for k, s in enumerate(signals if signals is not None else self.signals):
            if k != skip:
                out += keplerian(self.t, s["P"], s["K"], s["sec"], s["ses"], s["lam"], self.tref)
        return out

    def design(self):
        cols = [(self.inst == i).astype(float) for i in self.insts]
        cols.append((self.t - self.tref) / 365.25)
        return np.column_stack(cols)

    def _linear(self, y, sig):
        A = self.design() / sig[:, None]
        coef, *_ = np.linalg.lstsq(A, y / sig, rcond=None)
        return coef

    def residuals(self, signals=None):
        y = self.rv - self.kep_sum(signals)
        sig = self.sigma()
        coef = self._linear(y, sig)
        return y - self.design() @ coef, coef

    # -- fitting ---------------------------------------------------------
    def _pack(self):
        p, lo, hi = [], [], []
        T = float(self.t.max() - self.t.min())
        for s in self.signals:
            # the GLS peak is ~P^2/T wide: long periods need a wide bound
            dp = max(0.02 * s["P"], 0.5 * s["P"] ** 2 / T)
            p += [s["P"], s["K"], s["lam"]]
            lo += [max(s["P"] - dp, 1.0), 0.0, -np.inf]
            hi += [s["P"] + dp, np.inf, np.inf]
            if self.eccentric:
                p += [s["sec"], s["ses"]]
                lo += [-0.9, -0.9]
                hi += [0.9, 0.9]
        return np.array(p), np.array(lo), np.array(hi)

    def _unpack(self, p):
        out, j = [], 0
        for _ in self.signals:
            s = {"P": p[j], "K": p[j + 1], "lam": p[j + 2], "sec": 0.0, "ses": 0.0}
            j += 3
            if self.eccentric:
                s["sec"], s["ses"] = p[j], p[j + 1]
                j += 2
            out.append(s)
        return out

    def fit(self, n_jitter_iter: int = 3):
        for _ in range(n_jitter_iter):
            sig = self.sigma()
            if self.signals:
                p0, lo, hi = self._pack()
                p0 = np.clip(p0, lo + 1e-9, hi - 1e-9)

                def fun(p):
                    r, _ = self.residuals(self._unpack(p))
                    return r / sig
                sol = least_squares(fun, p0, bounds=(lo, hi), x_scale="jac", max_nfev=400)
                self.signals = self._unpack(sol.x)
            r, self.lin = self.residuals()
            # jitter: make reduced chi^2 ~ 1 per instrument
            for i in self.insts:
                m = self.inst == i
                if m.sum() < 5:
                    continue
                s2 = np.mean(r[m] ** 2 - self.err[m] ** 2)
                self.jit[i] = float(np.sqrt(max(s2, 0.0)))
        for s in self.signals:
            if s["K"] < 0:
                s["K"], s["lam"] = -s["K"], s["lam"] + np.pi
        return self

    def loglike(self):
        r, _ = self.residuals()
        sig = self.sigma()
        return float(-0.5 * np.sum((r / sig) ** 2 + np.log(2 * np.pi * sig ** 2)))

    def n_params(self):
        return len(self.insts) + 1 + len(self.insts) + len(self.signals) * (5 if self.eccentric else 3)

    def bic(self):
        return -2 * self.loglike() + self.n_params() * np.log(len(self.t))

    def param_errors(self):
        """1-sigma errors on P and K from the Jacobian at the solution."""
        if not self.signals:
            return []
        sig = self.sigma()
        p0, lo, hi = self._pack()

        def fun(p):
            r, _ = self.residuals(self._unpack(p))
            return r / sig
        J = np.empty((len(self.t), len(p0)))
        f0 = fun(p0)
        for j in range(len(p0)):
            h = 1e-6 * max(abs(p0[j]), 1e-3)
            pp = p0.copy(); pp[j] += h
            J[:, j] = (fun(pp) - f0) / h
        try:
            cov = np.linalg.inv(J.T @ J)
        except np.linalg.LinAlgError:
            cov = np.full((len(p0), len(p0)), np.nan)
        step = 5 if self.eccentric else 3
        return [{"P_err": float(np.sqrt(abs(cov[k * step, k * step]))),
                 "K_err": float(np.sqrt(abs(cov[k * step + 1, k * step + 1])))}
                for k in range(len(self.signals))]

    def copy(self):
        m = RVModel(self.t, self.rv, self.err, self.inst)
        m.tref, m.jit, m.eccentric = self.tref, dict(self.jit), self.eccentric
        m.signals = [dict(s) for s in self.signals]
        m.lin = self.lin.copy()
        return m


# --------------------------------------------------------------------------
# Periodograms
# --------------------------------------------------------------------------
def gls(t, y, dy, pmin: float = 1.1, pmax: float | None = None, ofac: float = 10):
    """Generalised Lomb-Scargle (floating mean). Returns freq, power, ls."""
    from astropy.timeseries import LombScargle
    T = t.max() - t.min()
    pmax = pmax or T
    fmin, fmax = 1 / pmax, 1 / pmin
    freq = np.arange(fmin, fmax, 1 / (ofac * T))
    ls = LombScargle(t, y, dy, fit_mean=True, center_data=True)
    return freq, ls.power(freq), ls


def fap(ls, power, fmax) -> float:
    try:
        return float(ls.false_alarm_probability(power, maximum_frequency=fmax, method="baluev"))
    except Exception:
        return float("nan")


def window_peaks(t, n: int = 5):
    """Strongest peaks of the spectral window (sampling periodicities)."""
    from astropy.timeseries import LombScargle
    T = t.max() - t.min()
    freq = np.arange(1 / (3 * T), 1.2, 1 / (10 * T))
    w = LombScargle(t, np.ones_like(t), fit_mean=False, center_data=False).power(freq)
    idx = np.argsort(w)[::-1]
    out = []
    for i in idx:
        if all(abs(freq[i] - f) > 2 / T for f, _ in out):
            out.append((float(freq[i]), float(w[i])))
        if len(out) >= n:
            break
    return out


# --------------------------------------------------------------------------
# Pre-whitening search
# --------------------------------------------------------------------------
def prewhiten(n: pd.DataFrame, fap_max: float = 1e-3, max_signals: int = 6,
              pmin: float = 1.1, pmax_factor: float = 1.0) -> tuple[RVModel, list[dict]]:
    t, rv, err, inst = (n[c].to_numpy(float) for c in ("bjd", "rv", "err", "inst"))
    model = RVModel(t, rv, err, inst.astype(int)).fit()
    T = t.max() - t.min()
    steps = []
    for k in range(max_signals):
        r, _ = model.residuals()
        freq, power, ls = gls(t, r, model.sigma(), pmin=pmin, pmax=pmax_factor * T)
        i = int(np.argmax(power))
        f_best, p_best = float(freq[i]), float(power[i])
        fa = fap(ls, p_best, freq.max())
        steps.append({"k": k + 1, "period": 1 / f_best, "power": p_best, "fap": fa,
                      "freq": freq, "pow": power})
        if not (fa < fap_max):
            break
        # The GLS peak can be a 1-day (or 1-year) alias of the true period:
        # with nightly sampling, f and |f - 1/day| fit almost equally well.
        # Fit each alias as the new signal and keep the best likelihood.
        cands = {"peak": f_best}
        for name, fa in (("1-day", 1 / SIDEREAL_DAY), ("1-year", 1 / YEAR)):
            for sgn, lab in ((-1, "-"), (1, "+")):
                fc = abs(f_best + sgn * fa)
                if 1 / (pmax_factor * T) < fc < 1 / pmin:
                    cands[f"{name} alias ({lab})"] = fc
        trials = {}
        for name, fc in cands.items():
            trials[name] = _add_signal(model, r, fc)
        ll = {name: m.loglike() for name, m in trials.items()}
        best = max(ll, key=ll.get)
        model = try_eccentric(trials[best])
        steps[-1]["alias_choice"] = best
        steps[-1]["alias_dlnL"] = {name: float(ll[best] - v) for name, v in ll.items() if name != best}
        steps[-1]["alias_periods"] = {name: float(1 / fc) for name, fc in cands.items()}
    return model, steps


def _add_signal(model: RVModel, r: np.ndarray, f: float) -> RVModel:
    """Copy of `model` with one more circular signal at frequency f, refit."""
    t = model.t
    X = np.column_stack([np.cos(2 * np.pi * f * (t - model.tref)),
                         np.sin(2 * np.pi * f * (t - model.tref))])
    sig = model.sigma()
    c, *_ = np.linalg.lstsq(X / sig[:, None], (r - np.average(r, weights=1 / sig ** 2)) / sig,
                            rcond=None)
    trial = model.copy()
    trial.signals.append({"P": 1 / f, "K": max(float(np.hypot(*c)), 0.1), "sec": 0.0, "ses": 0.0,
                          "lam": float(np.arctan2(-c[1], c[0]))})
    return trial.fit()


def try_eccentric(model: RVModel) -> RVModel:
    """Refit with free eccentricities; keep it only if BIC improves by > 10."""
    if not model.signals:
        return model
    ecc = model.copy()
    ecc.eccentric = True
    for s in ecc.signals:
        if s["sec"] == 0 and s["ses"] == 0:
            s["sec"], s["ses"] = 0.1, 0.1
    try:
        ecc.fit()
    except Exception:
        return model
    if model.eccentric:          # already eccentric: keep whichever fits better
        return ecc if ecc.bic() < model.bic() else model
    return ecc if ecc.bic() < model.bic() - 10 else model


def min_mass_earth(K, P, mstar, e=0.0):
    """Minimum mass (Earth masses), planet mass << star mass."""
    return float(K * (P * DAY / (2 * np.pi * G)) ** (1 / 3) * (mstar * MSUN) ** (2 / 3)
                 * np.sqrt(1 - e ** 2) / MEARTH)


# Dwarf mass from Teff (Pecaut & Mamajek 2013 table, rounded). Crude but fine
# for minimum masses; giants are flagged separately by their absolute mag.
_TEFF = [2800, 3100, 3400, 3600, 3850, 4100, 4400, 4700, 5000, 5300, 5550, 5770, 6000, 6300, 6600]
_MASS = [0.10, 0.16, 0.30, 0.45, 0.58, 0.64, 0.70, 0.75, 0.80, 0.87, 0.93, 1.00, 1.08, 1.20, 1.35]


def mass_from_teff(teff: float) -> float:
    return float(np.interp(teff, _TEFF, _MASS)) if np.isfinite(teff) else 1.0


# --------------------------------------------------------------------------
# Vetting
# --------------------------------------------------------------------------
def near(f1, f2, T, k: float = 1.5) -> bool:
    return abs(f1 - f2) < k / T


def alias_flags(P, earlier_periods, T):
    f = 1 / P
    flags = []
    for name, P0 in (("1 yr", YEAR), ("1/2 yr", YEAR / 2), ("lunar month", LUNAR)):
        if abs(P / P0 - 1) < 0.04:
            flags.append(f"period at the {name} (sampling/moonlight systematic?)")
    if P > 0.6 * T:
        flags.append("period longer than ~half the data span (trend / cycle?)")
    for Pj in earlier_periods:
        fj = 1 / Pj
        for name, fa in (("1-day", 1 / SIDEREAL_DAY), ("1-year", 1 / YEAR), ("lunar", 1 / LUNAR)):
            if near(f, abs(fj - fa), T) or near(f, fj + fa, T):
                flags.append(f"{name} alias of the stronger {Pj:.2f} d signal")
        for h in (2, 3, 0.5, 1 / 3):
            if abs(P / (Pj * h) - 1) < 0.02:
                flags.append(f"harmonic ({h:g}x) of the {Pj:.2f} d signal")
    return flags


def indicator_periodograms(n: pd.DataFrame, T: float) -> dict:
    """GLS of every activity indicator (fibre-upgrade jump removed), once per star."""
    t = n["bjd"].to_numpy(float)
    out = {}
    for c in INDICATORS:
        y = n[c].to_numpy(float).copy()
        ok = np.isfinite(y)
        if ok.sum() < 20 or np.nanstd(y[ok]) == 0:
            continue
        for i in (0, 1):
            m = ok & (n["inst"].to_numpy() == i)
            if m.any():
                y[m] -= np.median(y[m])
        freq, power, ls = gls(t[ok], y[ok], None, pmin=1.1, pmax=T)
        thr = float(ls.false_alarm_level(0.01, maximum_frequency=freq.max(), method="baluev"))
        out[c] = {"ok": ok, "y": y, "freq": freq, "power": power, "thr": thr,
                  "best_period": float(1 / freq[np.argmax(power)]),
                  "best_fap": fap(ls, float(power.max()), freq.max())}
    return out


def activity_check(pg: dict, P: float, rv_signal: np.ndarray, T: float) -> dict:
    """Does any activity indicator vary with the signal's period (or 2x, 1/2x)?
    Does the RV signal correlate with an indicator? `pg` from
    indicator_periodograms()."""
    from scipy import stats
    out, flags = {}, []
    for c, d in pg.items():
        freq, power, thr, ok = d["freq"], d["power"], d["thr"], d["ok"]
        hits = []
        for h in (1.0, 2.0, 0.5):
            sel = np.abs(freq - 1 / (P * h)) < 1.0 / T
            if sel.any() and power[sel].max() > thr:
                hits.append("same period" if h == 1 else f"{h:g}x")
        r, p = stats.pearsonr(rv_signal[ok], d["y"][ok])
        out[c] = {"peak_hits": hits, "r": float(r), "p": float(p),
                  "best_period": d["best_period"], "best_fap": d["best_fap"]}
        if hits:
            flags.append(f"{c} has a significant periodogram peak at {'/'.join(hits)} of the signal period")
        if p < 1e-3 and abs(r) > 0.3:
            flags.append(f"RV signal correlates with {c} (r = {r:+.2f})")
    return {"indicators": out, "flags": flags}


def coherence_check(model: RVModel, k: int) -> dict:
    """Fit signal k's amplitude separately in the first and second half of
    the nights (everything else fixed). A planet keeps the same K."""
    t = model.t
    order = np.argsort(t)
    half = np.zeros(len(t), bool)
    half[order[len(t) // 2:]] = True
    s = model.signals[k]
    others = model.kep_sum(skip=k)
    r_lin = model.rv - others - model.design() @ model.lin
    sig = model.sigma()
    res = []
    for m in (~half, half):
        ph = 2 * np.pi * (t[m] - model.tref) / s["P"]
        X = np.column_stack([np.cos(ph), np.sin(ph), np.ones(m.sum())])
        A = X / sig[m][:, None]
        c, *_ = np.linalg.lstsq(A, r_lin[m] / sig[m], rcond=None)
        cov = np.linalg.pinv(A.T @ A)
        K = float(np.hypot(c[0], c[1]))
        Ke = float(np.sqrt((c[0] ** 2 * cov[0, 0] + c[1] ** 2 * cov[1, 1]) / max(K ** 2, 1e-12)))
        res.append((K, Ke))
    (K1, e1), (K2, e2) = res
    diff = abs(K1 - K2) / np.hypot(e1, e2)
    ok = diff < 3 and K1 > 2 * e1 and K2 > 2 * e2
    return {"K_first": K1, "K_first_err": e1, "K_second": K2, "K_second_err": e2,
            "diff_sigma": float(diff), "coherent": bool(ok)}


# --------------------------------------------------------------------------
# Known planets (NASA Exoplanet Archive)
# --------------------------------------------------------------------------
def load_known_planets(max_age_days: float = 3.0) -> pd.DataFrame:
    """pscomppars (one row per confirmed planet) + every period ever listed in
    the `ps` table (alternative published solutions), with host identifiers."""
    path = config.CATALOG_DIR / "exoarchive_rvhosts.csv"
    if not (path.exists() and time.time() - path.stat().st_mtime < max_age_days * 86400):
        q1 = ("select pl_name,hostname,hd_name,hip_name,ra,dec,pl_orbper,pl_bmasse,"
              "pl_rvamp,discoverymethod,disc_year from pscomppars")
        q2 = "select pl_name,hostname,pl_orbper from ps where pl_orbper is not null"
        dfs = []
        for q in (q1, q2):
            url = f"{config.EXOARCHIVE_TAP}?query={urllib.parse.quote(q)}&format=csv"
            with urllib.request.urlopen(url, timeout=300) as r:
                dfs.append(pd.read_csv(r))
        comp, ps = dfs
        extra = ps.merge(comp.drop(columns=["pl_orbper"]).drop_duplicates("pl_name"),
                         on=["pl_name", "hostname"], how="left")
        allp = pd.concat([comp.assign(source="pscomppars"), extra.assign(source="ps")],
                         ignore_index=True)
        allp.to_csv(path, index=False)
    allp = pd.read_csv(path)
    eu = load_exoplanet_eu()
    return pd.concat([allp, eu], ignore_index=True) if len(eu) else allp


EU_URL = "https://exoplanet.eu/catalog/csv/"


def load_exoplanet_eu(max_age_days: float = 7.0) -> pd.DataFrame:
    """The Encyclopaedia of Exoplanetary Systems (exoplanet.eu), confirmed
    planets. It lists some RV planets the NASA archive does not (e.g.
    HD 157172 b). Returned in the NASA column names used here."""
    path = config.CATALOG_DIR / "exoplanet_eu.csv"
    if not (path.exists() and time.time() - path.stat().st_mtime < max_age_days * 86400):
        try:
            tmp = path.with_suffix(".part")
            urllib.request.urlretrieve(EU_URL, tmp)
            tmp.rename(path)
        except Exception as exc:
            print(f"[exoplanet.eu] download failed: {exc!r}", flush=True)
    if not path.exists():
        return pd.DataFrame()
    e = pd.read_csv(path, low_memory=False, on_bad_lines="skip")

    def pick(alt, prefix):
        for a in str(alt).split(","):
            a = a.strip()
            if a.upper().startswith(prefix):
                return a
        return np.nan
    alt = e["star_alternate_names"].fillna("")
    return pd.DataFrame({
        "pl_name": e["name"], "hostname": e["star_name"],
        "hd_name": [pick(a, "HD ") for a in alt], "hip_name": [pick(a, "HIP ") for a in alt],
        "ra": pd.to_numeric(e["ra"], errors="coerce"), "dec": pd.to_numeric(e["dec"], errors="coerce"),
        "pl_orbper": pd.to_numeric(e["orbital_period"], errors="coerce"),
        "pl_bmasse": pd.to_numeric(e["mass_sini"], errors="coerce") * 317.8,
        "pl_rvamp": np.nan, "discoverymethod": e["detection_type"],
        "disc_year": pd.to_numeric(e["discovered"], errors="coerce"), "source": "exoplanet.eu"})


def norm_name(s) -> str:
    s = str(s).upper().strip()
    for pre in ("NAME ", "V* ", "* ", "** "):
        if s.startswith(pre):
            s = s[len(pre):]
    s = s.replace(" ", "").replace("_", "")
    for a, b in (("GLIESE", "GJ"), ("GL", "GJ")):
        if s.startswith(a) and not s.startswith("GJ"):
            s = b + s[len(a):]
    return s


def simbad_ids(names: list[str]) -> dict[str, list[str]]:
    """All SIMBAD identifiers of each star (HD, HIP, GJ, proper names...),
    cached in cache/catalogs/simbad_ids.csv. Needed because the archive names
    hosts differently (GJ 699 = Barnard's star, GJ 17 = zet Tuc, ...)."""
    path = config.CATALOG_DIR / "simbad_ids.csv"
    have = pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["name", "ids", "ra", "dec"])
    todo = [n for n in names if n not in set(have["name"])]
    if todo:
        from astroquery.simbad import Simbad
        sb = Simbad()
        sb.add_votable_fields("ids")
        rows = []
        for i in range(0, len(todo), 100):
            chunk = todo[i:i + 100]
            try:
                r = sb.query_objects(chunk).to_pandas()
            except Exception as exc:
                print(f"[simbad] {exc!r}", flush=True)
                continue
            for k, nm in enumerate(chunk):
                m = r[r["object_number_id"] == k + 1] if "object_number_id" in r else r.iloc[k:k + 1]
                if len(m):
                    rows.append({"name": nm, "ids": str(m["ids"].iloc[0]) + "|" + str(m["main_id"].iloc[0]),
                                 "ra": m["ra"].iloc[0], "dec": m["dec"].iloc[0]})
                else:
                    rows.append({"name": nm, "ids": nm, "ra": np.nan, "dec": np.nan})
        have = pd.concat([have, pd.DataFrame(rows)], ignore_index=True)
        have.to_csv(path, index=False)
    return {r.name: [norm_name(x) for x in str(r.ids).split("|")] + [norm_name(r.name)]
            for r in have.itertuples() if r.name in set(names)}


def host_planets(star: pd.Series, known: pd.DataFrame, ids: list[str] | None = None,
                 radius_arcsec: float = 60.0) -> pd.DataFrame:
    """Known planets whose host matches this RVBank star by any SIMBAD
    identifier or by position."""
    idset = set(ids or []) | {norm_name(star["Name"])}
    import re
    strip_a = lambda x: re.sub(r"(?<=\d)A$", "", x)    # "HD 215456 A" -> "HD215456"
    names = pd.concat([known[c].map(norm_name).map(strip_a) for c in ("hostname", "hd_name", "hip_name")],
                      axis=1)
    by_name = names.isin(idset).any(axis=1)
    dra = (known["ra"] - star["ra"]) * np.cos(np.radians(star["dec"]))
    sep = np.hypot(dra, known["dec"] - star["dec"]) * 3600
    by_pos = sep < radius_arcsec
    hosts = set(known.loc[by_name | by_pos, "hostname"].dropna())
    return known[known["hostname"].isin(hosts)]


def match_known(P: float, planets: pd.DataFrame, tol: float = 0.05) -> str | None:
    for r in planets.dropna(subset=["pl_orbper"]).itertuples():
        for h, note in ((1, ""), (2, " (at 2x)"), (0.5, " (at 1/2)")):
            if abs(P / (h * r.pl_orbper) - 1) < tol:
                return f"{r.pl_name} (P = {r.pl_orbper:.3f} d){note}"
        # 1-day aliases of a known period
        fk = 1 / r.pl_orbper
        for fa in (1 / SIDEREAL_DAY,):
            if abs(1 / P - abs(fk - fa)) < 0.002 or abs(1 / P - (fk + fa)) < 0.002:
                return f"1-day alias of {r.pl_name} (P = {r.pl_orbper:.3f} d)"
    return None


# --------------------------------------------------------------------------
# Calibration: how many known planets on the searched stars came back?
# --------------------------------------------------------------------------
def recovery_table(results: pd.DataFrame, stars_meta: list[dict], known: pd.DataFrame,
                   ids: dict, star_list: pd.DataFrame, k_min: float = 2.0) -> pd.DataFrame:
    """One row per known planet (NASA pscomppars) on a searched star with
    period 1.1 d .. data span and listed RV semi-amplitude >= k_min m/s:
    was it recovered (a signal within 5% of its period, or 2x / 1/2)?"""
    lst = star_list.set_index("Name")
    comp = known[known["source"] == "pscomppars"]
    rows = []
    for m in stars_meta:
        nm, T = m["name"], m.get("T_days", np.nan)
        if nm not in lst.index or not np.isfinite(T):
            continue
        st = lst.loc[nm]
        pl = host_planets(pd.Series({"Name": nm, "ra": st["ra"], "dec": st["dec"]}), comp, ids.get(nm))
        sig = results[results["star"] == nm]
        for p in pl.drop_duplicates("pl_name").itertuples():
            if not (np.isfinite(p.pl_orbper) and 1.1 < p.pl_orbper < T):
                continue
            if not (np.isfinite(p.pl_rvamp) and p.pl_rvamp >= k_min):
                continue
            tol = np.clip(3 * sig["period_err"].fillna(0) / sig["period"], 0.05, 0.25)
            hit = sig[(np.abs(sig["period"] / p.pl_orbper - 1) < tol)]
            rows.append({"star": nm, "planet": p.pl_name, "P": p.pl_orbper, "K_lit": p.pl_rvamp,
                         "n_nights": m.get("n_nights"), "recovered": len(hit) > 0,
                         "K_found": float(hit["K"].iloc[0]) if len(hit) else np.nan})
    return pd.DataFrame(rows)
