"""Diagnostic figures. One function -> one 5-panel figure per signal:

  (a) normalised PDCSAP light curve + removed trend, transit times marked
  (b) BLS periodogram, peak and its 1/2, 2x aliases marked
  (c) phase-folded flattened light curve, binned points + box model
  (d) odd vs even transits side by side (EB check)
  (e) phase 0.5 window (secondary-eclipse check)
"""
from __future__ import annotations

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE, GREY, INK, MUTED = "#2a78d6", "#eb6834", "#a3a29b", "#0b0b0b", "#52514e"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": "#ecebe7",
    "grid.linewidth": 0.6, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "axes.titlelocation": "left", "axes.formatter.useoffset": False,
})


def _phase(t, period, t0):
    return (t - t0 + 0.5 * period) % period - 0.5 * period


def _binned(x, y, nbins, lo, hi):
    edges = np.linspace(lo, hi, nbins + 1)
    idx = np.digitize(x, edges) - 1
    ok = (idx >= 0) & (idx < nbins)
    cnt = np.bincount(idx[ok], minlength=nbins)
    s = np.bincount(idx[ok], y[ok], minlength=nbins)
    mid = 0.5 * (edges[1:] + edges[:-1])
    good = cnt > 2
    return mid[good], s[good] / cnt[good]


def candidate_figure(lc, bls: dict, st: dict, title: str):
    P, t0, D = st["period"], st["t0"], st["duration_h"] / 24
    depth = st["depth_ppm"] * 1e-6
    t, flux, trend, flat = lc.time, lc.flux, lc.trend, lc.flat
    ph = _phase(t, P, t0)

    fig = plt.figure(figsize=(9.2, 7.4), dpi=95, layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1])
    fig.suptitle(title, x=0.01, ha="left", fontsize=10.5, fontweight="bold", color=INK)

    # (a) raw light curve
    ax = fig.add_subplot(gs[0, :])
    ax.plot(t, flux, ".", ms=1, color=GREY, rasterized=True)
    ax.plot(t, trend, "-", lw=1, color=ORANGE, label="trend removed")
    n0 = int(np.ceil((t.min() - t0) / P)); n1 = int(np.floor((t.max() - t0) / P))
    for n in range(n0, n1 + 1):
        ax.axvline(t0 + n * P, color=BLUE, lw=0.8, alpha=0.6)
    lo, hi = np.percentile(flux, [0.5, 99.5]); pad = (hi - lo) * 0.3
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_title("(a) PDCSAP flux (blue ticks = predicted transits)")
    ax.set_xlabel("Time (BTJD = BJD − 2457000)"); ax.set_ylabel("Norm. flux")
    ax.legend(loc="upper right", frameon=False, fontsize=8)

    # (b) periodogram
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(bls["periods"], bls["power"], "-", lw=0.7, color=MUTED)
    ax.axvline(P, color=BLUE, lw=1.6, alpha=0.5)
    for r in (0.5, 2):
        if bls["periods"].min() < P * r < bls["periods"].max():
            ax.axvline(P * r, color=BLUE, lw=0.8, ls=":", alpha=0.7)
    ax.set_xscale("log")
    ax.set_title(f"(b) BLS periodogram  P={P:.4f} d  SDE={st['sde']:.1f}")
    ax.set_xlabel("Period (d)"); ax.set_ylabel("BLS power")

    # (c) phase fold
    ax = fig.add_subplot(gs[1, 1])
    w = max(3 * D, 0.25)
    sel = np.abs(ph) < w
    x = ph[sel] * 24
    ax.plot(x, flat[sel], ".", ms=1.5, color=GREY, rasterized=True)
    xb, yb = _binned(x, flat[sel], 60, -w * 24, w * 24)
    ax.plot(xb, yb, "o", ms=3.5, color=BLUE, label="binned")
    level = np.median(flat[~sel]) if (~sel).any() else 1.0
    bx = np.array([-w, -D / 2, -D / 2, D / 2, D / 2, w]) * 24
    by = level - np.array([0, 0, depth, depth, 0, 0])
    ax.plot(bx, by, "-", lw=1.5, color=ORANGE, label="box model")
    _ylim_around(ax, yb, depth, level)
    ax.set_title(f"(c) Folded  depth={st['depth_ppm']:.0f} ppm  dur={st['duration_h']:.2f} h  SNR={st['snr']:.1f}")
    ax.set_xlabel("Hours from mid-transit"); ax.set_ylabel("Flattened flux")
    ax.legend(loc="lower right", frameon=False, fontsize=8)

    # (d) odd vs even
    ax = fig.add_subplot(gs[2, 0])
    epoch = np.round((t - t0) / P).astype(int)
    off = 2.5 * w * 24
    allb = []
    for parity, color, shift, name in ((1, BLUE, -off / 2, "odd"), (0, ORANGE, off / 2, "even")):
        s2 = sel & (epoch % 2 == parity)
        if s2.sum() < 5:
            continue
        xx = ph[s2] * 24
        ax.plot(xx + shift, flat[s2], ".", ms=1.2, color=GREY, rasterized=True)
        xb2, yb2 = _binned(xx, flat[s2], 40, -w * 24, w * 24)
        ax.plot(xb2 + shift, yb2, "o", ms=3, color=color, label=name)
        allb.append(yb2)
        d = st["depth_odd_ppm" if parity else "depth_even_ppm"] * 1e-6
        ax.hlines(level - d, shift - D * 12, shift + D * 12, color=color, lw=1.5)
    if allb:
        _ylim_around(ax, np.concatenate(allb), depth, level)
    ax.set_xticks([])
    ax.set_title(f"(d) Odd | even transits  Δ={st['oe_sigma']:.1f}σ")
    ax.legend(loc="lower center", frameon=False, fontsize=8, ncol=2)

    # (e) secondary eclipse window
    ax = fig.add_subplot(gs[2, 1])
    ph2 = _phase(t, P, t0 + 0.5 * P)
    s3 = np.abs(ph2) < w
    ax.plot(ph2[s3] * 24, flat[s3], ".", ms=1.5, color=GREY, rasterized=True)
    xb3, yb3 = _binned(ph2[s3] * 24, flat[s3], 60, -w * 24, w * 24)
    ax.plot(xb3, yb3, "o", ms=3.5, color=BLUE)
    ax.axvspan(-D * 12, D * 12, color=ORANGE, alpha=0.08)
    _ylim_around(ax, yb3, depth, level)
    ax.set_title(f"(e) Phase 0.5 (secondary?)  {st['depth_sec_ppm']:.0f} ppm, {st['sec_sigma']:.1f}σ")
    ax.set_xlabel("Hours from phase 0.5")
    return fig


def _ylim_around(ax, yb, depth, level):
    if len(yb) == 0:
        return
    lo = min(np.min(yb), level - depth)
    hi = max(np.max(yb), level)
    span = max(hi - lo, 3e-4)
    ax.set_ylim(lo - 0.5 * span, hi + 0.5 * span)


def fig_to_jpeg(fig, quality: int = 72) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="jpeg", pil_kwargs={"quality": quality, "optimize": True})
    plt.close(fig)
    return buf.getvalue()


# --------------------------------------------------------------------------
# Method A (stacked FFI transits) and method B (RV) figures
# --------------------------------------------------------------------------
def stacked_transit_figure(z, row: dict, title: str):
    """(a) season-summed BLS periodogram, (b) phase-folded transit from all
    sectors (binned), (c) depth per observing season, (d) phase 0.5."""
    P, t0, D = row["period"], row["t0"], row["duration_h"] / 24
    t, y, dy = z["t"], z["flat"].astype(float), z["dy"].astype(float)
    depth = row["depth_ppm"] * 1e-6
    ph = _phase(t, P, t0)
    fig = plt.figure(figsize=(9.2, 6.4), dpi=95, layout="constrained")
    gs = fig.add_gridspec(2, 2)
    fig.suptitle(title, x=0.01, ha="left", fontsize=10.5, fontweight="bold", color=INK)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(z["periods"], z["power"], "-", lw=0.6, color=MUTED)
    ax.axvline(P, color=BLUE, lw=1.6, alpha=0.4)
    ax.set_xscale("log")
    ax.set_title(f"(a) BLS power summed over seasons  SDE={row['sde']:.1f}")
    ax.set_xlabel("Period (d)"); ax.set_ylabel("Σ log-likelihood")

    ax = fig.add_subplot(gs[0, 1])
    w = max(3 * D, 0.15)
    sel = np.abs(ph) < w
    x = ph[sel] * 24
    ax.plot(x, y[sel], ".", ms=0.8, color=GREY, alpha=0.4, rasterized=True)
    nb = int(np.clip(2 * w / (D / 6), 30, 90))
    edges = np.linspace(-w * 24, w * 24, nb + 1)
    idx = np.digitize(x, edges) - 1
    wt = 1 / dy[sel] ** 2
    mids, ys, es = [], [], []
    for i in range(nb):
        m = idx == i
        if m.sum() > 3:
            mids.append(0.5 * (edges[i] + edges[i + 1]))
            ys.append(np.sum(wt[m] * y[sel][m]) / wt[m].sum())
            es.append(1 / np.sqrt(wt[m].sum()))
    mids, ys, es = map(np.array, (mids, ys, es))
    ax.errorbar(mids, ys, es, fmt="o", ms=3.5, color=BLUE, lw=0.8, label="binned, all sectors")
    bx = np.array([-w, -D / 2, -D / 2, D / 2, D / 2, w]) * 24
    ax.plot(bx, 1 - np.array([0, 0, depth, depth, 0, 0]), "-", lw=1.4, color=ORANGE, label="box model")
    if len(ys):
        lo, hi = min(ys.min(), 1 - depth), max(ys.max(), 1)
        span = max(hi - lo, 1e-3)
        ax.set_ylim(lo - 0.6 * span, hi + 0.6 * span)
    ax.set_title(f"(b) {row['n_transits']} transits folded  P={P:.4f} d  depth={row['depth_ppm']:.0f} ppm"
                 f"  {row['duration_h']:.2f} h  SNR={row['snr']:.0f}")
    ax.set_xlabel("Hours from mid-transit"); ax.set_ylabel("Flattened flux")
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)

    ax = fig.add_subplot(gs[1, 0])
    seas = row.get("season_depths") or []
    if seas:
        tm = [s["t_mid"] for s in seas]
        ax.errorbar(tm, [s["depth_ppm"] for s in seas], [s["err_ppm"] for s in seas], fmt="o",
                    color=BLUE, ms=5)
        ax.axhline(row["depth_ppm"], color=ORANGE, lw=1)
        ax.axhline(0, color=GREY, lw=0.8)
    ax.set_title("(c) Depth in each observing season (a planet's is constant)")
    ax.set_xlabel("Time (BTJD)"); ax.set_ylabel("Depth (ppm)")

    ax = fig.add_subplot(gs[1, 1])
    ph2 = _phase(t, P, t0 + 0.5 * P)
    s3 = np.abs(ph2) < w
    xb3, yb3 = _binned(ph2[s3] * 24, y[s3], 50, -w * 24, w * 24)
    ax.plot(ph2[s3] * 24, y[s3], ".", ms=0.8, color=GREY, alpha=0.4, rasterized=True)
    ax.plot(xb3, yb3, "o", ms=3.5, color=BLUE)
    ax.axvspan(-D * 12, D * 12, color=ORANGE, alpha=0.08)
    _ylim_around(ax, yb3, depth, 1.0)
    ax.set_title(f"(d) Phase 0.5 (secondary eclipse?)  {row['sec_sigma']:.1f}σ")
    ax.set_xlabel("Hours from phase 0.5")
    return fig


def rv_system_figure(z, title: str, labels: list[str] | None = None, max_panels: int = 6):
    """GLS periodogram before each pre-whitening step + each signal phase-folded
    with the other signals, offsets and trend removed."""
    import json
    from .rv import keplerian
    sig = json.loads(str(z["signals"]))
    t, rv, err, inst = z["t"], z["rv"], z["sigma"], z["inst"]
    tref, lin = float(z["tref"]), z["lin"]
    insts = sorted(set(inst.tolist()))
    base = sum(lin[k] * (inst == i) for k, i in enumerate(insts)) + lin[-1] * (t - tref) / 365.25
    n = min(len(sig), max_panels)
    steps = sorted(k for k in z.files if k.startswith("pg") and k != "pg_freq")
    fig = plt.figure(figsize=(9.2, 2.0 + 1.9 * max(n, 1)), dpi=95, layout="constrained")
    gs = fig.add_gridspec(max(n, 1), 2, width_ratios=[1.25, 1])
    fig.suptitle(title, x=0.01, ha="left", fontsize=10.5, fontweight="bold", color=INK)
    f = z["pg_freq"]
    for k in range(max(n, 1)):
        ax = fig.add_subplot(gs[k, 0])
        key = f"pg{k + 1}"
        if key in z.files:
            ax.plot(1 / f, z[key], "-", lw=0.6, color=MUTED)
            if k < len(sig):
                ax.axvline(sig[k]["P"], color=BLUE, lw=1.5, alpha=0.4)
            ax.set_xscale("log")
        ax.set_ylabel("GLS power")
        ax.set_title(f"({chr(97 + 2 * k)}) residual periodogram, step {k + 1}")
        if k == n - 1:
            ax.set_xlabel("Period (d)")
        if k >= len(sig):
            continue
        s = sig[k]
        others = sum(keplerian(t, o["P"], o["K"], o["sec"], o["ses"], o["lam"], tref)
                     for j, o in enumerate(sig) if j != k) if len(sig) > 1 else 0
        r = rv - base - others
        phs = ((t - tref) / s["P"]) % 1
        ax = fig.add_subplot(gs[k, 1])
        ax.errorbar(phs, r, err, fmt=".", ms=3, color=GREY, lw=0.4, alpha=0.7)
        xb, yb = _binned(phs, r, 10, 0, 1)
        ax.plot(xb, yb, "o", ms=5, color=BLUE)
        pp = np.linspace(0, 1, 200)
        tt = tref + pp * s["P"]
        ax.plot(pp, keplerian(tt, s["P"], s["K"], s["sec"], s["ses"], s["lam"], tref), "-",
                color=ORANGE, lw=1.6)
        lim = max(3.5 * s["K"], np.percentile(np.abs(r), 90))
        ax.set_ylim(-lim, lim)
        lab = labels[k] if labels and k < len(labels) else ""
        ax.set_title(f"({chr(98 + 2 * k)}) P = {s['P']:.3f} d, K = {s['K']:.2f} m/s {lab}")
        ax.set_ylabel("RV (m/s)")
        if k == n - 1:
            ax.set_xlabel("Orbital phase")
    return fig
