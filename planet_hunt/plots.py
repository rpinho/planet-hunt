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
