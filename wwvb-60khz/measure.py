#!/usr/bin/env python3
"""Measure the WWVB grid: NIST's 60 kHz time code, one bit per second.

WWVB (Fort Collins, Colorado) transmits a 60 kHz carrier whose POWER
drops ~17 dB at the top of every second and restores after a duration
that IS the data:
    0.2 s low  ->  bit 0
    0.5 s low  ->  bit 1
    0.8 s low  ->  frame marker
Markers land on seconds 9, 19, 29, 39, 49, 59 (plus second 0), so the
marker comb ticks every 10 s and a full BCD date/time frame takes 60 s.

This script derives, from a raw IQ capture near 60 kHz:
  1. the 1 Hz second comb (from envelope-spectrum harmonics)
  2. the three pulse-width classes (matched-template per second)
  3. the 10 s marker comb and its coherence

Usage: python measure.py --iq capture.cs16 --fs 250000
"""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--secs", type=float, default=178.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    del raw

    # narrowband the carrier at the dial (+-30 Hz) and take the power
    # envelope at 100 Hz — the grid lives entirely in carrier POWER
    d1 = int(a.fs / 2500)
    z = lfilter(firwin(255, 1000, fs=a.fs), 1.0, x)[::d1]
    fs1 = a.fs / d1
    z = lfilter(firwin(301, 30, fs=fs1), 1.0, z)
    d2 = int(fs1 / 100)
    fs2 = fs1 / d2                       # 100 Hz envelope rate
    z = z[::d2][int(4 * fs2):]           # drop filter settle
    pw = np.abs(z) ** 2
    k = int(0.03 * fs2)
    pws = np.convolve(pw, np.ones(k) / k, "same")

    # --- 1) the second comb: envelope spectrum, harmonics 1..8
    e = pw - pw.mean()
    Nfft = 1 << 22
    E = np.abs(np.fft.rfft(e, Nfft))
    fax = np.fft.rfftfreq(Nfft, 1 / fs2)
    m1 = (fax > 0.97) & (fax < 1.03)
    f1 = fax[m1][np.argmax(E[m1])]
    floor = np.median(E[(fax > 2.2) & (fax < 4.8)])
    line_db = 20 * np.log10(E[m1].max() / floor)
    ests = []
    for h in range(1, 9):
        mh = (fax > h * f1 - 0.02) & (fax < h * f1 + 0.02)
        ests.append(fax[mh][np.argmax(E[mh])] / h)
    f_sec = np.mean(ests)
    print(f"1 Hz second comb: line {line_db:+.1f} dB above floor, "
          f"harmonic-avg {f_sec:.6f} Hz -> period {1/f_sec:.6f} s")

    # --- 2) find the second boundary by folding at 1 s and sliding a
    # bit-0 template (low first 0.2 s, high after)
    L = int(round(fs2))
    m = len(pws) // L
    prof = pws[: m * L].reshape(m, L).mean(0)
    tpl0 = np.ones(L)
    tpl0[: int(0.2 * L)] = 0.02
    cc = [np.dot(np.roll(prof, -s), tpl0) for s in range(L)]
    s0 = int(np.argmax(cc))

    hi = np.percentile(pws, 85)
    lo = np.percentile(pws, 8)
    r = lo / hi
    print(f"carrier power drop (envelope p85 -> p8): "
          f"{10*np.log10(hi/lo):.1f} dB")

    # --- 3) matched-template classification of every second
    T = {}
    for lab, wd in (("0", 0.2), ("1", 0.5), ("M", 0.8)):
        t = np.ones(L)
        t[: int(wd * L)] = r
        T[lab] = t / t.mean()
    mid = np.sqrt(hi * lo)
    kinds, margins, restore = [], [], []
    i = s0
    while i + L <= len(pws):
        seg = pws[i:i + L] / max(pws[i:i + L].mean(), 1e-9)
        scores = {lab: -np.mean((seg - t) ** 2) for lab, t in T.items()}
        order = sorted(scores, key=scores.get, reverse=True)
        kinds.append(order[0])
        margins.append(scores[order[0]] - scores[order[1]])
        # restore time: first sustained crossing back above mid-level
        above = pws[i:i + L] > mid
        rr = np.nan
        for j in range(int(0.10 * L), L - 5):
            if above[j:j + 5].all():
                rr = j / L
                break
        restore.append(rr)
        i += L
    kinds = np.array(kinds)
    margins = np.array(margins)
    restore = np.array(restore)
    gate = margins > np.percentile(margins, 30)   # confident 70 %
    print(f"confident seconds: {gate.sum()}/{len(kinds)}  "
          f"{dict(Counter(kinds[gate]))}")
    for lab, target in (("0", 0.2), ("1", 0.5), ("M", 0.8)):
        sel = restore[gate & (kinds == lab) & ~np.isnan(restore)]
        sel = sel[np.abs(sel - target) < 0.15]
        if len(sel):
            print(f"  width '{lab}': n={len(sel):3d}  "
                  f"mean {np.mean(sel):.3f} s  "
                  f"sd {np.std(sel)*1000:.0f} ms  "
                  f"(published {target:.1f} s)")

    # --- 4) the 10 s marker comb
    t_m = (s0 / fs2 + np.arange(len(kinds)))[gate & (kinds == "M")]
    if len(t_m) >= 3:
        ang = (t_m % 10) / 10 * 2 * np.pi
        R = np.abs(np.mean(np.exp(1j * ang)))
        mu = np.angle(np.mean(np.exp(1j * ang))) / (2 * np.pi) * 10 % 10
        print(f"marker comb: {len(t_m)} confident markers, "
              f"phase {mu:.2f} s mod 10, coherence R = {R:.3f}")
    else:
        print("marker comb: too few confident markers — inconclusive")

    # --- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    tt = np.arange(L) / fs2
    for lab, color in (("0", "#0aa"), ("1", "#d80"), ("M", "#d33")):
        sel = np.where(gate & (kinds == lab))[0]
        if not len(sel):
            continue
        segs = np.array([pws[s0 + kk * L: s0 + kk * L + L] for kk in sel])
        a1.plot(tt, segs.mean(0) / hi, color=color, lw=1.6,
                label=f"'{lab}' seconds (n={len(sel)})")
    for wd in (0.2, 0.5, 0.8):
        a1.axvline(wd, color="#888", ls=":", lw=.8)
    a1.set(title="WWVB: average second, by pulse class\n"
                 "(dotted: published 0.2 / 0.5 / 0.8 s restore times)",
           xlabel="time into second (s)", ylabel="carrier power (rel.)")
    a1.legend(loc="lower right", fontsize=8)
    ok = gate & ~np.isnan(restore)
    a2.hist(restore[ok], bins=np.arange(0.10, 1.0, 0.025),
            color="#0aa", edgecolor="none")
    for wd in (0.2, 0.5, 0.8):
        a2.axvline(wd, color="#d33", ls=":", lw=1)
    a2.set(title="restore-time histogram, all confident seconds",
           xlabel="carrier restore time (s)", ylabel="seconds")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "wwvb.png", dpi=110)
    print(f"figure -> {out / 'wwvb.png'}")

    # fold the whole capture at 10 s: ten per-second notches, one of
    # them (the marker) visibly longer — the 10 s comb in one picture
    L10 = 10 * L
    m10 = (len(pws) - s0) // L10
    fold10 = pws[s0: s0 + m10 * L10].reshape(m10, L10).mean(0)
    fig2, ax = plt.subplots(figsize=(11.5, 3.4))
    ax.plot(np.arange(L10) / fs2, fold10 / hi, lw=1.2, color="#0aa")
    for kk in range(10):
        ax.axvline(kk, color="#ccc", lw=.5)
    ax.set(title=f"capture folded at 10 s ({m10} folds): a notch every "
                 "second, and one longer notch — the frame marker",
           xlabel="time mod 10 s", ylabel="carrier power (rel.)")
    fig2.tight_layout()
    fig2.savefig(out / "wwvb_fold10.png", dpi=110)
    print(f"figure -> {out / 'wwvb_fold10.png'}")


if __name__ == "__main__":
    main()
