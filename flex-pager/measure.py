#!/usr/bin/env python3
"""Measure the FLEX pager grid off a 929 MHz band capture.

FLEX (Motorola, 1993 — still carrying hospital pages today) divides
time into 1.875 s frames, 128 frames per 4-minute cycle. Each frame:
sync (sent robust) + frame info + 11 data blocks, at 1600 / 3200 /
6400 bps. High-speed frames use 4-FSK with deviations +-4800 Hz
(outer symbols) and +-1600 Hz (inner symbols) at 3200 symbols/s.

This script derives, from raw IQ of a paging channel:
  1. the frame comb — burst start times fold at exactly 1.875 s
  2. the 4-FSK deviation grid (+-1600 / +-4800 Hz)
  3. the 3200 Hz symbol clock (level-change spectrum)
and runs one negative control: the POCSAG sync word 0x7CD215D8 at
2400 Bd, which this band does NOT contain.

Usage: python measure.py --iq capture.cs16 --fs 1000000
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

HERE = Path(__file__).resolve().parent
POCSAG_SYNC = 0x7CD215D8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--f0", type=float, default=None,
                    help="channel offset from dial, Hz (default: auto)")
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16)
    x = raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)
    del raw
    n_ax = np.arange(len(x), dtype=np.float64)

    # --- pick the busiest channel by mean PSD (25 kHz smoothing)
    if a.f0 is None:
        N = 4096
        m = len(x) // N
        S = (np.abs(np.fft.fft(x[: m * N].reshape(m, N)[::7],
                               axis=1)) ** 2).mean(0)
        S = np.fft.fftshift(S)
        fax = np.fft.fftshift(np.fft.fftfreq(N, 1 / a.fs))
        kw = max(1, int(20e3 / (a.fs / N)))
        Ssm = np.convolve(S, np.ones(kw) / kw, "same")
        f0 = fax[np.argmax(Ssm)]
    else:
        f0 = a.f0
    print(f"channel offset from dial: {f0:+.0f} Hz")

    # --- channelize: +-12 kHz, decimate to 50 kHz
    bb = x * np.exp(-2j * np.pi * f0 * n_ax / a.fs).astype(np.complex64)
    del x
    dec = int(a.fs / 50e3)
    bb = lfilter(firwin(255, 12e3, fs=a.fs).astype(np.float32),
                 1.0, bb)[::dec]
    fsd = a.fs / dec

    # --- burst detection on channel power
    p = np.convolve(np.abs(bb) ** 2, np.ones(500) / 500, "same")
    on = p > np.median(p) * 8
    d_on = np.diff(on.astype(int))
    st = np.where(d_on == 1)[0] + 1
    en = np.where(d_on == -1)[0] + 1
    if on[0]:
        st = np.insert(st, 0, 0)
    if on[-1]:
        en = np.append(en, len(on))
    bursts = []
    cs, ce = st[0], en[0]
    for s, e in zip(st[1:], en[1:]):
        if s - ce < 0.05 * fsd:
            ce = e
        else:
            bursts.append((cs, ce))
            cs, ce = s, e
    bursts.append((cs, ce))
    bursts = [(s, e) for s, e in bursts if e - s > 0.15 * fsd]
    starts = np.array([s / fsd for s, e in bursts])
    print(f"{len(bursts)} bursts; starts {starts.round(2).tolist()}")

    # --- 1) the frame comb: scan fold period, require phase coherence
    # (skip any burst clipped at the capture edge — its start is fake)
    ok = starts > 0.1
    best = (0.0, 0.0)
    for P in np.linspace(1.5, 2.5, 10001):
        ang = (starts[ok] % P) / P * 2 * np.pi
        R = np.abs(np.mean(np.exp(1j * ang)))
        if R > best[0]:
            best = (R, P)
    R, P = best
    d = np.diff(starts)
    print(f"frame comb: best period {P:.4f} s (R = {R:.3f}) — "
          f"published FLEX frame = 1.8750 s")
    print(f"  burst spacings / 1.875: "
          f"{(d / 1.875).round(3).tolist()}")

    # --- 2) 4-FSK deviations from the longest burst
    disc = np.angle(bb[1:] * np.conj(bb[:-1])) * fsd / (2 * np.pi)
    s, e = max(bursts, key=lambda b: b[1] - b[0])
    dd = disc[s + 100:e - 100]
    ctr = (np.percentile(dd, 99.5) + np.percentile(dd, 0.5)) / 2
    dc = dd - ctr
    h, edges = np.histogram(dc[np.abs(dc) < 8000], bins=160)
    hs = np.convolve(h, np.ones(5) / 5, "same")
    pk = [i for i in range(2, 158)
          if hs[i] == hs[i - 2:i + 3].max() and hs[i] > 0.12 * hs.max()]
    lobes = np.array([(edges[i] + edges[i + 1]) / 2 for i in pk])
    print(f"FSK lobes (Hz, carrier removed): {lobes.round(0).tolist()}")
    if len(lobes) == 4:
        outer = (lobes[3] - lobes[0]) / 2
        inner = (lobes[2] - lobes[1]) / 2
        print(f"  4-FSK deviation: outer +-{outer:.0f} Hz "
              f"(published 4800), inner +-{inner:.0f} Hz "
              f"(published 1600)")

    # --- 3) symbol clock: spectrum of level-change events
    tr = np.abs(np.diff(dc))
    tr[tr < 1500] = 0.0
    Nf = 1 << 20
    T = np.abs(np.fft.rfft(tr - tr.mean(), Nf))
    fx = np.fft.rfftfreq(Nf, 1 / fsd)
    mm = (fx > 500) & (fx < 8000)
    fl = np.median(T[mm])
    i0 = np.argmax(T[mm])
    f_sym = fx[mm][i0]
    print(f"symbol clock: {f_sym:.2f} Hz "
          f"({20*np.log10(T[mm][i0]/fl):+.1f} dB) — "
          f"published 3200 sym/s for 6400 bps 4-FSK")

    # --- negative control: POCSAG sync at 2400 Bd, all bursts
    pat = (2 * np.array([int(c) for c in f"{POCSAG_SYNC:032b}"],
                        np.int8) - 1).astype(np.float32)
    hits = 0
    for s, e in bursts:
        db = disc[s:e]
        db = db - (np.percentile(db, 85) + np.percentile(db, 15)) / 2
        sps = fsd / 2400.0
        best_ph, best_q = 0.0, -1.0
        for ph in np.linspace(0, sps, 24, endpoint=False):
            idx = (np.arange(int((len(db) - ph - 1) / sps)) * sps
                   + ph).astype(int)
            q = np.abs(db[idx]).mean()
            if q > best_q:
                best_q, best_ph = q, ph
        idx = (np.arange(int((len(db) - best_ph - 1) / sps)) * sps
               + best_ph).astype(int)
        bits = np.sign(db[idx]).astype(np.float32)
        c = np.correlate(bits, pat)
        hits += int((np.abs(c) >= 24).sum())   # >= 28/32 bits, either pol
    print(f"POCSAG sync control: {hits} hits in {len(bursts)} bursts "
          f"(this band is FLEX, not POCSAG)")

    # --- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    cx = (edges[:-1] + edges[1:]) / 2
    a1.fill_between(cx, h, color="#0aa", step="mid")
    for v in (-4800, -1600, 1600, 4800):
        a1.axvline(v, color="#d33", ls=":", lw=1)
    a1.set(title="instantaneous frequency, longest burst\n"
                 "(dotted: published FLEX 4-FSK +-1600 / +-4800 Hz)",
           xlabel="Hz from carrier", ylabel="samples")
    mm2 = (fx > 2800) & (fx < 3600)
    a2.semilogy(fx[mm2], T[mm2] + 1e-9, lw=.7, color="#0aa")
    a2.axvline(3200, color="#d33", ls=":", lw=1)
    a2.set(title=f"level-change spectrum\nsymbol clock "
                 f"{f_sym:.2f} Hz (dotted: published 3200)",
           xlabel="Hz", ylabel="amplitude")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "flex.png", dpi=110)
    print(f"figure -> {out / 'flex.png'}")

    fig2, ax = plt.subplots(figsize=(11.5, 3.2))
    dp = max(1, int(fsd / 1000))
    tt = np.arange(0, len(p), dp) / fsd
    ax.plot(tt, 10 * np.log10(p[::dp] / np.median(p)), lw=.5,
            color="#0aa")
    ph0 = np.angle(np.mean(np.exp(1j * (starts[ok] % P) / P
                                  * 2 * np.pi))) / (2 * np.pi) * P % P
    for k in range(int(tt[-1] / P) + 1):
        ax.axvline(ph0 + k * P, color="#d33", lw=.4, alpha=.5)
    ax.set(title=f"channel power vs time — every burst starts on the "
                 f"{P:.4f} s frame comb (red)",
           xlabel="capture time (s)", ylabel="dB over floor")
    fig2.tight_layout()
    fig2.savefig(out / "flex_frames.png", dpi=110)
    print(f"figure -> {out / 'flex_frames.png'}")


if __name__ == "__main__":
    main()
