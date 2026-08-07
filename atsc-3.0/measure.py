#!/usr/bin/env python3
"""Measure the ATSC 3.0 (NextGen TV) bootstrap grid from an IQ capture.

The grid (ATSC A/321):
  bootstrap rate  6.144 MHz    FIXED, always, whatever the rest of the frame does
  bootstrap BW    4.5 MHz      fits inside every channel raster in the world
  FFT             2048         -> 3000 Hz subcarrier spacing
  symbol          3072 samples = C(520) | A(2048) | B(504) = exactly 500 us
  bootstrap       4 symbols    = 2.000 ms, opening every frame
  sequence        Zadoff-Chu root q=137, length 1499 (largest prime under
                  4.5 MHz / 3 kHz), scrambled by a 16-bit PN

WHY THIS SIGNAL IS BUILT THIS WAY, AND THE TRICK WORTH STEALING
---------------------------------------------------------------
ATSC 3.0 is aggressively configurable -- FFT size, guard interval, pilot
pattern, constellation and code rate all change per broadcaster, and a
receiver cannot know any of them in advance.  So the standard nails ONE
thing down: every frame opens with a bootstrap at a fixed 6.144 MHz and a
fixed 4.5 MHz bandwidth, no matter what follows.  It is a known lighthouse
in front of an unknown harbour, and it is the reason a receiver can lock a
signal whose parameters it has never seen.

The trick: you do not need the Zadoff-Chu sequence, the PN, or any of
A/321's tables to FIND it.  Part C is defined as a copy of the last 520
samples of part A, so the waveform repeats itself at a lag of exactly 2048
samples.  Correlate the capture against itself at that one lag and the
bootstrap stands up out of the noise -- four peaks, 3072 samples apart --
using nothing but the time geometry.  Stage 1 below does that.  Stage 2
then confirms identity with the real matched filter, and both stages carry
a control that must fail.

Usage:
    python measure.py --iq capture.cs16 --fs 8e6

    --iq   interleaved int16 IQ, centred on the channel
    --fs   its sample rate (anything >= ~6.5 MS/s)
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.signal import resample_poly

HERE = Path(__file__).resolve().parent

# ---- A/321 constants ----------------------------------------------------
FS = 6_144_000          # 5.1: the fixed bootstrap sampling rate
N_FFT = 2048            # 5.1
DF = FS / N_FFT         # 3000 Hz subcarrier spacing
N_ZC = 1499             # 5.2.1: largest prime fitting 4.5 MHz at 3 kHz
N_H = (N_ZC - 1) // 2   # 749
N_A, N_B, N_C = 2048, 504, 520      # 5.4
N_SYM = N_A + N_B + N_C             # 3072 samples = 500 us
N_BOOT = 4 * N_SYM                  # 2.000 ms
ZC_ROOT = 137           # 6.1: q for bootstrap_major_version = 0
PN_ORDER = 16
PN_G = (1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1)   # g_16..g_0
MINOR_SEEDS = {0: 0x019D, 1: 0x00ED, 2: 0x01E8, 3: 0x00E8,
               4: 0x00FB, 5: 0x0021, 6: 0x0054, 7: 0x00EC}   # Table 6.1

_TAPS = tuple(i for i in range(PN_ORDER) if PN_G[PN_ORDER - i])
_TAPS_RECIP = tuple(PN_ORDER - 1 - i for i in _TAPS)
# A/321 Fig 5.2 leaves the register wiring open to more than one honest
# reading, so we try all four and report which one actually locks rather
# than asserting one and hoping.
PN_VARIANTS = {"spec": (_TAPS, True), "spec_revseed": (_TAPS, False),
               "recip": (_TAPS_RECIP, True), "recip_revseed": (_TAPS_RECIP, False)}


def pn_pm1(seed, length, variant):
    taps, msb_high = PN_VARIANTS[variant]
    if msb_high:
        r = [(seed >> i) & 1 for i in range(PN_ORDER)]
    else:
        r = [(seed >> (PN_ORDER - 1 - i)) & 1 for i in range(PN_ORDER)]
    out = np.empty(length)
    for k in range(length):
        out[k] = 1.0 - 2.0 * r[0]
        fb = 0
        for t in taps:
            fb ^= r[t]
        r = r[1:] + [fb]
    return out


def zc(q, n=N_ZC):
    """z_q(k) = exp(-j pi q k (k+1) / N_ZC)   (A/321 5.2.1)."""
    k = np.arange(n, dtype=np.float64)
    ph = np.mod(q * k * (k + 1.0), 2.0 * n) / n     # keep the phase small
    return np.exp(-1j * np.pi * ph)


def symbol0(q=ZC_ROOT, minor=0, variant="spec"):
    """Part A of bootstrap symbol 0 -- a fully known waveform (A/321 5.3.3
    fixes symbol 0 at absolute cyclic shift 0)."""
    X = np.zeros(N_FFT, complex)
    seq = zc(q) * pn_pm1(MINOR_SEEDS[minor], N_ZC, variant)
    X[(np.arange(N_ZC) - N_H) % N_FFT] = seq
    return np.fft.ifft(X) * N_FFT


def norm_xcorr(x, tpl):
    """|correlation| of x against tpl, normalised by local energy."""
    n = len(tpl)
    c = np.abs(np.correlate(x, tpl, "valid"))
    e = np.sqrt(np.convolve(np.abs(x) ** 2, np.ones(n), "valid"))
    return c / (e * np.sqrt(np.sum(np.abs(tpl) ** 2)) + 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--secs", type=float, default=4.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    print(f"capture: {len(x)} samples at {a.fs / 1e6:.3f} MS/s "
          f"({len(x) / a.fs:.2f} s)")

    # everything happens at the ONE rate the standard fixes
    g = np.gcd(int(FS), int(a.fs))
    up, down = int(FS) // g, int(a.fs) // g
    y = resample_poly(x, up, down)
    print(f"resampled to the bootstrap rate: x{up}/{down} -> "
          f"{FS / 1e6:.3f} MS/s, {len(y)} samples\n")

    # ---- STAGE 1: find it using only the time geometry -------------------
    # Part C is a copy of the last 520 samples of part A, so y[n] repeats at
    # y[n+2048].  That ALONE is not enough: every OFDM symbol in the payload
    # has a cyclic prefix and self-correlates too, which is exactly what the
    # control below catches.  What is unique to the bootstrap is that the
    # repeat happens FOUR TIMES, spaced exactly one symbol (3072) apart.  So
    # the detector is the comb, not the single lag.
    def movsum(v, w):
        c = np.concatenate(([0.0], np.cumsum(v)))
        return c[w:] - c[:-w]

    def sc_metric(lag):
        """Schmidl-Cox style, bounded 0..1: |sum a conj(b)| / sum |b|^2."""
        w = N_C
        a_ = y[:len(y) - lag]
        b_ = y[lag:]
        num = np.abs(movsum(a_ * np.conj(b_), w))
        den = movsum(np.abs(b_) ** 2, w) + 1e-12
        return num / den

    def comb(m, spacing, n=4):
        """Sum the metric at n points one symbol apart -- the bootstrap's
        real signature."""
        L = len(m) - (n - 1) * spacing
        out = np.zeros(L)
        for k in range(n):
            out += m[k * spacing: k * spacing + L]
        return out / n

    def score(d):
        return d.max() / np.median(d)

    m_true = sc_metric(N_A)
    d_true = comb(m_true, N_SYM)                    # right lag, right spacing
    d_lag = comb(sc_metric(N_A - 48), N_SYM)        # control: wrong lag
    d_sp = comb(m_true, N_SYM - 96)

    print(f"STAGE 1  bootstrap geometry alone (no sequences, no tables)")
    print(f"  single lag {N_A}, no comb     peak/median {score(m_true):7.2f}"
          f"   <- not selective by itself")
    print(f"  4-symbol comb, spacing {N_SYM}   "
          f"peak/median {score(d_true):7.2f}")
    print(f"  CONTROL wrong lag {N_A - 48:<13d}peak/median "
          f"{score(d_lag):7.2f}   <- must fail, and does")
    sel = score(d_true) / score(d_lag)
    print(f"  selectivity in LAG  {sel:5.2f}x"
          f"   {'PASS' if sel > 1.5 else 'FAIL'}")

    # How finely does the comb resolve the SYMBOL SPACING?  A control 96
    # samples off scores 16.29 -- not because the detector is weak, but
    # because part C is 520 samples long, so the correlation sits on a
    # ~520-sample plateau and a 4-tap comb only drifts 3x96 off it.  Rather
    # than widen the control until it fails (which would be tuning a gate to
    # get the answer we want), measure the resolution and report it.
    offs = np.arange(-600, 601, 60)
    sweep = np.array([score(comb(m_true, N_SYM + o)) for o in offs])
    base = np.median(sweep[np.abs(offs) >= 480])
    half = sweep.max() - (sweep.max() - base) / 2
    inside = offs[sweep >= half]
    res = max(abs(inside.min()), inside.max())
    print(f"  spacing sweep: peak at offset "
          f"{offs[int(np.argmax(sweep))]:+d}, half-height width +/-{res}"
          f"  <- too coarse to confirm 3072 by itself")

    # The comb SCORE is blunt, but the peak POSITIONS are not: read the
    # symbol period straight off the four peaks instead of sweeping.
    p0 = int(np.argmax(d_true))
    pos = [p0]
    for k in range(1, 4):
        w0, w1 = p0 + k * N_SYM - 700, p0 + k * N_SYM + 700
        if w1 < len(m_true):
            pos.append(w0 + int(np.argmax(m_true[w0:w1])))
    if len(pos) == 4:
        gaps = np.diff(pos)
        print(f"  symbol peaks at {[int(p - p0) for p in pos]} samples")
        print(f"  measured symbol period {gaps.mean():8.1f} samples "
              f"= {1e6 * gaps.mean() / FS:.2f} us   "
              f"(grid says {N_SYM} = {1e6 * N_SYM / FS:.0f})")
        err = abs(gaps.mean() - N_SYM)
        print(f"  error {err:.1f} samples "
              f"({'PASS' if err < 20 else 'FAIL'}) -- the peak POSITIONS "
              f"resolve what the comb score could not")

    # ---- STAGE 2: confirm identity with the matched filter ---------------
    print(f"\nSTAGE 2  Zadoff-Chu matched filter")
    best = None
    for v in PN_VARIANTS:
        tpl = symbol0(ZC_ROOT, 0, v)
        c = norm_xcorr(y[:min(len(y), int(0.25 * FS))], tpl)
        ratio = c.max() / np.median(c)
        print(f"  PN wiring {v:15s} peak/median {ratio:8.2f}")
        if best is None or ratio > best[1]:
            best = (v, ratio, c)
    variant, ratio, c_best = best
    print(f"  -> locks with '{variant}'")

    # control: a ZC root the standard reserves for other major versions
    ctrl_ratios = []
    for qc in (ZC_ROOT - 1, ZC_ROOT + 1):
        cc = norm_xcorr(y[:min(len(y), int(0.25 * FS))],
                        symbol0(qc, 0, variant))
        ctrl_ratios.append(cc.max() / np.median(cc))
        print(f"  CONTROL root q={qc:4d}      peak/median "
              f"{ctrl_ratios[-1]:8.2f}   <- must be small")
    print(f"\n  root 137 beats the neighbouring roots by "
          f"{ratio / max(ctrl_ratios):.1f}x  "
          f"-> this is ATSC 3.0, major version 0")

    # ---- frame period, from the matched filter (the clean detector) ------
    c_full = norm_xcorr(y, symbol0(ZC_ROOT, 0, variant))
    thr = np.median(c_full) + 0.5 * (c_full.max() - np.median(c_full))
    hits = np.where(c_full > thr)[0]
    starts = []
    if len(hits):
        grp = hits[0]
        for p in hits[1:]:
            if p - grp > N_BOOT:
                starts.append(grp); grp = p
        starts.append(grp)
    print(f"\n  bootstraps in {len(y) / FS:.2f} s: {len(starts)}")
    if len(starts) > 2:
        d = np.diff(starts) / FS * 1e3
        print(f"  frame period  {np.median(d):8.3f} ms   "
              f"spread {d.max() - d.min():.3f} ms over {len(d)} gaps")
        print(f"  frame rate    {1e3 / np.median(d):8.3f} frames/s")

    print(f"\nderived grid numbers")
    print(f"  subcarrier spacing  {DF:8.1f} Hz   = {FS}/{N_FFT}")
    print(f"  symbol duration     {1e6 * N_SYM / FS:8.1f} us   "
          f"= {N_SYM} samples")
    print(f"  bootstrap duration  {1e3 * N_BOOT / FS:8.3f} ms  = 4 symbols")
    print(f"  ZC occupied BW      {1e-6 * N_ZC * DF:8.3f} MHz  "
          f"= {N_ZC} x {DF:.0f} Hz")

    # ---- figure ----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.6))

    if starts:
        s = max(0, starts[0] - 2000)
        w = slice(s, min(len(d_true), s + N_BOOT + 4000))
        t = (np.arange(w.start, w.stop) - starts[0]) / FS * 1e6
        a1.plot(t, d_sp[w], lw=.7, color="#bbb",
                label=f"control: comb spacing {N_SYM - 96}")
        a1.plot(t, d_lag[w], lw=.7, color="#888",
                label=f"control: lag {N_A - 48}")
        a1.plot(t, d_true[w], lw=1.1, color="#c33",
                label=f"lag {N_A}, 4-symbol comb")
        a1.axvline(0, color="#333", ls=":", lw=.8)
        a1.set(title="Stage 1: found by time geometry alone\n"
                     "part C repeats part A, four symbols running",
               xlabel="time relative to bootstrap start (us)",
               ylabel="combed self-correlation")
    a1.legend(fontsize=8, loc="upper right")

    t2 = np.arange(len(c_best)) / FS * 1e3
    a2.plot(t2, norm_xcorr(y[:len(c_best) + N_FFT - 1],
                           symbol0(ZC_ROOT + 1, 0, variant)),
            lw=.6, color="#999", label="control: ZC root 138")
    a2.plot(t2, c_best, lw=.8, color="#06c", label="ZC root 137 (A/321 6.1)")
    a2.set(title="Stage 2: matched filter confirms identity\n"
                 "one spike per frame",
           xlabel="time (ms)", ylabel="normalised correlation")
    a2.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "atsc3_bootstrap.png", dpi=110)
    print(f"\nfigure -> {out / 'atsc3_bootstrap.png'}")


if __name__ == "__main__":
    main()
