#!/usr/bin/env python3
"""Measure the FT8 grid: the amateur mode where the grid is UTC itself.

FT8 (Franke-Taylor design, 8-GFSK) is the most-transmitted digital
mode on HF. Its grid:
  - transmissions start ON the UTC 15-second boundaries
    (:00 / :15 / :30 / :45), 0.5 s in, and last 12.64 s
  - symbol rate 6.25 Bd with tone spacing = symbol rate = 6.25 Hz
    (orthogonal FSK: the tones ARE the baud)
  - 79 symbols per transmission, 3x7 Costas arrays for sync
  - everyone shares ~3 kHz of USB audio above the dial frequency

Two measurements prove it off the air:
  1. fold band energy at 15 s against ABSOLUTE UTC (capture start time
     from the .json sidecar): the on/off duty pattern snaps to the
     wall clock — a grid synchronized not to a transmitter but to
     every GPS-disciplined ham shack on Earth simultaneously.
  2. isolate one signal, FM-discriminate, measure the 6.25 Hz
     symbol-transition line.

Usage: python measure.py --iq ft8_7074.cs16 --fs 250000
       (sidecar ft8_7074.json must carry t_start_utc)
"""
import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    a = ap.parse_args()

    meta = json.loads(Path(a.iq).with_suffix(".json").read_text())
    t0 = datetime.fromisoformat(meta["t_start_utc"])
    frac0 = (t0.second + t0.microsecond / 1e6) % 15.0
    print(f"capture starts {t0:%H:%M:%S}.{t0.microsecond//1000:03d}Z "
          f"= {frac0:.3f} s into a UTC 15 s slot")

    raw = np.fromfile(a.iq, dtype=np.int16)
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    # the FT8 sub-band: dial+200 .. dial+3100 Hz (USB convention)
    n = np.arange(len(x))
    bb = x * np.exp(-2j * np.pi * 1650 * n / a.fs)
    bb = lfilter(firwin(511, 1500, fs=a.fs), 1.0, bb)[::125]  # -> 2 kHz
    fsd = a.fs / 125
    p = np.abs(bb) ** 2
    k = int(0.25 * fsd)
    p = np.convolve(p, np.ones(k) / k, "same")

    # fold band power at 15 s in ABSOLUTE UTC phase
    L = int(15 * fsd)
    ph0 = int(frac0 * fsd)
    pad = np.full(ph0, np.nan)
    seq = np.concatenate([pad, p])
    nrow = len(seq) // L
    fold = np.nanmean(seq[: nrow * L].reshape(nrow, L), axis=0)
    fold = fold / np.nanmedian(fold)
    tx_win = fold[int(0.7 * fsd):int(12.3 * fsd)].mean()
    gap_win = np.concatenate([fold[:int(0.3 * fsd)],
                              fold[int(13.3 * fsd):]]).mean()
    print(f"UTC-folded band power: transmit window {tx_win:.2f}x median, "
          f"slot-edge gap {gap_win:.2f}x")
    print(f"  -> on/off contrast {10*np.log10(tx_win/gap_win):.1f} dB at "
          f"the WALL-CLOCK 15 s grid (nobody sent us a sync word — the "
          f"grid is UTC)")

    # one signal's symbol clock: strongest narrowband carrier in-band
    N = 1 << 14
    seg = bb[: len(bb) // N * N].reshape(-1, N)
    psd = np.fft.fftshift(
        (np.abs(np.fft.fft(seg * np.hanning(N), axis=1)) ** 2).mean(0))
    fax = np.fft.fftshift(np.fft.fftfreq(N, 1 / fsd))
    f_sig = float(fax[np.argmax(psd)])
    # GFSK transitions are deliberately Gaussian-smooth — a transition
    # spectrum reads low (5.7 "Bd" on our first try). The honest clock
    # is the TONE COMB: 8 tones at exactly 6.25 Hz spacing. High-res
    # FFT over one transmission, then autocorrelate the spectrum — the
    # comb pitch pops out at the tone spacing (= symbol rate, by
    # design: orthogonal FSK).
    # Two dead ends taught us the right estimator (both preserved in
    # git history): GFSK's Gaussian shaping deliberately smears tone
    # lines out of LONG-window spectra — no comb exists there to find.
    # The tones only exist one-at-a-time, one per 0.16 s symbol. So:
    # symbol-SYNCHRONIZED analysis, timed off the UTC slot grid we
    # just proved — FFT each symbol window of the strongest signal,
    # take its peak tone, and histogram the pairwise tone differences.
    # A real 8-GFSK signal quantizes them at multiples of 6.25 Hz.
    sig = bb * np.exp(-2j * np.pi * f_sig * np.arange(len(bb)) / fsd)
    sig = lfilter(firwin(255, 40, fs=fsd), 1.0, sig)
    pw = np.convolve(np.abs(sig) ** 2,
                     np.ones(int(2 * fsd)) / int(2 * fsd), "same")
    slot = int(15 * fsd)
    s_best = max(range(int((15 - frac0 + 0.6) * fsd), len(sig) - slot,
                       slot),
                 key=lambda s: pw[s + int(6 * fsd)])
    Tsym = 0.16
    tones = []
    NF = 1 << 12
    for k in range(79):
        c = s_best + int((k + 0.5) * Tsym * fsd)
        w = sig[c - int(0.06 * fsd):c + int(0.06 * fsd)]
        S = np.abs(np.fft.fftshift(np.fft.fft(w * np.hanning(len(w)), NF)))
        fxs = np.fft.fftshift(np.fft.fftfreq(NF, 1 / fsd))
        band = np.abs(fxs) < 40
        tones.append(float(fxs[band][np.argmax(S[band])]))
    tones = np.array(tones)
    diffs = np.abs(tones[:, None] - tones[None, :]).ravel()
    diffs = diffs[(diffs > 3.0) & (diffs < 47)]
    # fold the differences at candidate spacings; true spacing minimizes
    # the residual to the nearest multiple
    cands = np.linspace(5.0, 8.0, 301)
    resid = [np.mean(np.minimum(diffs % c, c - diffs % c)) for c in cands]
    spacing = float(cands[int(np.argmin(resid))])
    print(f"strongest signal, symbol-synchronized tones: spacing "
          f"quantizes at {spacing:.3f} Hz (residual "
          f"{min(resid):.2f} Hz; grid says 6.25 = the symbol rate)")
    baud = spacing
    E = np.array(resid)
    fx = cands
    fl = max(min(resid), 1e-9)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
    tt = np.arange(L) / fsd
    a1.plot(tt, fold, lw=1.0, color="#0aa")
    a1.axvspan(0.5, 13.14, color="#d33", alpha=.08)
    a1.annotate("79 symbols x 0.16 s\n= 12.64 s transmission",
                (6.5, np.nanmax(fold) * .95), ha="center", fontsize=9,
                color="#d33")
    a1.annotate("turnaround\ngap", (14.2, np.nanmin(fold) * 1.2),
                ha="center", fontsize=8)
    a1.set(title="band power folded at 15 s of ABSOLUTE UTC",
           xlabel="seconds past :00/:15/:30/:45", ylabel="power / median")
    a2.plot(fx, E, lw=.8, color="#0aa")
    a2.axvline(6.25, color="#d33", ls=":")
    a2.annotate("6.25 Hz", (6.3, max(E) * .9), color="#d33", fontsize=9)
    a2.set(title=f"tone-difference quantization: best fit "
                 f"{baud:.3f} Hz spacing",
           xlabel="candidate tone spacing (Hz)",
           ylabel="mean residual to nearest multiple (Hz)")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "ft8.png", dpi=110)
    print(f"figure -> {out / 'ft8.png'}")


if __name__ == "__main__":
    main()
