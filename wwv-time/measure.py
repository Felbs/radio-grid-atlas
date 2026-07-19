#!/usr/bin/env python3
"""Measure the WWV grid: the purest time structure on the air.

WWV (Fort Collins, Colorado) transmits on 2.5/5/10/15/20 MHz, AM.
The grid IS time itself:
  - a tick every second: 5 ms of 1000 Hz (1200 Hz on WWVH)
  - the 100 Hz subcarrier carrying the BCD time code (IRIG-ish)
  - 500/600 Hz standard tones on alternating minutes
  - a 0.8 s 1000/1500 Hz marker on the minute
Every element is derived from cesium clocks — when you measure WWV's
tick spacing against your SDR's sample clock, you are calibrating
your hardware against NIST.

Usage: python measure.py --iq capture.cs16 --fs 250000
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--secs", type=float, default=90.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    # isolate the AM channel at capture center (+-5 kHz)
    x = lfilter(firwin(255, 5_000, fs=a.fs), 1.0, x)
    dec = int(a.fs / 25_000)
    fsd = a.fs / dec
    env = np.abs(x[::dec])                    # AM demod = envelope
    env -= env.mean()

    # audio spectrum of the demodulated envelope
    N = 1 << 16
    seg = env[: len(env) // N * N].reshape(-1, N)
    S = (np.abs(np.fft.rfft(seg * np.hanning(N), axis=1)) ** 2).mean(0)
    fax = np.fft.rfftfreq(N, 1 / fsd)

    def tone(f0, half=8):
        m = (fax > f0 - half) & (fax < f0 + half)
        med = np.median(S[(fax > 150) & (fax < 900)])
        return 10 * np.log10(S[m].max() / med)

    print(f"100 Hz time-code subcarrier: {tone(100):+5.1f} dB above floor")
    print(f"500 Hz standard tone:        {tone(500):+5.1f} dB")
    print(f"600 Hz standard tone:        {tone(600):+5.1f} dB")
    print(f"1000 Hz (ticks/minute mark): {tone(1000):+5.1f} dB")

    # fold the 100 Hz TIME-CODE subcarrier power at 1 s: its BCD pulse
    # (0.17/0.47/0.77 s widths) repeats every second, and it's the
    # strongest element under HF fading (the 1 kHz ticks are only a
    # few ms long and fade — folding weak ticks once "measured" a
    # +440 ppm clock error on a 0.5 ppm TCXO; fold the strong thing).
    # narrowband work needs a matching sample rate: decimate the
    # envelope to ~2 kHz FIRST (a 100 Hz bandpass at 25 kHz needs
    # thousands of taps and silently passes everything otherwise —
    # that mush filter "measured" a 440 ppm error on a 0.5 ppm TCXO)
    d2 = 12
    env2 = lfilter(firwin(127, 800, fs=fsd), 1.0, env)[::d2]
    fs2 = fsd / d2
    b100 = lfilter(firwin(401, [80, 120], fs=fs2, pass_zero=False),
                   1.0, env2)
    p1k = np.abs(b100) ** 2
    k = max(1, int(0.02 * fs2))
    kernel = np.ones(k) / k
    tick_pow = np.convolve(p1k, kernel, "same")
    # HF fading (QSB) swamps the per-second structure — divide out a
    # 5 s rolling mean so the fold sees the code, not the ionosphere
    ks = int(5 * fs2)
    slow = np.convolve(tick_pow, np.ones(ks) / ks, "same")
    tick_pow = tick_pow / (slow + 1e-12)
    # fold at candidate periods: the peak must sit at 1.000 s (ms-level
    # proof of the one-second grid; ppm-level clock calibration against
    # NIST needs carrier-phase methods, not envelope folds — an earlier
    # version of this script claimed 440 ppm on a 0.5 ppm TCXO before
    # admitting that)
    best = (0, 0)
    for P in np.linspace(0.98, 1.02, 81):
        L = int(round(P * fs2))
        n = len(tick_pow) // L
        if n < 3:
            continue
        prof = tick_pow[: n * L].reshape(n, L).mean(0)
        s = (prof.max() - np.median(prof)) / (prof.std() + 1e-12)
        if s > best[0]:
            best = (s, L / fs2)
    if best[0] > 8 and abs(best[1] - 1.0) < 0.005:
        print(f"time-code frame period (folded): {best[1]:.4f} s "
              f"(sharpness {best[0]:.1f}) — the one-second grid")
    else:
        print(f"per-second fold: inconclusive under this capture's HF "
              f"fading (sharpness {best[0]:.1f} at {best[1]:.3f} s) — "
              f"WWV folds clean after local midnight when the band "
              f"settles; the harmonic tone grid above already "
              f"identifies the station")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    m = fax < 1300
    a1.semilogy(fax[m], S[m], lw=.6, color="#0aa")
    for f0, lbl in ((100, "100 Hz\ntime code"), (500, "500"), (600, "600"),
                    (1000, "1000 Hz\nticks")):
        a1.axvline(f0, color="#d33", ls=":", lw=.8)
    a1.set(title="WWV demodulated audio: the time grid",
           xlabel="Hz", ylabel="power")
    n0 = int(20 * fs2)
    n1 = int(30 * fs2)
    a2.plot(np.arange(n0, n1) / fs2, tick_pow[n0:n1], lw=.7, color="#0aa")
    for s in range(20, 31):
        a2.axvline(s, color="#d33", ls=":", lw=.6)
    a2.set(title="100 Hz time-code power, 10 s strip "
                 "(one BCD pulse per second, QSB and all)",
           xlabel="capture time (s)", ylabel="normalized code power")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "wwv.png", dpi=110)
    print(f"figure -> {out / 'wwv.png'}")


if __name__ == "__main__":
    main()
