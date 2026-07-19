#!/usr/bin/env python3
"""Measure the RDS grid: 57 kHz BPSK at exactly 1187.5 bps.

The magic numbers are not arbitrary: 57 kHz = 3 x the 19 kHz stereo
pilot (so it's free to regenerate), and 1187.5 bps = 57,000 / 48 (so
the bit clock is phase-locked to the subcarrier itself — 48 subcarrier
cycles per bit). This script recovers the subcarrier, demodulates it,
and measures the symbol rate from the signal's own zero crossings.

Usage: python measure.py --iq capture.cs16 --fs 2976750
Output: figures/rds.png (constellation + eye) + printed rates.
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
    ap.add_argument("--secs", type=float, default=10.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    dec = max(1, int(a.fs / 300_000))
    fsd = a.fs / dec
    x = lfilter(firwin(91, 120_000, fs=a.fs), 1.0, x)[::dec]
    c = np.angle(x[1:] * np.conj(x[:-1]))

    # lock the 19 kHz pilot, cube its phasor -> exact 57 kHz reference
    n = np.arange(len(c))
    p = c * np.exp(-2j * np.pi * 19_000 * n / fsd)
    a1 = float(np.exp(-2 * np.pi * 400.0 / fsd))
    for _ in range(3):
        p = lfilter([1 - a1], [1, -a1], p)
    u = p / (np.abs(p) + 1e-9)
    pilot_hz = 19_000 + np.diff(np.unwrap(np.angle(u))).mean() * fsd / (2 * np.pi)
    print(f"pilot measured:   {pilot_hz:9.1f} Hz")
    print(f"57k reference:    {3 * pilot_hz:9.1f} Hz = 3 x pilot (by construction)")

    # demodulate RDS: mix by conj(u^3 * e^{j3w19}) and lowpass +-2.4k
    z = c * np.exp(-2j * np.pi * 57_000 * n / fsd) * np.conj(u) ** 3
    b = firwin(401, 2_400, fs=fsd)
    z = lfilter(b, 1.0, z)[len(b):]
    # carrier residue -> rotate so BPSK lands on the real axis
    ang = 0.5 * np.angle((z ** 2).mean())
    z = z * np.exp(-1j * ang)

    # symbol rate from the data itself: spectrum of |d/dt Re z|
    d = np.abs(np.diff(z.real))
    D = np.abs(np.fft.rfft(d - d.mean()))
    fax = np.fft.rfftfreq(len(d), 1 / fsd)
    m = (fax > 800) & (fax < 1600)
    baud = fax[m][np.argmax(D[m])]
    print(f"symbol rate:      {baud:9.2f} Bd   (grid says 1187.5 = 57000/48; "
          f"measured/57k = 1/{57_000 / baud:.2f})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1_, a2_) = plt.subplots(1, 2, figsize=(11, 4.6))
    sps = fsd / 1187.5
    k = (np.arange(0, len(z) - sps, sps / 2)).astype(int)[200:1800]
    a1_.plot(z.real[k], z.imag[k], ".", ms=2, alpha=.4, color="#0aa")
    a1_.set_title("RDS constellation (BPSK on the 57 kHz subcarrier)")
    a1_.set_xlabel("I")
    a1_.set_ylabel("Q")
    a1_.axhline(0, color="#888", lw=.5)
    a1_.axvline(0, color="#888", lw=.5)
    seg = int(2 * sps)
    tr = z.real[1000:1000 + 40 * seg]
    for i in range(0, len(tr) - seg, seg):
        a2_.plot(np.arange(seg) / sps, tr[i:i + seg], color="#0aa",
                 alpha=.25, lw=.8)
    a2_.set_title(f"eye diagram — measured {baud:.1f} Bd "
                  f"(57000/48 = 1187.5)")
    a2_.set_xlabel("bit periods")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "rds.png", dpi=110)
    print(f"figure -> {out / 'rds.png'}")


if __name__ == "__main__":
    main()
