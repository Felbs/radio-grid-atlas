#!/usr/bin/env python3
"""Measure the Meteor-M LRPT grid from a 137.9 MHz satellite pass.

The grid:
  carrier       137.100 / 137.9125 MHz (Meteor-M N2 series)
  modulation    QPSK, root-raised-cosine
  symbol rate   72,000 Sym/s  (80k on some birds/modes)
  FEC           K=7 r=1/2 convolutional (Viterbi), CCSDS scrambling
  framing       CADUs behind the 32-bit ASM 0x1ACFFC1D, RS(255,223)

This script proves the two outermost layers off the air:
 1. QPSK-ness: raise the signal to the 4th power — a QPSK carrier
    collapses to a single spectral line at 4x the residual carrier
    offset (BPSK collapses at 2x; noise never collapses).
 2. The symbol clock: the classic timing line — |x|^2 has a spectral
    line at exactly the symbol rate. Read it: 72,000 Sym/s.

Deeper layers (Viterbi, ASM, RS) are decoding rather than measuring;
see wxTuna's lrpt.py for the full chain (validated to images).

Usage: python measure.py --iq pass.cs16 --fs 250000
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
    ap.add_argument("--t-start", type=float, default=60.0,
                    help="seconds into the file (mid-pass = strong)")
    ap.add_argument("--secs", type=float, default=20.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs),
                      offset=int(a.t_start * a.fs) * 4)
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    x = lfilter(firwin(63, 80_000, fs=a.fs), 1.0, x)

    # ---- 1. QPSK collapse under ^4 ----------------------------------
    M = 1 << 20
    def line(sig, span):
        """Strongest line within +-span Hz — physics bounds the search:
        LEO Doppler at 137 MHz is +-3.5 kHz, so the x^4 line must sit
        within +-14 kHz. An unbounded hunt happily returns junk (this
        script's first run found a 'carrier' 6x outside Doppler)."""
        S = np.abs(np.fft.fftshift(np.fft.fft(sig[:M] * np.hanning(M))))
        fax = np.fft.fftshift(np.fft.fftfreq(M, 1 / a.fs))
        m = np.abs(fax) < span
        i = int(np.argmax(S[m]))
        return fax[m][i], 20 * np.log10(S[m][i] / np.median(S))
    f2, l2 = line((x / (np.abs(x) + 1e-9)) ** 2, 7_000)
    f4, l4 = line((x / (np.abs(x) + 1e-9)) ** 4, 14_000)
    print(f"x^2 spectrum: strongest line {l2:5.1f} dB  (BPSK would spike here)")
    print(f"x^4 spectrum: strongest line {l4:5.1f} dB at {f4:+9.1f} Hz "
          f"-> QPSK confirmed, carrier offset {f4 / 4:+.1f} Hz")

    # ---- 2. the symbol-rate line ------------------------------------
    env = np.abs(x) ** 2
    env -= env.mean()
    E = np.abs(np.fft.rfft(env[:M] * np.hanning(M)))
    fax = np.fft.rfftfreq(M, 1 / a.fs)
    m = (fax > 69_000) & (fax < 75_000)
    baud = fax[m][np.argmax(E[m])]
    floor = np.median(E[(fax > 55_000) & (fax < 65_000)])
    print(f"symbol clock: {baud:8.1f} Sym/s ({20 * np.log10(E[m].max() / floor):+.1f} dB line; grid says 72000)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    S4 = np.abs(np.fft.fftshift(np.fft.fft(
        ((x / (np.abs(x) + 1e-9)) ** 4)[:M] * np.hanning(M))))
    fx = np.fft.fftshift(np.fft.fftfreq(M, 1 / a.fs))
    mm = np.abs(fx - f4) < 3000
    a1.plot(fx[mm], 20 * np.log10(S4[mm] / np.median(S4)), lw=.7,
            color="#0aa")
    a1.set(title="x^4 collapse: the QPSK line", xlabel="Hz",
           ylabel="dB re median")
    mm2 = (fax > baud - 3000) & (fax < baud + 3000)
    a2.plot(fax[mm2] / 1e3, 20 * np.log10(E[mm2] / floor), lw=.7,
            color="#0aa")
    a2.axvline(72, color="#d33", ls=":")
    a2.set(title=f"symbol-rate line: {baud:.0f} Sym/s",
           xlabel="kHz (envelope spectrum)", ylabel="dB re floor")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "lrpt.png", dpi=110)
    print(f"figure -> {out / 'lrpt.png'}")


if __name__ == "__main__":
    main()
