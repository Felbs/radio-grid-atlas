#!/usr/bin/env python3
"""Measure the FM stereo composite grid from a raw IQ capture.

The grid: everything in an FM station's baseband hangs off ONE tone —
the 19 kHz pilot. The stereo subcarrier is exactly 2x it (38 kHz), RDS
is exactly 3x it (57 kHz). This script demodulates the composite and
measures all three relationships, plus recovers actual stereo audio via
pilot phase-doubling, printing the separation it achieves.

Usage:
  python measure.py --iq capture.cs16 --fs 2976750 [--f-offset 0]

Input: interleaved int16 IQ centered on the station. Any sample rate
>= 500 kHz works. Output: printed measurements + figures/composite.png
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter, welch

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--secs", type=float, default=8.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)

    # channel-select: +-120 kHz (Carson BW for +-75 kHz deviation FM),
    # then decimate so the discriminator only sees ONE station
    dec = max(1, int(a.fs / 300_000))
    fsd = a.fs / dec
    x = lfilter(firwin(91, 120_000, fs=a.fs), 1.0, x)[::dec]

    # FM discriminator: instantaneous frequency = the composite baseband
    c = np.angle(x[1:] * np.conj(x[:-1])) * fsd / (2 * np.pi)  # in Hz

    # composite spectrum
    f, psd = welch(c, fs=fsd, nperseg=1 << 14)
    db = 10 * np.log10(psd + 1e-12)

    def peak_near(f0, half=400):
        m = (f > f0 - half) & (f < f0 + half)
        return f[m][np.argmax(db[m])], db[m].max()

    def centroid(f0, half):
        """Modulated subcarriers are BANDS, not tones — a peak-pick
        lands on whichever sideband is loudest right now; the power
        centroid finds the (suppressed) carrier underneath."""
        m = (f > f0 - half) & (f < f0 + half)
        w = psd[m] - psd[m].min()
        return float((f[m] * w).sum() / w.sum()), 10 * np.log10(
            psd[m].max() + 1e-12)

    fp, ap_ = peak_near(19_000)          # the pilot IS a tone
    fs38, a38 = centroid(38_000, 3000)   # DSB-SC band around 2x
    f57, a57 = centroid(57_000, 2600)    # RDS band around 3x
    floor = np.median(db[(f > 20_000) & (f < 22_000)])

    print(f"pilot:        {fp:9.1f} Hz   ({ap_ - floor:+5.1f} dB above "
          f"guard-band floor)")
    print(f"38k region:   {fs38:9.1f} Hz   = {fs38 / fp:.4f} x pilot "
          f"(grid says 2)")
    print(f"57k RDS:      {f57:9.1f} Hz   = {f57 / fp:.4f} x pilot "
          f"(grid says 3)")

    # stereo via pilot phase-doubling: (L-R) lands in the IMAGINARY part
    # of composite * conj(pilot_phasor^2) — sin-phased subcarrier vs sin
    # pilot. (We lost an afternoon to the 90-degree convention; the
    # measurement below settles it on any capture.)
    n = np.arange(len(c))
    p = c * np.exp(-2j * np.pi * fp * n / fsd)
    a1 = float(np.exp(-2 * np.pi * 500.0 / fsd))
    for _ in range(3):
        p = lfilter([1 - a1], [1, -a1], p)
    u = p / (np.abs(p) + 1e-9)
    z = c * np.exp(-2j * np.pi * 2 * fp * n / fsd) * np.conj(u) ** 2
    b15 = firwin(305, 15_000, fs=fsd)
    zf = lfilter(b15, 1.0, z)[len(b15):]
    m = lfilter(b15, 1.0, c)[len(b15):]
    s_re = 10 * np.log10((zf.real ** 2).mean() / (m ** 2).mean())
    s_im = 10 * np.log10((zf.imag ** 2).mean() / (m ** 2).mean())
    print(f"L-R energy:   Re branch {s_re:+5.1f} dB rel M | Im branch "
          f"{s_im:+5.1f} dB rel M  -> stereo lives in "
          f"{'Im' if s_im > s_re else 'Re'}")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 5))
    m_ = f < 80_000
    ax.plot(f[m_] / 1e3, db[m_], lw=0.7, color="#0aa")
    for freq, label in ((fp, "19 kHz pilot\n(the clock)"),
                        (2 * fp, "38 kHz = 2x pilot\n(L-R) DSB-SC"),
                        (3 * fp, "57 kHz = 3x pilot\nRDS")):
        ax.axvline(freq / 1e3, color="#d33", ls=":", lw=1)
        ax.annotate(label, (freq / 1e3, db[m_].max() - 3),
                    ha="center", fontsize=9, color="#d33")
    ax.annotate("L+R mono audio", (7, db[m_].max() - 8), fontsize=9)
    ax.set_xlabel("composite frequency (kHz)")
    ax.set_ylabel("PSD (dB)")
    ax.set_title("FM stereo composite, demodulated off-air: "
                 "one pilot, everything else locked to it")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "composite.png", dpi=110)
    print(f"figure -> {out / 'composite.png'}")


if __name__ == "__main__":
    main()
