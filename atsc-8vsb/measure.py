#!/usr/bin/env python3
"""Measure the ATSC 1.0 (8-VSB) grid from an IQ capture of a TV channel.

The grid:
  symbol rate   10.762238 MSym/s  (= 4.5 MHz x 684/286)
  segment       832 symbols (4 sync + 828 data)  -> 12,935.4 seg/s
  field         313 segments (1 field sync + 312 data) -> 41.32 fields/s
  pilot         a real carrier remnant 309.441 kHz above the lower
                channel edge (the receiver's first handhold)
  field sync    PN511 + 3x PN63 training sequences (the equalizer's
                free lunch, sent 41 times a second)

The trick of this script: you do NOT need to sample at 10.76+ MS/s to
prove the grid. The 4-symbol segment sync repeats at exactly
symbol_rate/832 = 12,935.4 Hz, which stamps a cyclostationary line at
that frequency into the squared envelope |x|^2 — visible at any
capture rate. Find the line, multiply by 832, and you've measured the
symbol clock to ppm precision without ever slicing a symbol.

Usage: python measure.py --iq capture.cs16 --fs 8000000
       (IQ centered mid-channel; any rate >= ~6 MS/s)
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SEG_HZ = 10_762_238.0 / 832          # 12,935.38 Hz
FIELD_HZ = SEG_HZ / 313              # 41.33 Hz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--secs", type=float, default=4.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)

    # ---- 1. the VSB spectrum + pilot --------------------------------
    N = 1 << 15
    seg = x[: len(x) // N * N].reshape(-1, N)
    psd = np.fft.fftshift(
        (np.abs(np.fft.fft(seg * np.hanning(N), axis=1)) ** 2).mean(0))
    fax = np.fft.fftshift(np.fft.fftfreq(N, 1 / a.fs))
    db = 10 * np.log10(psd + 1e-9)
    shelf = np.median(db[np.abs(fax) < 2e6])
    pk = int(np.argmax(db - (np.abs(fax) > 3.2e6) * 100))
    pilot_off = fax[pk]
    print(f"pilot tone:   {pilot_off / 1e6:+.4f} MHz from capture center, "
          f"{db[pk] - shelf:+.1f} dB above the data shelf")
    print(f"              (grid: 309.441 kHz above the lower channel "
          f"edge = -2.6906 MHz from channel center)")

    # ---- 2. fold the envelope at the segment period -----------------
    # The 4-symbol segment sync recurs every 832 symbols. Fold the
    # envelope at a candidate period and average: at the WRONG period
    # the sync smears to nothing; at the right one (618.46 samples at
    # 8 MS/s) it stands up out of the noise. The sharpest fold IS the
    # measurement of the symbol clock (period x 832 / fs).
    env = np.abs(x[: int(3.0 * a.fs)]) ** 2
    env -= env.mean()
    P0 = a.fs / SEG_HZ
    n_seg = int(len(env) / (P0 * 1.01)) - 2
    n_seg = min(n_seg, 3000)
    j = np.arange(int(P0))
    cands = P0 * (1 + np.linspace(-3e-3, 3e-3, 121))
    sharp = np.empty(len(cands))
    profs = []
    base = np.arange(n_seg)[:, None]
    for ci, P in enumerate(cands):
        pos = base * P + j[None, :]
        i0 = pos.astype(np.int64)
        fr = pos - i0
        prof = (env[i0] * (1 - fr) + env[i0 + 1] * fr).mean(axis=0)
        prof -= np.median(prof)
        sharp[ci] = prof.max() / (prof.std() + 1e-12)
        profs.append(prof)
    bi = int(np.argmax(sharp))
    P_best = cands[bi]
    seg_hz = a.fs / P_best
    sym_rate = seg_hz * 832
    print(f"fold sweep:   sharpest at period {P_best:.3f} samples "
          f"(sharpness {sharp[bi]:.1f} vs {np.median(sharp):.1f} off-grid)")
    print(f"segment rate: {seg_hz:10.2f} Hz   (grid says {SEG_HZ:.2f})")
    print(f"symbol clock: {sym_rate / 1e6:10.6f} MSym/s = 832 x segment "
          f"rate (grid says 10.762238)")
    print(f"field rate:   {seg_hz / 313:10.3f} Hz = segment/313 "
          f"(grid says {FIELD_HZ:.3f}; PN511 training rides here)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    a1.plot(fax / 1e6, db, lw=.5, color="#0aa")
    a1.annotate("pilot\n(the handhold)", (pilot_off / 1e6, db[pk]),
                fontsize=9, color="#d33", ha="right")
    a1.annotate("8-VSB data shelf\n(10.76 MSym/s of it)",
                (0.5, shelf + 4), fontsize=9)
    a1.set(title="8-VSB channel: flat data shelf + pilot",
           xlabel="MHz from capture center", ylabel="PSD (dB)")
    prof = profs[bi]
    off_p = profs[10]
    t_us = j / a.fs * 1e6
    a2.plot(t_us, off_p / (np.abs(prof).max() + 1e-12), lw=.7,
            color="#888", label="folded 0.25% off-period: nothing")
    a2.plot(t_us, prof / (np.abs(prof).max() + 1e-12), lw=.9,
            color="#0aa", label=f"folded at {P_best:.2f} samples: "
                                f"the segment sync appears")
    a2.legend(fontsize=8, loc="upper right")
    a2.set(title=f"envelope folded at the segment period "
                 f"({n_seg} segments averaged)",
           xlabel="time within segment (us)", ylabel="normalized")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "atsc_grid.png", dpi=110)
    print(f"figure -> {out / 'atsc_grid.png'}")


if __name__ == "__main__":
    main()
