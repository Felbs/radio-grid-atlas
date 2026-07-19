#!/usr/bin/env python3
"""Measure the ADS-B / Mode S grid: the 1 microsecond pulse ruler.

Finds preambles by matched correlation against the 0/1.0/3.5/4.5 us
pulse pattern, PPM-slices the following 112 bits, verifies CRC-24
(poly 0xFFF409), and plots the pulse timeline of a verified frame.
The CRC pass count is the honest detection statistic — impulse noise
fakes preambles at 1090 MHz constantly.

Usage: python measure.py --iq capture.cs16 --fs 2400000
       (30 s centered on 1090 MHz; cu8 rtl-style also fine via --cu8)
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CRC_POLY = 0xFFF409


def crc24(bits112):
    """Mode S parity: shift-register division of the whole frame by the
    FULL 25-bit generator 0x1FFF409 — the x^24 term is not optional.
    (Writing the generator as '0xFFF409' and hand-aligning shifts cost
    this script two wrong versions; the leading term is what cancels
    the bit you aligned to. Remainder 0 == valid DF17/18.)"""
    reg = 0
    for b in bits112:
        reg = (reg << 1) | int(b)
        if reg & (1 << 24):
            reg ^= 0x1FFF409
    return (reg & 0xFFFFFF) == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--cu8", action="store_true")
    ap.add_argument("--secs", type=float, default=30.0)
    a = ap.parse_args()

    if a.cu8:
        raw = np.fromfile(a.iq, dtype=np.uint8,
                          count=int(2 * a.secs * a.fs)).astype(np.float32)
        x = (raw[0::2] - 127.5) + 1j * (raw[1::2] - 127.5)
    else:
        raw = np.fromfile(a.iq, dtype=np.int16,
                          count=int(2 * a.secs * a.fs))
        x = raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)
    env = np.abs(x)
    spb = a.fs / 1e6                 # samples per microsecond
    hp = int(round(spb / 2))         # half-microsecond pulse

    # matched preamble: +1 where pulses live, -1 in the quiet
    tpl = np.full(int(9 * spb), -0.5, np.float32)
    for t0 in (0.0, 1.0, 3.5, 4.5):
        tpl[int(t0 * spb):int(t0 * spb) + hp] = 1.0
    corr = np.correlate(env - env.mean(), tpl, mode="valid")
    # robust threshold: median + k*MAD, NOT mean + k*std — on a short
    # or busy capture the frames themselves inflate the std until the
    # threshold climbs above its own peaks (found via synthetic TX)
    med = np.median(corr)
    mad = np.median(np.abs(corr - med)) + 1e-9
    th = med + 12 * mad
    cand_all = np.where(corr > th)[0]
    cand = []
    if len(cand_all):
        # keep each cluster's correlation ARGMAX, not its first sample —
        # a noise shoulder can start a cluster microseconds before the
        # true peak and push the frame outside the alignment scan
        # (found by differential debug against a reference decoder:
        # slicer read 112/112 at truth, detector never offered truth)
        splits = np.where(np.diff(cand_all) > 2 * spb)[0] + 1
        for grp in np.split(cand_all, splits):
            cand.append(int(grp[np.argmax(corr[grp])]))
    cand = np.array(cand, dtype=np.int64)

    # PPM slicing that works at ANY sample rate: integrate envelope
    # energy over each half-microsecond window via a cumulative sum
    # with fractional edges (at 2.4 MS/s a half-bit is 1.2 samples —
    # single-sample slicing straddles edges and fails every CRC).
    F = np.concatenate(([0.0], np.cumsum(env, dtype=np.float64)))

    def F_at(pos):
        i = pos.astype(np.int64)
        fr = pos - i
        return F[i] * (1 - fr) + F[i + 1] * fr

    def window_energy(a, b):
        return F_at(b) - F_at(a)

    ok = 0
    keep = None
    bit_starts = np.arange(112) * spb
    # correlation peaks land within ~1 us of true preamble start; scan
    # sub-microsecond alignments and let the CRC pick the right one
    subs = np.arange(-1.0, 1.01, 0.25) * spb
    for c in cand[:20000]:
        if c + 10 * spb + 113 * spb > len(env):
            break
        for off in subs:
            a = c + off + 8 * spb + bit_starts
            e1 = window_energy(a, a + spb / 2)
            e2 = window_energy(a + spb / 2, a + spb)
            bits = (e1 > e2).astype(np.int8)
            if crc24(bits):
                ok += 1
                if keep is None:
                    keep = (int(c + off), bits)
                break
    print(f"preamble candidates: {len(cand)}   CRC-24 verified frames: {ok}")
    if keep is None:
        print("no verified frames — busier sky or better antenna needed")
        return
    c, bits = keep
    icao = int("".join(str(b) for b in bits[8:32]), 2)
    print(f"first verified frame: DF={int(''.join(str(b) for b in bits[:5]), 2)} "
          f"ICAO {icao:06X}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 4))
    n0 = int(c - 2 * spb)
    n1 = int(c + 30 * spb)
    t = (np.arange(n0, n1) - c) / spb
    ax.plot(t, env[n0:n1], lw=.8, color="#0aa")
    for t0 in (0.0, 1.0, 3.5, 4.5):
        ax.axvline(t0, color="#d33", ls=":", lw=1)
    ax.annotate("preamble pulses at 0, 1.0, 3.5, 4.5 us", (0, env[n0:n1].max()),
                fontsize=9, color="#d33")
    ax.annotate("112 PPM bits (1 us each) ->", (9, env[n0:n1].max() * .8),
                fontsize=9)
    ax.set(title=f"a CRC-verified extended squitter off the air "
                 f"(ICAO {icao:06X})", xlabel="microseconds from preamble",
           ylabel="|envelope|")
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "adsb.png", dpi=110)
    print(f"figure -> {out / 'adsb.png'}")


if __name__ == "__main__":
    main()
