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
    reg = 0
    for b in bits112[:88]:
        reg = ((reg << 1) | int(b)) & 0xFFFFFF
        if reg & 0x800000:
            reg ^= CRC_POLY & 0x7FFFFF
    # standard Mode S: remainder over the first 88 bits vs parity field
    r = 0
    data = int("".join(str(int(b)) for b in bits112), 2)
    for i in range(112 - 24):
        if data & (1 << (111 - i)):
            data ^= CRC_POLY << (88 - 1 - i)
    return (data & 0xFFFFFF) == 0


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
    th = corr.mean() + 6 * corr.std()
    cand = np.where(corr > th)[0]
    cand = cand[np.insert(np.diff(cand) > 2 * spb, 0, True)]

    ok = 0
    keep = None
    for c in cand[:20000]:
        s0 = c + int(8 * spb)
        need = int(112 * spb) + hp
        if s0 + need > len(env):
            break
        first = env[s0 + (np.arange(112) * spb).astype(int)]
        second = env[s0 + (np.arange(112) * spb).astype(int) + hp]
        bits = (first > second).astype(np.int8)
        if crc24(bits):
            ok += 1
            if keep is None:
                keep = (c, bits)
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
