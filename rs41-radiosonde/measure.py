#!/usr/bin/env python3
"""Measure the Vaisala RS41 radiosonde grid from a 400 MHz IQ capture.

The grid: 4800 Bd GFSK, one 320-byte frame per second. Every frame
begins with the same 8-byte sync — but the ENTIRE frame (sync included)
is whitened by a fixed LFSR mask, so the bits on the air are
  0x86 35 F4 40 93 DF 1A 60
(the whitened form — correlate for THAT, not for the datasheet's
de-whitened bytes; that distinction cost us a day of zero decodes).
After de-whitening, the payload is protected by Reed-Solomon (255,231)
over GF(256) with polynomial 0x11D, with the frame zero-padded to 518
bytes before RS — pad wrong and the GPS coordinates decode to garbage.

This script: FSK-discriminates, recovers the bit clock from the
signal's own transitions, correlates the 64-bit whitened sync, and
reports the measured baud rate + frame period + sync hits.

Usage: python measure.py --iq capture.cs16 --fs 250000
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

HERE = Path(__file__).resolve().parent
# The 64 sync bits AS TRANSMITTED (whitened, on-air bit order). Note
# this is NOT the datasheet hex rendered MSB-first — byte bit-order
# transforms bite hard here. Constants must be proven against live
# signal, not derived on paper (we did both wrong before doing it
# right; see wxTuna's sonde ledger).
SYNC_ONAIR = ("00001000011011010101001110001000"
              "01000100011010010100100000011111")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--secs", type=float, default=20.0)
    a = ap.parse_args()

    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)

    # find the sonde by its FINGERPRINT: FSK shows two spectral lobes
    # exactly one deviation-pair apart (~4.8 kHz for RS41). Hunting the
    # strongest peak instead grabs whatever interferer is loudest —
    # this script's first two runs proved it twice.
    N = 8192
    seg = x[: len(x) // N * N][: 400 * N].reshape(-1, N)
    psd = np.fft.fftshift(
        (np.abs(np.fft.fft(seg * np.hanning(N), axis=1)) ** 2).mean(0))
    fax = np.fft.fftshift(np.fft.fftfreq(N, 1 / a.fs))
    db = 10 * np.log10(psd + 1e-12)
    db -= np.median(db)
    c = N // 2
    db[c - 2:c + 3] = 0                     # ignore DC spike
    binw = a.fs / N
    pk_i = [i for i in range(1, N - 1)
            if db[i] > 5 and db[i] >= db[i - 1] and db[i] >= db[i + 1]
            and abs(fax[i]) < 60_000]
    best = None
    for i in pk_i:
        for j in pk_i:
            if j <= i:
                continue
            if abs((j - i) * binw - 4800) < 700:
                if best is None or db[i] + db[j] > best[2]:
                    best = (i, j, db[i] + db[j])
    if best is None:
        raise SystemExit("no FSK twin-lobe pair found — is a sonde in "
                         "this capture?")
    f0 = float((fax[best[0]] + fax[best[1]]) / 2)
    print(f"FSK pair found: lobes {fax[best[0]]:+.0f} / "
          f"{fax[best[1]]:+.0f} Hz -> carrier {f0:+.0f} Hz from center")
    x = x * np.exp(-2j * np.pi * f0 * np.arange(len(x)) / a.fs)
    x = lfilter(firwin(63, 12_000, fs=a.fs), 1.0, x)
    d = np.angle(x[1:] * np.conj(x[:-1]))
    d -= np.median(d)                       # residual offset

    # baud from the transition spectrum
    tr = np.abs(np.diff(np.sign(d)))
    T = np.abs(np.fft.rfft(tr - tr.mean()))
    fax = np.fft.rfftfreq(len(tr), 1 / a.fs)
    m = (fax > 4000) & (fax < 5600)
    baud = fax[m][np.argmax(T[m])]
    print(f"symbol rate:  {baud:8.1f} Bd   (grid says 4800)")

    # bit-slice at the measured baud and correlate the whitened sync
    sps = a.fs / baud
    idx = (np.arange(int(len(d) / sps) - 2) * sps + sps / 2).astype(int)
    bits = (d[idx] > 0).astype(np.int8)
    pat = np.array([int(b) for b in SYNC_ONAIR], np.int8)
    bb, pp = 2 * bits - 1, 2 * pat - 1
    corr = np.correlate(bb.astype(np.float32), pp.astype(np.float32))
    # FSK mark/space sense depends on your tuner's sideband sense —
    # a fully inverted correlation is the same lock, mirrored
    if np.abs(corr.min()) > corr.max():
        corr = -corr
        print("(polarity inverted — bits flipped, same grid)")
    th = 58                                  # >=58/64 bits match
    hits = np.where(corr >= 2 * th - 64)[0]
    if len(hits):
        hits = hits[np.insert(np.diff(hits) > 100, 0, True)]
    print(f"sync hits:    {len(hits)} frames (64-bit whitened sync, "
          f">= {th}/64 bits)")
    if len(hits) > 2:
        gaps = np.diff(hits) / baud
        print(f"frame period: {np.median(gaps):8.3f} s  "
              f"(one 320-byte frame per second; inter-frame gap is NOT "
              f"an integer bit count — resync per frame!)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    a1.plot(np.arange(len(corr))[:int(6 * baud)] / baud,
            (corr[:int(6 * baud)] + 64) / 128, lw=.5, color="#0aa")
    a1.set(title="whitened-sync correlation (peaks = frame starts)",
           xlabel="time (bit periods -> ~seconds/4800)",
           ylabel="fraction of 64 bits matching")
    a1.axhline(th / 64, color="#d33", ls=":", lw=1)
    seg = int(2 * sps)
    off = int(hits[0] * sps) if len(hits) else 0
    tr2 = d[off:off + 60 * seg]
    for i in range(0, len(tr2) - seg, seg):
        a2.plot(np.arange(seg) / sps, tr2[i:i + seg], color="#0aa",
                alpha=.2, lw=.7)
    a2.set(title=f"GFSK eye at measured {baud:.0f} Bd",
           xlabel="bit periods")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "rs41.png", dpi=110)
    print(f"figure -> {out / 'rs41.png'}")


if __name__ == "__main__":
    main()
