#!/usr/bin/env python3
"""Measure LoRa symbol lock: split timing from frequency with the preamble's own geometry.

A LoRa preamble is 8 up-chirps, 2 sync symbols, then 2.25 DOWN-chirps. De-chirping a window
that starts tau chips off, from a transmitter cfo bins off, puts the peak at
    up-chirp:    fu = cfo - tau
    down-chirp:  fd = cfo + tau
Two unknowns, two measurements:  cfo = (fu + fd) / 2,  tau = (fd - fu) / 2.
This script does only that: find the preamble, read fu and fd, solve, re-cut the stream with the
answer, and show that both peaks now sit at bin 0 and the sync word lands on integers.

Usage: python measure.py --iq cap.cs16 --fs 1e6 --offset-khz -200 --sf 7 --bw 125e3
       [--sweep CFO+5k.cs16:5 CFO-10k.cs16:-10 ...]   captures with a known transmitter detune (kHz)
Output: figures/lora.png + printed bins/chips/Hz.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.signal import resample_poly

HERE = Path(__file__).resolve().parent


def chirps(sf):
    N = 1 << sf
    k = np.arange(N)
    up = np.exp(1j * np.pi * (k * k / N - k)).astype(np.complex64)   # -bw/2 .. +bw/2 in N chips
    return N, up, np.conj(up)


def dechirp(seg, ref, zoom=1):
    return np.abs(np.fft.fft(seg * ref, len(ref) * zoom))


def fine_peak(seg, ref, N, zoom=16):
    """Signed fractional bin (-N/2, N/2] of the dechirp peak, parabolic interpolation, + sharpness."""
    f = dechirp(seg, ref, zoom)
    k = int(f.argmax()); M = len(f)
    l, c, r = f[(k - 1) % M], f[k], f[(k + 1) % M]
    den = l - 2 * c + r
    b = (k + (0.5 * (l - r) / den if den != 0 else 0.0)) / zoom
    return (b + N / 2) % N - N / 2, float(c / (f.mean() + 1e-12))


def wrap(v, N):
    return (v + N / 2) % N - N / 2


def sync_symbols(syncword, N):
    def sym(nib):
        v = ((nib * 8 + 64) % 128) - 64          # nibble x 8, signed (measured on air, see numpy-lora)
        return v % N
    return sym(syncword >> 4), sym(syncword & 0xF)


def load(path, fs, secs):
    raw = np.fromfile(path, dtype=np.int16, count=int(2 * secs * fs) if secs else -1)
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)).astype(np.complex64)


def cut(x, fs, offset_hz, dec, k0, nsamp, shift=0.0):
    """Full-rate samples from k0 (+shift) shifted to 0 Hz and decimated to one sample per chip."""
    k0 = int(round(k0 + shift)); k1 = k0 + nsamp * dec
    seg = x[max(k0, 0):k1]
    n = np.arange(len(seg), dtype=np.float64) + max(k0, 0)
    seg = seg * np.exp(-2j * np.pi * offset_hz * n / fs)
    return resample_poly(seg, 1, dec).astype(np.complex64)


def find_preamble(y, N, dn, min_up=5, thresh=3.0):
    """First window beginning a run of min_up up-chirps whose dechirp bins agree within +-1.
    Noise almost never repeats a bin five times in a row — bin agreement does the detecting."""
    nsym = len(y) // N
    W = y[:nsym * N].reshape(nsym, N) * dn
    F = np.abs(np.fft.fft(W, axis=1))
    bins = F.argmax(axis=1); mags = F.max(axis=1) / (F.mean(axis=1) + 1e-12)
    for i in range(nsym - min_up):
        if mags[i] <= thresh: continue
        d = (bins[i:i + min_up] - bins[i] + N // 2) % N - N // 2
        if np.all(np.abs(d) <= 1) and np.all(mags[i:i + min_up] > thresh):
            return i * N
    return None


def lock(x, fs, offset_hz, sf, bw, syncword, start_chip, verbose=True):
    """Returns dict(fu, fd, cfo, tau, residual up/down, sync symbols before/after, spectra)."""
    N, up, dn = chirps(sf)
    dec = int(round(fs / bw))
    s0, s1 = sync_symbols(syncword, N)
    back = 2 * N
    y = cut(x, fs, offset_hz, dec, (start_chip - back) * dec, 24 * N)
    pb, _ = fine_peak(y[back:back + N], dn, N)
    # --- walk in 1/8-symbol steps until the two sync symbols appear relative to the preamble bin
    q = N // 8; cands = []
    for w in range(96):                        # up to 12 symbols of preamble (16-symbol preambles exist)
        o = back + w * q
        if o + 4 * N > len(y): break
        a, ma = fine_peak(y[o:o + N], dn, N); b, mb = fine_peak(y[o + N:o + 2 * N], dn, N)
        ref = pb - (o - back)                  # moving the window by d chips moves every bin by -d
        ra, rb = (a - ref) % N, (b - ref) % N
        if min(abs(ra - s0), N - abs(ra - s0)) < 2.5 and min(abs(rb - s1), N - abs(rb - s1)) < 2.5:
            cands.append((ma + mb, o))
        elif cands:
            break
    if not cands:
        raise RuntimeError("sync word not found")
    ks = max(cands)[1]
    # --- THE MEASUREMENT: up-chirp bin from the preamble, down-chirp bin 2 windows after sync
    fu = float(np.median([fine_peak(y[ks - k * N:ks - (k - 1) * N], dn, N)[0] for k in (1, 2, 3, 4)]))
    fd = float(np.mean([fine_peak(y[ks + k * N:ks + (k + 1) * N], up, N)[0] for k in (2, 3)]))
    sync_before = [fine_peak(y[ks + j * N:ks + (j + 1) * N], dn, N)[0] % N for j in (0, 1)]
    spec_up_before = dechirp(y[ks - N:ks], dn); spec_dn_before = dechirp(y[ks + 2 * N:ks + 3 * N], up)
    # --- solve; the sum is ambiguous by N/2 (bins wrap), the sign of tau by symmetry: try both,
    #     keep whichever puts the two preamble bins, the sync word and the down-chirp where they belong
    first_chip = start_chip - back + ks - 6 * N
    best = None
    for c_off in (0.0, N / 2):
        cfo = wrap((fu + fd) / 2 + c_off, N); tau = wrap((fd - fu) / 2 + c_off, N)
        if abs(cfo) > N / 4 - 1: continue                     # LoRa chips specify +-25 % of BW
        for sgn in (+1, -1):
            y2 = cut(x, fs, offset_hz, dec, first_chip * dec, 13 * N, shift=sgn * tau * dec)
            k = np.arange(len(y2)); y2 = y2 * np.exp(-2j * np.pi * cfo * k / N).astype(np.complex64)
            e = 0.0
            for j, want in ((-2, 0), (-1, 0), (0, s0), (1, s1)):
                v, _ = fine_peak(y2[6 * N + j * N:6 * N + (j + 1) * N], dn, N); e += abs(wrap(v - want, N))
            d, _ = fine_peak(y2[8 * N:9 * N], up, N); e += abs(wrap(d, N))
            if best is None or e < best[0]: best = (e, cfo, sgn * tau, y2)
    e, cfo, tau, y2 = best
    ru = float(np.mean([fine_peak(y2[6 * N + j * N:6 * N + (j + 1) * N], dn, N)[0] for j in (-4, -3, -2, -1)]))
    rd = float(np.mean([fine_peak(y2[6 * N + j * N:6 * N + (j + 1) * N], up, N)[0] for j in (2, 3)]))
    sync_after = [fine_peak(y2[6 * N + j * N:6 * N + (j + 1) * N], dn, N)[0] % N for j in (0, 1)]
    hz_per_bin = bw / N
    if verbose:
        print(f"  up-chirp bin  fu = {fu:+7.2f}   down-chirp bin fd = {fd:+7.2f}")
        print(f"  cfo = (fu+fd)/2 = {cfo:+7.2f} bins = {cfo * hz_per_bin:+8.1f} Hz")
        print(f"  tau = (fd-fu)/2 = {tau:+7.2f} chips = {tau / bw * 1e6:+8.1f} us")
        print(f"  after re-cut: up {ru:+.3f}  down {rd:+.3f} bins   sync {sync_before[0]:.2f},{sync_before[1]:.2f}"
              f" -> {sync_after[0]:.2f},{sync_after[1]:.2f}  (want {s0},{s1})")
    # lock curve: up-chirp sharpness vs window offset around the found alignment
    offs = np.arange(-N // 2, N // 2 + 1, 2); sharp = []
    for o in offs:
        y3 = cut(x, fs, offset_hz, dec, first_chip * dec, 6 * N, shift=(tau + o) * dec)
        k = np.arange(len(y3)); y3 = y3 * np.exp(-2j * np.pi * cfo * k / N)
        sharp.append(fine_peak(y3[2 * N:3 * N], dn, N)[1])
    return dict(fu=fu, fd=fd, cfo=cfo, tau=tau, ru=ru, rd=rd, s_before=sync_before, s_after=sync_after,
                want=(s0, s1), up_b=spec_up_before, dn_b=spec_dn_before,
                up_a=dechirp(y2[5 * N:6 * N], dn), dn_a=dechirp(y2[8 * N:9 * N], up),
                offs=offs, sharp=np.array(sharp), hz_per_bin=hz_per_bin, N=N, err=e)


def first_frame(path, fs, offset_hz, sf, bw, syncword, secs, verbose=True):
    N, up, dn = chirps(sf)
    dec = int(round(fs / bw))
    x = load(path, fs, secs)
    # coarse scan at one sample per chip, in 2 s slabs
    hop = int(2 * bw); n = 0
    while n * N < len(x) // dec - 8 * N:
        y = cut(x, fs, offset_hz, dec, n * N * dec, hop)
        p = find_preamble(y, N, dn)
        if p is not None:
            try:
                return lock(x, fs, offset_hz, sf, bw, syncword, n * N + p // N * N, verbose)
            except RuntimeError:
                pass
        n += hop // N - 8
    raise RuntimeError(f"no LoRa preamble in {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--offset-khz", type=float, default=0.0, help="LoRa channel minus capture centre")
    ap.add_argument("--sf", type=int, default=7)
    ap.add_argument("--bw", type=float, default=125e3)
    ap.add_argument("--sync", type=lambda s: int(s, 0), default=0x34)
    ap.add_argument("--secs", type=float, default=0.0, help="seconds to read (0 = all)")
    ap.add_argument("--sweep", nargs="*", default=[], help="path:detune_khz captures, same SF/BW/offset")
    a = ap.parse_args()
    off = a.offset_khz * 1e3

    print(f"{a.iq}: SF{a.sf} BW{a.bw / 1e3:g}k, {a.bw / (1 << a.sf):.1f} Hz per bin, {1e6 / a.bw:.1f} us per chip")
    r = first_frame(a.iq, a.fs, off, a.sf, a.bw, a.sync, a.secs)

    sweep = []
    for item in a.sweep:
        path, khz = item.rsplit(":", 1)
        print(f"{path}: transmitter detuned {float(khz):+g} kHz")
        try:
            s = first_frame(path, a.fs, off, a.sf, a.bw, a.sync, a.secs)
            sweep.append((float(khz), s["cfo"] * s["hz_per_bin"] / 1e3))
        except RuntimeError as ex:
            print(f"  {ex}")
    sweep.append((0.0, r["cfo"] * r["hz_per_bin"] / 1e3))
    if len(sweep) > 1:
        d = np.array(sweep)
        print("\ncommanded detune vs measured cfo (kHz):")
        for k, m in sorted(sweep): print(f"  {k:+6.1f}  ->  {m:+7.2f}   diff {m - k:+.2f}")
        print(f"  slope {np.polyfit(d[:, 0], d[:, 1], 1)[0]:.4f} (1.0000 = the receiver reads the detune exactly)")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    N = r["N"]; b = wrap(np.arange(N), N); order = np.argsort(b)
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    for a_, before, after, name, bin_ in ((ax[0, 0], r["up_b"], r["up_a"], "up-chirp (preamble)", r["fu"]),
                                          (ax[0, 1], r["dn_b"], r["dn_a"], "down-chirp (after sync)", r["fd"])):
        a_.plot(b[order], before[order] / before.max(), "--", color="0.5", label="raw window")
        a_.plot(b[order], after[order] / after.max(), color="C0", label="after lock")
        a_.axvline(bin_, color="C3", lw=0.8); a_.axvline(0, color="C0", lw=0.8)
        a_.set_title(f"{name}: peak {bin_:+.2f} → {(r['ru'] if 'up' in name else r['rd']):+.3f} bins")
        a_.set_xlabel("de-chirped FFT bin (signed)"); a_.legend(); a_.set_xlim(-N / 2, N / 2)
    ax[1, 0].plot(r["offs"], r["sharp"], ".-")
    ax[1, 0].axvline(0, color="C3", lw=0.8)
    ax[1, 0].set_title(f"preamble sharpness barely cares about timing (τ = {r['tau']:+.2f} chips)")
    ax[1, 0].set_xlabel("chips from the solved alignment"); ax[1, 0].set_ylabel("peak / mean")
    if len(sweep) > 1:
        d = np.array(sorted(sweep))
        ax[1, 1].plot(d[:, 0], d[:, 1], "o", label="measured cfo")
        lim = max(abs(d).max() * 1.2, 1)
        ax[1, 1].plot([-lim, lim], [-lim, lim], "--", color="0.6", label="y = x")
        ax[1, 1].set_xlabel("transmitter detune commanded (kHz)"); ax[1, 1].set_ylabel("cfo read from the preamble (kHz)")
        ax[1, 1].set_title("frequency offset: measured vs commanded"); ax[1, 1].legend()
    else:
        ax[1, 1].axis("off")
    fig.suptitle(f"LoRa SF{a.sf} BW{a.bw / 1e3:g}k: fu = cfo − τ, fd = cfo + τ   →   "
                 f"cfo {r['cfo'] * r['hz_per_bin']:+.0f} Hz, τ {r['tau'] / a.bw * 1e6:+.1f} µs")
    fig.tight_layout()
    (HERE / "figures").mkdir(exist_ok=True)
    fig.savefig(HERE / "figures" / "lora.png", dpi=120)
    print(f"\nwrote {HERE / 'figures' / 'lora.png'}")


if __name__ == "__main__":
    main()
