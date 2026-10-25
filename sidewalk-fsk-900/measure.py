#!/usr/bin/env python3
"""Measure the Amazon Sidewalk SubG-FSK beacon grid from raw IQ.

Every Sidewalk gateway (an Echo, a Ring bridge) sends a Beacon every 10.08 s on
a hopping channel, with a preamble long enough (spec: 140 ms) that an endpoint
scanning 69 channels cannot miss it.  The PHY is IEEE 802.15.4g SUN FSK:

  GFSK 50 kbit/s, deviation 25 kHz (h = 1), channels 902.2 + 0.2 n MHz (n = 0..68)
  preamble 0x55.., SFD 0x904E, PHR (mode switch, FCS type, whitening, 11-bit length)
  PN9 data whitening, CRC-16 ITU-T (0x1021, init 0) over the transmitted bits

This script derives, from the IQ alone:
  1. the bit rate, from the run lengths of the discriminator inside the preamble
  2. the deviation, from the two tone levels
  3. the preamble duration and the SFD position, per burst
  4. the PHR fields; the PN9 whitening (found as the sequence phase that zeroes
     the CRC); CRC-verified frames
  5. the channel grid: every verified burst's centre against 902.2 + 0.2 n MHz
  6. the beacon period: gaps between verified beacons that carry the same
     (masked) gateway ID, against 10.08 s; and the 63 ms slot grid
PRIVACY: gateway IDs are masked to two bytes' worth of hash; nothing about which
beacon is ours is printed or kept.  Beacon payloads are unencrypted by design
(cloud status, hop seed, sequence number) and are not printed either.

Usage: python measure.py --iq sidewalk_907p4_2000k.cs16 --fs 2000000 --fc 907.4e6 [--png figures/sidewalk.png]
"""
import argparse
import collections
import hashlib
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

RATE = 50e3
SFD = "1001000001001110"           # 0x904E as transmitted
CH0, CH_STEP = 902.2e6, 200e3
BEACON_S, SLOT_S = 10.080, 0.063


def crc16_bits(bitstr, poly=0x1021, init=0):
    c = init
    for ch in bitstr:
        fb = ((c >> 15) & 1) ^ int(ch)
        c = ((c << 1) & 0xFFFF) ^ (poly if fb else 0)
    return c


def pn9():
    st, out = 0x1FF, []
    for _ in range(511):
        bit = st & 1
        fb = ((st >> 5) ^ st) & 1
        st = (st >> 1) | (fb << 8)
        out.append(bit)
    return out


PN9 = pn9() * 3


def burst_scan(path, fs, win_s=1e-3, min_s=0.1, blk=4_000_000):
    """Wideband envelope at 1 ms -> bursts longer than min_s (preamble + frame is ~155 ms)."""
    m = int(win_s * fs)
    mm = np.memmap(path, np.int16, mode="r")
    env = []
    for a in range(0, mm.size, 2 * blk):
        raw = np.asarray(mm[a:a + 2 * blk], dtype=np.float32)
        p = raw[0::2] ** 2 + raw[1::2] ** 2
        env.append(np.sqrt(p[:p.size // m * m].reshape(-1, m).mean(1)))
    e = np.concatenate(env)
    thr = np.median(e) * 3
    idx = np.flatnonzero(e > thr)
    segs = np.split(idx, np.flatnonzero(np.diff(idx) > 2) + 1) if idx.size else []
    return [(s[0] * win_s, (s[-1] - s[0] + 1) * win_s) for s in segs if (s[-1] - s[0] + 1) * win_s >= min_s], mm.size // 2 / fs


def load(path, fs, t0, dur):
    raw = np.fromfile(path, dtype=np.int16, offset=max(int(t0 * fs), 0) * 4, count=int(dur * fs) * 2)
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0


def centres(x, fs, fc):
    """A long FSK burst may overlap another on a different channel: return the channel-grid centres with power."""
    n = 1 << 14
    F = np.abs(np.fft.fftshift(np.fft.fft(x[:n] * np.hanning(n), n))) ** 2
    fax = np.fft.fftshift(np.fft.fftfreq(n, 1 / fs)) + fc
    out = []
    for k in range(69):
        f0 = CH0 + k * CH_STEP
        w = np.abs(fax - f0) < 60e3
        if w.any():
            out.append((float(F[w].sum()), f0))
    out.sort(reverse=True)
    return out


def demod(x, fs, f_off):
    z = lfilter(firwin(129, 60e3, fs=fs), 1.0, x * np.exp(-2j * np.pi * f_off / fs * np.arange(x.size)))
    d = np.angle(z[1:] * np.conj(z[:-1])) * fs / (2 * np.pi)
    p = np.abs(z[1:]) ** 2
    on = p > 0.3 * np.percentile(p, 95)
    return d, on


def analyse(x, fs, f_off):
    d, on = demod(x, fs, f_off)
    if on.sum() < 1000:
        return None
    dd = d[on]
    lo, hi = np.percentile(dd, 10), np.percentile(dd, 90)
    thr = 0.5 * (lo + hi)
    i0, i1 = np.flatnonzero(on)[0], np.flatnonzero(on)[-1]
    lv = (d > thr).astype(np.int8)
    tr = np.flatnonzero(np.diff(lv[i0:i1]) != 0)
    runs = np.diff(tr)
    T0 = fs / RATE
    k = np.round(runs / T0)
    good = (k >= 1) & (np.abs(runs / np.maximum(k, 1) - T0) < 0.25 * T0)
    if good.sum() < 100:
        return None
    T = runs[good].sum() / k[good].sum()
    best = None
    for ph in np.arange(0, T, T / 8):
        ii = (i0 + ph + T / 2 + np.arange(int((i1 - i0 - ph - T) / T)) * T).astype(int)
        q = np.abs(d[ii] - thr).mean()
        if best is None or q > best[0]:
            best = (q, (d[ii] > thr).astype(np.uint8), ph)
    bits = "".join(map(str, best[1]))
    res = dict(T=T, dev=(hi - lo) / 2, n_bits=len(bits), frame=None)
    j = bits.find(SFD)
    if j < 0:
        return res
    res["sfd_sample"] = i0 + best[2] + j * T          # where the sync word starts, in samples of the demodulated burst
    # preamble length: alternating bits before the SFD
    p = j
    while p >= 2 and bits[p - 1] != bits[p - 2]:
        p -= 1
    res["preamble_bits"] = j - p
    phr = bits[j + 16:j + 32]
    if len(phr) < 16:
        return res
    fcs16, whit = int(phr[3]), int(phr[4])
    length = int(phr[5:16], 2)
    k0 = j + 32
    if k0 + 8 * length > len(bits):
        res["frame"] = dict(ok=False, length=length, note="truncated")
        return res
    raw = bits[k0:k0 + 8 * length]
    # whitening: the PN9 phase that zeroes the CRC (searched once, then fixed)
    ok, ph_used, dw = False, None, raw
    for ph in (9,) + tuple(range(511)):
        w = PN9[ph:ph + 8 * length]
        cand = "".join(str(int(a) ^ b) for a, b in zip(raw, w))
        if crc16_bits(cand) == 0:
            ok, ph_used, dw = True, ph, cand
            break
    octets = bytes(int(dw[m:m + 8][::-1], 2) for m in range(0, 8 * length, 8))
    gid = hashlib.sha256(octets[3:7]).hexdigest()[:4] if length >= 7 else "????"
    res["frame"] = dict(ok=ok, length=length, fcs16=fcs16, whit=whit, pn9_phase=ph_used, gid=gid, head=octets[:3].hex())
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--fc", type=float, required=True)
    ap.add_argument("--png")
    a = ap.parse_args()
    fs, fc = a.fs, a.fc
    bursts, total = burst_scan(a.iq, fs)
    print(f"bursts >= 100 ms: {len(bursts)} in {total:.0f} s  ({a.iq})")
    rows = []
    for t0, dur in bursts:
        x = load(a.iq, fs, t0 - 2e-3, dur + 4e-3)
        for pw, f0 in centres(x, fs, fc)[:2]:
            r = analyse(x, fs, f0 - fc)
            if r is None:
                continue
            r.update(t0=t0, dur=dur, f0=f0, ch=int(round((f0 - CH0) / CH_STEP)))
            rows.append(r)
            if r["frame"] and r["frame"]["ok"]:
                break
    ver = [r for r in rows if r["frame"] and r["frame"]["ok"]]
    print(f"FSK bursts analysed: {len(rows)}; with SFD 0x904E: {sum(1 for r in rows if 'preamble_bits' in r)}; CRC-verified frames: {len(ver)}")
    if not ver:
        return
    T = np.array([r["T"] for r in ver]) / fs
    print(f"bit period: {np.median(T)*1e6:.3f} us = {1/np.median(T):.1f} bit/s (spread {np.std(T)*1e6:.3f} us); "
          f"deviation ±{np.median([r['dev'] for r in ver])/1e3:.1f} kHz (spec: 50 kbit/s, 25 kHz, h = 1)")
    pre = np.array([r["preamble_bits"] for r in ver]) * np.median(T)
    print(f"preamble before the SFD: median {np.median(pre)*1e3:.1f} ms, 10-90 % {np.percentile(pre,10)*1e3:.1f}-{np.percentile(pre,90)*1e3:.1f} ms (spec beacon preamble 140 ms)")
    print(f"PHR: lengths {collections.Counter(r['frame']['length'] for r in ver)}, FCS16 flag {collections.Counter(r['frame']['fcs16'] for r in ver)}, "
          f"whitening flag {collections.Counter(r['frame']['whit'] for r in ver)}, PN9 phase that zeroes the CRC {collections.Counter(r['frame']['pn9_phase'] for r in ver)}")
    print(f"first three dewhitened octets: {collections.Counter(r['frame']['head'] for r in ver)}")
    offs = np.array([r["f0"] for r in ver])
    print(f"channels used (n in 902.2 + 0.2 n MHz): {sorted(collections.Counter(r['ch'] for r in ver).items())}")
    byid = collections.defaultdict(list)
    for r in ver:
        byid[r["frame"]["gid"]].append(r["t0"])
    print(f"distinct gateway IDs (hashed): {len(byid)}; seen more than once: {sum(1 for v in byid.values() if len(v) > 1)}")
    gaps = []
    for gid, ts in byid.items():
        ts = sorted(ts)
        for g in np.diff(ts):
            gaps.append(g)
    gaps = np.array(gaps)
    if gaps.size:
        k = np.round(gaps / BEACON_S)
        per = gaps[k >= 1] / k[k >= 1]
        print(f"same-gateway gaps: {np.round(gaps, 3).tolist()} s -> per beacon {np.median(per):.3f} s (spec 10.080 s), "
              f"residual vs 10.080 s: {np.round((gaps - k*BEACON_S)*1e3, 0).tolist()} ms")
    ts_all = np.array(sorted(r["t0"] for r in ver))
    print(f"beacons per second in this {total:.0f} s, {fs/1e6:.0f} MHz window ({len(set(r['ch'] for r in ver))} of 69 channels in view): {len(ver)/total:.2f}")
    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        r = max(ver, key=lambda r: r["dev"])
        x = load(a.iq, fs, r["t0"] - 2e-3, r["dur"] + 4e-3)
        d, on = demod(x, fs, r["f0"] - fc)
        Tb = r["T"]; s0 = r["sfd_sample"]
        a0, a1 = int(s0 - 1.0e-3 * fs), int(s0 + 3.0e-3 * fs)
        t = (np.arange(a0, a1) - s0) / fs * 1e3
        ax[0].plot(t, d[a0:a1] / 1e3, lw=0.6)
        for kb in range(-50, 150):
            ax[0].axvline(kb * Tb / fs * 1e3, color="0.88", lw=0.5, zorder=0)
        ax[0].axvspan(0, 16 * Tb / fs * 1e3, color="C1", alpha=0.25, label="sync word 0x904E")
        ax[0].set_xlim(t[0], t[-1])
        ax[0].set_ylim(-80, 80)
        ax[0].set_xlabel("ms from the sync word")
        ax[0].set_ylabel("instantaneous frequency, kHz")
        ax[0].set_title(f"end of {r['preamble_bits']*Tb/fs*1e3:.0f} ms of 0x55, sync, frame — 20 µs bit grid")
        ax[0].legend(loc="upper right", fontsize=8)
        ax[1].hist(offs / 1e6, bins=np.arange(CH0 / 1e6 - 0.1, (CH0 + 69 * CH_STEP) / 1e6, 0.2), color="C0")
        ax[1].set_xlim((fc - fs / 2) / 1e6, (fc + fs / 2) / 1e6)
        ax[1].set_xlabel("MHz (bins = 902.2 + 0.2 n)")
        ax[1].set_ylabel("verified beacons")
        ax[1].set_title("every beacon on the 200 kHz grid")
        cols = {}
        for r2 in ver:
            cols.setdefault(r2["frame"]["gid"], f"C{len(cols) % 10}")
            ax[2].plot(r2["t0"], r2["ch"], "o", ms=4, color=cols[r2["frame"]["gid"]])
        ax[2].set_xlabel("time, s")
        ax[2].set_ylabel("channel n")
        ax[2].set_title(f"{len(ver)} beacons from {len(byid)} gateways (one colour each)")
        fig.tight_layout()
        Path(a.png).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(a.png, dpi=110)
        print(f"figure: {a.png}")


if __name__ == "__main__":
    main()
