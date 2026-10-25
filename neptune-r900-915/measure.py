#!/usr/bin/env python3
"""Measure the Neptune R900 water-meter grid off a 912.6 MHz capture.

Neptune's R900 radio modules sit on water meters and share the 910-920 MHz
ISM band with Itron ERT gas/electric modules (see ../ert-amr-915).  Per the
open decoders (rtl_433 neptune_r900.c, rtlamr):

  modulation   on-off keying, ~30 us chips  (rtl_433: OOK_PULSE_PCM, 30 us)
  preamble     alternating chips, then the sync AB 52 CC D2 + 1100101 (as rtl_433 aligns it)
  line code    168 chips after the sync, read in groups of 4: each group is
               one of the six weight-2 patterns 0011 0101 0110 1100 1010 1001
               = a base-6 digit; pairs of digits (base 36) -> 5 bits each
               -> 105 payload bits (ID, flags, consumption)
  integrity    NO CRC.  The only check is the code itself: every group of
               four chips must have exactly two ones.

This script derives, from raw IQ:
  1. which bursts are two-level (OOK) and which are not, by envelope
  2. the chip clock, by folding every level transition of every OOK burst
     modulo a trial period and taking the sharpest fold (no decoding)
  3. frames: preamble found, 42 groups all weight-2 ("code-valid"), ID masked
  4. how strong the integrity check is: the chance a random 168-chip string
     passes it, and how many bursts that pass are then consistent
  5. hop channels: carrier offset of every code-valid burst, comb search
  6. optional cross-check against rtl_433 on the same file

PRIVACY: meter IDs are masked to their last two digits and consumption
values are never printed.  Neighbours' water meters are not ours to publish.

Usage: python measure.py --iq r900_912p6_2400k.cs16 --fs 2400000 [--referee rtl433.jsonl] [--png figures/r900.png]
"""
import argparse
import collections
import json
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

# rtl_433 searches a zero-based bit stream for 55 55 55 A9 66 69 65; in chip terms (drop its leading 0) that is
# ...AA AB 52 CC D2 followed by 1100101: the sync runs 7 chips past 0xD2, and the data start there.
PRE_CHIPS = [int(c) for c in "".join(f"{v:08b}" for v in (0x55, 0x55, 0x55, 0xA9, 0x66, 0x69, 0x65))[1:]]
SIX = {0b0011: 0, 0b0101: 1, 0b0110: 2, 0b1100: 3, 0b1010: 4, 0b1001: 5}
T_CHIP_NOM = 30e-6
N_DATA = 168


# ---------------------------------------------------------------- code
def chips_to_payload(chips):
    """168 chips -> 105 bits via the 2-of-4 code, or None if any group is invalid."""
    digits = []
    for i in range(0, N_DATA, 4):
        g = (chips[i] << 3) | (chips[i + 1] << 2) | (chips[i + 2] << 1) | chips[i + 3]
        if g not in SIX:
            return None
        digits.append(SIX[g])
    bits = []
    for i in range(0, 42, 2):
        v = 6 * digits[i] + digits[i + 1]          # base 36, five bits
        bits += [(v >> k) & 1 for k in (4, 3, 2, 1, 0)]
    return bits


def payload_to_chips(bits):
    """Inverse of chips_to_payload for the synthetic self-test."""
    inv = {v: k for k, v in SIX.items()}
    chips = []
    for i in range(0, 105, 5):
        v = int("".join(map(str, bits[i:i + 5])), 2)
        for dgt in (v // 6, v % 6):
            g = inv[dgt]
            chips += [(g >> 3) & 1, (g >> 2) & 1, (g >> 1) & 1, g & 1]
    return chips


def meter_id(bits):
    """32-bit ID exactly as rtl_433 builds it: the first four bytes of the 105-bit payload, big-endian."""
    return int("".join(map(str, bits[:32])), 2)


def ook_iq(chips, fs, f_off, rng, snr_db=20, lead=40):
    chips = [0] * lead + list(chips) + [0] * lead
    n = int(len(chips) * T_CHIP_NOM * fs)
    idx = np.minimum((np.arange(n) / (T_CHIP_NOM * fs)).astype(int), len(chips) - 1)
    x = np.asarray(chips, float)[idx] * np.exp(2j * np.pi * f_off / fs * np.arange(n))
    noise = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    return x + noise * 10 ** (-snr_db / 20)


# ---------------------------------------------------------------- receiver
def burst_scan(path, fs, blk_s=5.0, win_s=100e-6):
    m = int(win_s * fs)
    env = []
    with open(path, "rb") as fh:
        while True:
            raw = np.fromfile(fh, dtype=np.int16, count=int(blk_s * fs) * 2)
            if raw.size < 2 * m:
                break
            p = raw[0::2].astype(np.float32) ** 2 + raw[1::2].astype(np.float32) ** 2
            env.append(np.sqrt(p[:p.size // m * m].reshape(-1, m).mean(1)))
    e = np.concatenate(env)
    med = np.median(e)
    thr = med + 6 * 1.4826 * np.median(np.abs(e - med))
    idx = np.flatnonzero(e > thr)
    segs = np.split(idx, np.flatnonzero(np.diff(idx) > 3) + 1) if idx.size else []
    return [(s[0] * win_s, (s[-1] - s[0] + 1) * win_s) for s in segs if (s[-1] - s[0] + 1) * win_s >= 2.5e-3], e.size * win_s


def load(path, fs, t0, dur):
    raw = np.fromfile(path, dtype=np.int16, offset=max(int(t0 * fs), 0) * 4, count=int(dur * fs) * 2)
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0


def carrier(x, fs):
    n = 1 << 15
    F = np.abs(np.fft.fft(x[:n] * np.hanning(min(n, x.size)), n))
    k = int(np.argmax(F))
    d = 0.5 * (F[k - 1] - F[k + 1]) / (F[k - 1] - 2 * F[k] + F[k + 1] + 1e-12)
    return ((k + d) / n * fs + fs / 2) % fs - fs / 2


def envelope(x, fs, f_off, taps):
    z = x * np.exp(-2j * np.pi * f_off / fs * np.arange(x.size))
    return np.abs(lfilter(taps, 1.0, z))


def two_level(env):
    hi = np.median(env[env > np.percentile(env, 75)])
    lo = np.median(env[env < np.percentile(env, 25)])
    return hi / max(lo, 1e-9), hi, lo


def transitions(env, hi, lo):
    lv = (env > 0.5 * (hi + lo)).astype(np.int8)
    on = np.flatnonzero(lv)
    if on.size < 8:
        return None, lv
    tr = np.flatnonzero(np.diff(lv[on[0]:on[-1] + 1]) != 0).astype(float) + on[0] + 0.5
    return tr, lv


def fold_sharpness(tr_s, period):
    """Circular coherence of transition times folded modulo one period."""
    return float(np.abs(np.mean(np.exp(2j * np.pi * tr_s / period))))


def slice_chips(lv, tr, sps):
    ph = np.angle(np.mean(np.exp(2j * np.pi * (tr % sps) / sps)))
    start = (ph / (2 * np.pi) * sps + sps / 2) % sps
    idx = np.arange(start, lv.size, sps)
    return lv[idx.astype(int)]


def find_frames(chips):
    s = "".join(map(str, chips))
    pat = "".join(map(str, PRE_CHIPS))
    out = []
    i = s.find(pat)
    while i >= 0:
        j = i + len(pat)
        if j + N_DATA <= len(s):
            bits = chips_to_payload([int(c) for c in s[j:j + N_DATA]])
            if bits is not None:
                out.append((j, bits))
        i = s.find(pat, i + 1)
    return out


def mask(i):
    return f"..{i % 100:02d}"


# ================================================================= main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--referee")
    ap.add_argument("--png")
    a = ap.parse_args()
    fs = a.fs
    dec_taps = firwin(129, 150e3, fs=fs)
    rng = np.random.default_rng(11)

    # ---- how strong is "all 42 groups weight-2" as an integrity check?
    p_group = 6 / 16
    print(f"integrity: a random 4-chip group is code-valid with p = 6/16; 42 groups -> p = {p_group**42:.1e} "
          f"(a random 168-chip string passes once in {1/p_group**42:.1e})")

    # ---- synthetic self-test
    ok = 0
    for k in range(4):
        bits = list(rng.integers(0, 2, 105))
        chips = PRE_CHIPS + payload_to_chips(bits)
        x = ook_iq(chips, fs, rng.uniform(-800e3, 800e3), rng)
        f = carrier(x, fs)
        env = envelope(x, fs, f, dec_taps)
        ratio, hi, lo = two_level(env)
        tr, lv = transitions(env, hi, lo)
        c2 = slice_chips(lv, tr, T_CHIP_NOM * fs)
        fr = find_frames(c2)
        ok += bool(fr and fr[0][1] == bits)
    print(f"synthetic self-test: {ok}/4 frames recovered bit-exact (random carrier offsets up to ±800 kHz, 20 dB SNR)")

    bursts, total_s = burst_scan(a.iq, fs)
    print(f"\nbursts >= 2.5 ms: {len(bursts)} in {total_s:.0f} s")
    rows = []
    for t0, dur in bursts:
        x = load(a.iq, fs, t0 - 0.3e-3, dur + 0.6e-3)
        f = carrier(x, fs)
        env = envelope(x, fs, f, dec_taps)
        ratio, hi, lo = two_level(env)
        tr, lv = transitions(env, hi, lo)
        rows.append(dict(t0=t0, dur=dur, f=f, ratio=ratio, tr=tr, lv=lv))
    ook = [r for r in rows if r["ratio"] >= 4 and r["tr"] is not None and r["tr"].size >= 20]
    print(f"two-level (OOK) envelopes: {len(ook)}   other (FSK or too weak): {len(rows) - len(ook)}")

    # ---- chip clock by folding, pooled over all OOK bursts, decoding-free
    periods = np.arange(26e-6, 36e-6, 0.01e-6)
    best = None
    for r in ook:
        tr_s = r["tr"] / fs
        sharp = np.array([fold_sharpness(tr_s - tr_s[0], p) for p in periods])
        r["T_fold"] = periods[int(np.argmax(sharp))]
    Ts = np.array([r["T_fold"] for r in ook])
    # keep bursts whose own fold lands near 30 us (the R900 population); ERT bursts fold at 30.5 us too, so verify below
    print(f"per-burst chip period from the fold: median {np.median(Ts)*1e6:.2f} us, "
          f"10-90 % {np.percentile(Ts,10)*1e6:.2f}-{np.percentile(Ts,90)*1e6:.2f} us")

    # ---- frames
    ver = []
    for r in ook:
        chips = slice_chips(r["lv"], r["tr"], r["T_fold"] * fs)
        fr = find_frames(chips)
        if fr:
            pos, bits = fr[0]
            r["bits"] = bits
            r["id"] = meter_id(bits)
            r["n_chips_after_sync"] = chips.size - pos
            ver.append(r)
    byid = collections.defaultdict(list)
    for r in ver:
        byid[r["id"]].append(r)
    print(f"code-valid R900 frames: {len(ver)}  from {len(byid)} meters (IDs masked)")
    if ver:
        Tv = np.array([r["T_fold"] for r in ver])
        print(f"chip period over verified frames: {np.median(Tv)*1e6:.3f} us -> {1/np.median(Tv):.0f} chips/s "
              f"(rtl_433's 30 us = 33333/s; a 32768 Hz watch crystal would be 30.52 us)  spread {np.std(Tv)*1e6:.2f} us")
        print(f"burst length of verified frames: {np.median([r['dur'] for r in ver])*1e3:.1f} ms  "
              f"(preamble + 63-chip sync + 168 data chips ~ 7.0 ms at 30.52 us)")
        print(f"envelope hi/lo of verified frames: median {np.median([r['ratio'] for r in ver]):.0f} "
              f"(OOK)  vs non-OOK bursts: {np.median([r['ratio'] for r in rows if r['ratio'] < 4]):.1f}")
        for mid, rs in sorted(byid.items(), key=lambda kv: -len(kv[1]))[:12]:
            ts = np.array(sorted(r["t0"] for r in rs))
            print(f"  meter {mask(mid)}: {len(rs)} bursts, carrier offsets {sorted(np.round([r['f']/1e3 for r in rs], 0))} kHz"
                  + (f", gaps {np.round(np.diff(ts), 1)} s" if ts.size > 1 else ""))
        # hop comb
        offs = np.array([r["f"] for r in ver])
        if offs.size >= 6:
            spac = np.arange(20e3, 600e3, 1e3)
            R = np.array([np.abs(np.mean(np.exp(2j * np.pi * offs / s))) for s in spac])
            kbest = int(np.argmax(R))
            print(f"hop comb over {offs.size} verified bursts: best spacing {spac[kbest]/1e3:.1f} kHz (R = {R[kbest]:.2f}); "
                  f"occupied slots {sorted(set(np.round(offs/spac[kbest]).astype(int)))} x {spac[kbest]/1e3:.0f} kHz")
    if a.referee and Path(a.referee).exists():
        cnt = collections.Counter()
        ids = set()
        for line in open(a.referee):
            try:
                j = json.loads(line)
            except ValueError:
                continue
            cnt[j.get("model")] += 1
            if j.get("model", "").startswith("Neptune"):
                ids.add(str(j.get("id"))[-2:])
        ours = {mask(m)[2:] for m in byid}
        print(f"referee (rtl_433 on the same file): {sum(cnt.values())} frames {dict(cnt)}")
        print(f"  R900 IDs (last 2 digits): ours {len(ours)}, rtl_433 {len(ids)}, common {len(ours & ids)}, "
              f"ours-only {len(ours - ids)}, rtl_433-only {len(ids - ours)}")

    if a.png and ver:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        r = max(ver, key=lambda r: r["ratio"])
        x = load(a.iq, fs, r["t0"] - 0.3e-3, r["dur"] + 0.6e-3)
        env = envelope(x, fs, r["f"], dec_taps)
        t = np.arange(env.size) / fs * 1e3
        ax[0].plot(t, env / env.max(), lw=0.7)
        for k in range(0, 60):
            ax[0].axvline(0.3 + k * r["T_fold"] * 1e3, color="0.85", lw=0.5, zorder=0)
        ax[0].set_xlim(0.3, 0.3 + 60 * r["T_fold"] * 1e3)
        ax[0].set_xlabel("ms")
        ax[0].set_ylabel("envelope")
        ax[0].set_title(f"one R900 burst on the {r['T_fold']*1e6:.1f} µs chip grid (on-off keyed)")
        tr_all = np.concatenate([(r2["tr"] - r2["tr"][0]) / fs for r2 in ver])
        sharp = [fold_sharpness(tr_all, p) for p in periods]
        ax[1].plot(periods * 1e6, sharp)
        ax[1].axvline(30.0, color="C3", ls="--", lw=0.8, label="30 µs (rtl_433)")
        ax[1].axvline(1e6 / 32768, color="C2", ls=":", lw=0.8, label="30.52 µs (32768 Hz)")
        ax[1].set_xlabel("trial chip period, µs")
        ax[1].set_ylabel("fold sharpness, all verified bursts")
        ax[1].set_title("the chip clock, from level changes alone")
        ax[1].legend(fontsize=8)
        cols = {}
        for r2 in ver:
            cols.setdefault(r2["id"], f"C{len(cols) % 10}")
            ax[2].plot(r2["t0"], r2["f"] / 1e3, "o", ms=4, color=cols[r2["id"]])
        ax[2].set_xlabel("time, s")
        ax[2].set_ylabel("carrier offset from 912.6 MHz, kHz")
        ax[2].set_title(f"{len(ver)} verified bursts, one colour per meter")
        fig.tight_layout()
        Path(a.png).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(a.png, dpi=110)
        print(f"\nfigure: {a.png}")


if __name__ == "__main__":
    main()
