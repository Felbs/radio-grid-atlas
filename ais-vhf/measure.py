#!/usr/bin/env python3
"""Measure the AIS grid off a 162 MHz marine-VHF capture.

AIS (ITU-R M.1371) is the ship-reporting grid: two 25 kHz channels
(A = 161.975 MHz, B = 162.025 MHz), 9600 bps GMSK (BT 0.4, h = 0.5),
NRZI + HDLC framing (0x7E flags, bit stuffing, CRC-16/X.25), and a
self-organized TDMA grid of 2250 slots per minute — one slot =
26.6667 ms = 256 bits — phase-locked to UTC by every transmitter.

This script derives, from raw IQ centered at 162.000 MHz:
  1. CRC-verified HDLC frames (the gold standard of proof)
  2. the 9600 Hz symbol clock (unrestricted level-change spectrum)
  3. GMSK deviation, with the measuring chain's bias calibrated
     out on synthetic bursts (synthetic TX first — house law)
  4. the SOTDMA slot comb: start-flag times of independent
     transmitters folded at 26.6667 ms
  5. if a sidecar JSON gives the capture UTC start: base-station
     msg-4 time-of-day vs the wall clock

PRIVACY: prints message types, source classes and masked MMSIs only.
Decoded positions are never printed.

Usage: python measure.py --iq ais_162.cs16 --fs 250000 [--meta ais_162.json]
"""
import argparse
import collections
import json
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter, filtfilt

HERE = Path(__file__).resolve().parent
RB = 9600.0                    # bits per second
SLOT = 60.0 / 2250             # 26.6667 ms
CHAN = {"A": -25e3, "B": +25e3}   # offsets from a 162.000 MHz dial
FLAG = [0, 1, 1, 1, 1, 1, 1, 0]
GOOD_RESIDUAL = 0xF0B8         # CRC-16/X.25 magic remainder

# ---------------------------------------------------------------- HDLC / CRC
def crc16_bits(bits, crc=0xFFFF):
    for b in bits:
        crc = (crc >> 1) ^ 0x8408 if (crc ^ int(b)) & 1 else crc >> 1
    return crc


def make_fcs(data_bits):
    crc = crc16_bits(data_bits) ^ 0xFFFF
    return [(crc >> i) & 1 for i in range(16)]          # LSB first on wire


def stuff(bits):
    out, run = [], 0
    for b in bits:
        out.append(b)
        run = run + 1 if b else 0
        if run == 5:
            out.append(0)
            run = 0
    return out


def unstuff(bits):
    out, run, i = [], 0, 0
    while i < len(bits):
        if run == 5:
            if bits[i]:
                return None                              # 6 ones inside frame
            run = 0
            i += 1
            continue
        out.append(bits[i])
        run = run + 1 if bits[i] else 0
        i += 1
    return out


def nrzi_encode(bits, level=0):
    out = []
    for b in bits:
        level ^= (b == 0)
        out.append(level)
    return out


def nrzi_decode(levels):
    l = np.asarray(levels, int)
    return (l[1:] == l[:-1]).astype(int)                 # change = 0


# ---------------------------------------------------------------- receiver
def fm_freq(z, fs):
    return np.angle(z[1:] * np.conj(z[:-1])) * fs / (2 * np.pi)


def demod(z, fs, smooth_taps):
    f = filtfilt(smooth_taps, [1.0], fm_freq(z, fs))
    cfo = np.median(f)
    return f - cfo, cfo


def recover_bits(fsm, fs):
    sps = fs / RB
    s = np.sign(fsm)
    zc = np.flatnonzero(np.diff(s) != 0).astype(float)
    if zc.size < 4:
        idx = np.arange(int(sps / 2), fsm.size, sps).astype(int)
        return (fsm[idx] > 0).astype(int), sps / 2
    ph = np.angle(np.mean(np.exp(2j * np.pi * (zc % sps) / sps)))
    start = (ph / (2 * np.pi) * sps + sps / 2) % sps
    idx = np.arange(start, fsm.size, sps).astype(int)
    return (fsm[idx] > 0).astype(int), start


def hunt_frames(levels):
    bits = nrzi_decode(levels)
    bstr = "".join(map(str, bits))
    pos, i = [], bstr.find("01111110")
    while i != -1:
        pos.append(i)
        i = bstr.find("01111110", i + 1)
    out = []
    for a in pos:
        for b in pos:
            if not 24 <= b - (a + 8) <= 1100:
                continue
            un = unstuff([int(c) for c in bstr[a + 8:b]])
            if un and len(un) >= 40 and crc16_bits(un) == GOOD_RESIDUAL:
                out.append({"start_flag": a, "end_flag": b,
                            "data_bits": un[:-16]})
    return out, pos, bits


def ais_fields(data_bits):
    """AIS fields sit on the per-octet bit-reversed stream."""
    bits = []
    for k in range(len(data_bits) // 8):
        bits.extend(data_bits[8 * k:8 * k + 8][::-1])

    def u(a, n):
        v = 0
        for i in range(n):
            v = (v << 1) | bits[a + i]
        return v
    return u(0, 6), u(8, 30), u


# ---------------------------------------------------------------- synthetic TX
def gmsk_burst(bits_lv, fs, fd=2400.0, bt=0.4):
    """GMSK by phase integration on the fs grid (fractional sps exact)."""
    sps = fs / RB
    n = int(round(len(bits_lv) * sps)) + 1
    f = np.zeros(n)
    for k, b in enumerate(bits_lv):
        f[int(round(k * sps)):int(round((k + 1) * sps))] = 2 * b - 1
    t = np.arange(-int(3 * sps), int(3 * sps) + 1) / sps
    sig = np.sqrt(np.log(2)) / (2 * np.pi * bt)
    g = np.exp(-t ** 2 / (2 * sig ** 2))
    fr = np.convolve(f, g / g.sum(), "same") * fd
    return np.exp(2j * np.pi * np.cumsum(fr) / fs)


def make_ais_bits(data_bits):
    body = stuff(list(data_bits) + make_fcs(data_bits))
    return nrzi_encode([0] * 8 + [0, 1] * 12 + FLAG + body + FLAG + [0] * 24)


# ---------------------------------------------------------------- deviation fit
def fit_burst(fsm, levels, ph, fr, fs, smooth_taps, bt=0.4, spsi=25):
    """LS-fit a regenerated freq trajectory to the measured one -> fd."""
    sps = fs / RB
    a_bit, b_bit = max(fr["start_flag"] - 24, 0), fr["end_flag"] + 8
    lv = levels[a_bit:b_bit + 1]
    nrz = np.repeat(2 * np.asarray(lv, float) - 1, spsi)
    t = np.arange(-3 * spsi, 3 * spsi + 1) / spsi
    sig = np.sqrt(np.log(2)) / (2 * np.pi * bt)
    g = np.exp(-t ** 2 / (2 * sig ** 2))
    mf = filtfilt(smooth_taps, [1.0], np.convolve(nrz, g / g.sum(), "same"))
    meas = fsm[int(ph + a_bit * sps):int(ph + (b_bit + 1) * sps)]
    m = np.interp(np.linspace(0, 1, meas.size),
                  np.linspace(0, 1, mf.size), mf)
    best = None
    for sh in range(-8, 9):
        mm = np.roll(m, sh)
        fd = float(mm @ meas) / float(mm @ mm)
        err = float(((meas - fd * mm) ** 2).sum())
        if best is None or err < best[0]:
            best = (err, fd)
    return best[1]


# ================================================================= main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--meta", default=None,
                    help="sidecar JSON with t_start_utc (enables msg-4 check)")
    a = ap.parse_args()
    fs = a.fs
    dec = max(1, int(fs / 50e3))
    fsd = fs / dec
    chan_taps = firwin(129, 12.5e3, fs=fs)
    smooth = firwin(31, 0.75 * RB, fs=fsd)
    rng = np.random.default_rng(3)

    # ---- SYNTHETIC TX FIRST: prove the whole receiver before the capture
    ok = 0
    for _ in range(5):
        data = list(rng.integers(0, 2, 168))
        x = gmsk_burst(make_ais_bits(data), fs)
        nn = np.arange(x.size)
        y = x * np.exp(2j * np.pi * (-25e3 + rng.uniform(-500, 500)) / fs * nn)
        y = np.concatenate([np.zeros(3000), y, np.zeros(3000)])
        y = y + 0.05 * (rng.standard_normal(y.size)
                        + 1j * rng.standard_normal(y.size))
        z = y * np.exp(-2j * np.pi * -25e3 / fs * np.arange(y.size))
        z = lfilter(chan_taps, 1.0, z)[::dec]
        fsm, _ = demod(z, fsd, smooth)
        lv, ph = recover_bits(fsm, fsd)
        frames, _, _ = hunt_frames(lv)
        ok += any(f["data_bits"] == data for f in frames)
    print(f"synthetic self-test: {ok}/5 CRC-verified recoveries")
    assert ok == 5, "receiver failed its own synthetic gauntlet"

    # ---- deviation-ruler calibration: known fd=2400 through the same chain
    fits = []
    for _ in range(6):
        data = list(rng.integers(0, 2, 168))
        x = gmsk_burst(make_ais_bits(data), fs, fd=2400.0, bt=0.4)
        nn = np.arange(x.size)
        y = np.concatenate([np.zeros(3000),
                            x * np.exp(2j * np.pi * -25.15e3 / fs * nn),
                            np.zeros(3000)])
        y = y + 0.03 * (rng.standard_normal(y.size)
                        + 1j * rng.standard_normal(y.size))
        z = y * np.exp(-2j * np.pi * -25e3 / fs * np.arange(y.size))
        z = lfilter(chan_taps, 1.0, z)[::dec]
        fsm, _ = demod(z, fsd, smooth)
        lv, ph = recover_bits(fsm, fsd)
        frames, _, _ = hunt_frames(lv)
        if frames:
            fits.append(fit_burst(fsm, lv, ph, frames[0], fsd, smooth))
    bias = np.median(fits) / 2400.0
    print(f"deviation ruler: chain reads x{bias:.4f} on a known 2400 Hz burst")

    # ---- stream both channels: 1 ms envelope -> burst candidates
    nsamp = Path(a.iq).stat().st_size // 4
    env = {c: [] for c in CHAN}
    blk = int(30 * fs)
    with open(a.iq, "rb") as fh:
        off = 0
        while True:
            raw = np.fromfile(fh, dtype=np.int16, count=blk * 2)
            if raw.size == 0:
                break
            x = (raw[0::2].astype(np.float32)
                 + 1j * raw[1::2].astype(np.float32)) / 32768.0
            nn = np.arange(off, off + x.size)
            off += x.size
            for c, fo in CHAN.items():
                z = x * np.exp(-2j * np.pi * fo / fs * nn).astype(np.complex64)
                z = lfilter(chan_taps, 1.0, z)[::dec]
                m = int(1e-3 * fsd)
                e = np.abs(z[:z.size // m * m]).reshape(-1, m).mean(1)
                env[c].append(e.astype(np.float32))
    bursts = []
    for c in CHAN:
        e = np.concatenate(env[c])
        med = np.median(e)
        thr = med + 8 * 1.4826 * np.median(np.abs(e - med))
        idx = np.flatnonzero(e > thr)
        segs = (np.split(idx, np.flatnonzero(np.diff(idx) > 5) + 1)
                if idx.size else [])
        bursts += [(c, s[0] * 1e-3, (s[-1] - s[0] + 1) * 1e-3)
                   for s in segs if s[-1] - s[0] >= 3]
    print(f"burst candidates: "
          f"{collections.Counter(c for c, *_ in bursts)} in {nsamp/fs:.0f} s")

    # ---- decode every candidate
    def load(t0, dur):
        raw = np.fromfile(a.iq, dtype=np.int16, offset=int(t0 * fs) * 4,
                          count=int(dur * fs) * 2)
        return (raw[0::2].astype(np.float32)
                + 1j * raw[1::2].astype(np.float32)) / 32768.0

    NF = 1 << 16
    clk = np.zeros(NF // 2 + 1)
    devs = []
    recs = []
    fits_real = []
    for c, t0, dur in bursts:
        ta = max(t0 - 0.010, 0)
        x = load(ta, dur + 0.035)
        nn = np.arange(int(ta * fs), int(ta * fs) + x.size)
        z = x * np.exp(-2j * np.pi * CHAN[c] / fs * nn).astype(np.complex64)
        z = lfilter(chan_taps, 1.0, z)[::dec]
        fsm, cfo = demod(z, fsd, smooth)
        lv, ph = recover_bits(fsm, fsd)
        frames, flags, bits = hunt_frames(lv)
        if not frames:
            continue
        fr = frames[0]
        sps = fsd / RB
        t_flag = ta + (ph + (fr["start_flag"] + 1) * sps) / fsd
        mt, mmsi, u = ais_fields(fr["data_bits"])
        # training check: the 20 bits before the flag must alternate 0101
        # (the first ~2 of the 24 training bits ride the TX power ramp)
        pre = list(bits[max(fr["start_flag"] - 20, 0):fr["start_flag"]])
        train_ok = (len(pre) == 20
                    and all(pre[i] != pre[i + 1] for i in range(19)))
        recs.append({"ch": c, "t": t_flag, "mt": mt, "mmsi": mmsi,
                     "nbits": len(fr["data_bits"]), "cfo": cfo, "u": u,
                     "train": train_ok, "fsm": None})
        # symbol-clock impulses
        i0 = int(ph + max(fr["start_flag"] - 24, 0) * sps)
        i1 = int(ph + (fr["end_flag"] + 8) * sps)
        seg = fsm[i0:i1]
        imp = np.zeros(seg.size)
        imp[np.flatnonzero(np.diff(np.sign(seg)) != 0)] = 1.0
        clk += np.abs(np.fft.rfft(imp, NF)) ** 2
        devs.append(seg)
        fits_real.append((dur, fit_burst(fsm, lv, ph, fr, fsd, smooth)))
        if len(recs) == 1:      # keep one strong burst for the figure
            fig_burst = (fsm, ph, fr, sps)

    print(f"CRC-verified frames: {len(recs)} "
          f"(A {sum(r['ch']=='A' for r in recs)}, "
          f"B {sum(r['ch']=='B' for r in recs)})")
    print(f"training 0101 (last 20 bits before flag) bit-exact: "
          f"{sum(r['train'] for r in recs)}/{len(recs)}")
    print(f"payload sizes: "
          f"{dict(collections.Counter(r['nbits'] for r in recs))}")
    print(f"channel CFO (our TCXO + theirs): "
          f"A {np.median([r['cfo'] for r in recs if r['ch']=='A']):+.0f} Hz, "
          f"B {np.median([r['cfo'] for r in recs if r['ch']=='B']):+.0f} Hz")

    # ---- census (masked)
    def klass(m):
        s = f"{m:09d}"
        return ("base station" if s.startswith("00")
                else "AtoN" if s.startswith("99")
                else "Class B vessel" if s.startswith("33")
                else "vessel")
    per = collections.defaultdict(list)
    for r in recs:
        per[(r["mmsi"], r["mt"])].append(r["t"])
    print("sources (MMSIs masked):")
    for (m, mt), ts in sorted(per.items()):
        ts = np.sort(ts)
        iv = (f", median interval {np.median(np.diff(ts)):.3f} s"
              if ts.size > 2 else "")
        print(f"  ...{f'{m:09d}'[-3:]} ({klass(m)}) msg type {mt}: "
              f"n={len(ts)}{iv}")

    # ---- symbol clock (unrestricted search above 1 kHz)
    fx = np.fft.rfftfreq(NF, 1 / fsd)
    band = fx > 1000
    i = np.argmax(clk * band)
    d = 0.5 * (clk[i - 1] - clk[i + 1]) / (clk[i - 1] - 2 * clk[i] + clk[i + 1])
    f_sym = fx[i] + d * (fx[1] - fx[0])
    snr = 10 * np.log10(clk[i] / np.median(clk[band]))
    print(f"symbol clock: {f_sym:.2f} Hz ({snr:+.1f} dB) — published 9600")

    # ---- deviation (calibrated)
    fd_raw = np.median([f for _, f in fits_real])
    print(f"GMSK deviation: raw fit {fd_raw:+.0f} Hz, "
          f"calibrated {fd_raw / bias:+.0f} Hz "
          f"-> h = {2 * fd_raw / bias / RB:.3f} (published 2400 Hz, h = 0.5)")

    # ---- the slot comb
    groups = collections.defaultdict(list)
    for r in recs:
        groups[(r["mmsi"], r["ch"])].append(r["t"])
    print(f"SOTDMA slot comb (T = {SLOT*1e3:.4f} ms), start-flag time mod T:")
    phases = {}
    for k, ts in sorted(groups.items()):
        ts = np.array(ts)
        v = np.mean(np.exp(2j * np.pi * (ts % SLOT) / SLOT))
        phases[k] = (np.angle(v) / (2 * np.pi)) % 1.0
        print(f"  ...{f'{k[0]:09d}'[-3:]} ch {k[1]}: n={len(ts)} "
              f"R={np.abs(v):.4f} phase={phases[k]:.4f}")
    ks = list(phases)
    dmax = max((min((phases[p] - phases[q]) % 1, (phases[q] - phases[p]) % 1)
                for i, p in enumerate(ks) for q in ks[i + 1:]), default=0)
    resid = np.concatenate(
        [(((np.array(ts) % SLOT) / SLOT - phases[k] + .5) % 1) - .5
         for k, ts in groups.items()])
    print(f"  max pairwise phase spread: {dmax * SLOT * 1e3:.2f} ms; "
          f"residual jitter {resid.std() * SLOT * 1e6:.0f} us rms")

    # ---- msg 4 vs wall clock
    if a.meta:
        from datetime import datetime, timedelta, timezone
        t0u = datetime.fromisoformat(json.load(open(a.meta))["t_start_utc"])
        offs, dates = [], collections.Counter()
        for r in recs:
            if r["mt"] != 4:
                continue
            u = r["u"]
            dates[f"{u(38,14):04d}-{u(52,4):02d}-{u(56,5):02d}"] += 1
            tf = t0u + timedelta(seconds=r["t"])
            fsec = tf.hour * 3600 + tf.minute * 60 + tf.second \
                + tf.microsecond / 1e6
            offs.append((u(61, 5) * 3600 + u(66, 6) * 60 + u(72, 6)) - fsec)
        if offs:
            offs = np.array(offs)
            print(f"msg-4 UTC time-of-day vs capture clock: n={offs.size}, "
                  f"offset {offs.mean():+.3f} s, spread {offs.std()*1e3:.1f} ms")
            print(f"msg-4 date fields: {dict(dates)} "
                  f"(capture date {t0u.date()})")

    # ---- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = HERE / "figures"
    out.mkdir(exist_ok=True)

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13.5, 4.0))
    fsm, ph, fr, sps = fig_burst
    b0 = max(fr["start_flag"] - 26, 0)
    i0, i1 = int(ph + b0 * sps), int(ph + (fr["start_flag"] + 40) * sps)
    tt = (np.arange(i0, i1) - ph) / fsd * 1e3
    a1.plot(tt, fsm[i0:i1] / 1e3, lw=.8, color="#0aa")
    for lab, x0, x1, col in [
            ("training 0101", fr["start_flag"] - 24, fr["start_flag"], "#dd3"),
            ("flag 0x7E", fr["start_flag"], fr["start_flag"] + 8, "#d33")]:
        a1.axvspan(x0 / RB * 1e3, x1 / RB * 1e3, color=col, alpha=.18,
                   label=lab)
    a1.axhline(2.4, color="#d33", ls=":", lw=1)
    a1.axhline(-2.4, color="#d33", ls=":", lw=1)
    a1.legend(loc="lower right", fontsize=8)
    a1.set(title="one burst, demodulated frequency\n(dotted: published "
                 "+-2400 Hz)", xlabel="ms", ylabel="kHz from carrier")

    mm = (fx > 9400) & (fx < 9800)
    a2.semilogy(fx[mm], clk[mm] + 1e-9, lw=.7, color="#0aa")
    a2.axvline(9600, color="#d33", ls=":", lw=1)
    a2.set(title=f"level-change spectrum, all bursts\nsymbol clock "
                 f"{f_sym:.2f} Hz (dotted: 9600)", xlabel="Hz")

    dv = np.concatenate(devs)
    h, edges = np.histogram(dv, bins=300, range=(-4000, 4000))
    a3.fill_between((edges[:-1] + edges[1:]) / 2 / 1e3, h, color="#0aa",
                    step="mid")
    for v in (-2.4, 2.4):
        a3.axvline(v, color="#d33", ls=":", lw=1)
    a3.set(title="instantaneous frequency, all bursts\n(dotted: +-2400 Hz)",
           xlabel="kHz from carrier", ylabel="samples")
    fig.tight_layout()
    fig.savefig(out / "ais.png", dpi=110)
    print(f"figure -> {out / 'ais.png'}")

    fig2, (b1, b2) = plt.subplots(2, 1, figsize=(11.5, 5.6),
                                  height_ratios=[1, 1.4], sharex=True)
    styles = {}
    pal = ["#0aa", "#d33", "#a70", "#55c", "#7a0"]
    for r in recs:
        k = (r["mmsi"], r["ch"])
        if k not in styles:
            styles[k] = pal[len(styles) % len(pal)]
    for k, ts in groups.items():
        lab = f"...{f'{k[0]:09d}'[-3:]} {klass(k[0]).split()[0]} ch {k[1]}"
        b1.plot(ts, [1 if k[1] == "A" else 0] * len(ts), "|", ms=14,
                color=styles[k], label=lab)
        b2.plot(ts, (np.array(ts) % SLOT) * 1e3, ".", ms=5,
                color=styles[k])
    b1.set(yticks=[0, 1], yticklabels=["ch B", "ch A"],
           title="every CRC-verified burst in 300 s")
    b1.legend(loc="center right", fontsize=7, ncol=1)
    b2.axhline(np.mean([(t % SLOT) * 1e3 for ts in groups.values()
                        for t in ts]), color="#d33", lw=.6)
    b2.set(ylim=(0, SLOT * 1e3), xlabel="capture time (s)",
           ylabel="start-flag time mod 26.667 ms (ms)",
           title="three independent transmitters, one slot comb")
    fig2.tight_layout()
    fig2.savefig(out / "ais_slots.png", dpi=110)
    print(f"figure -> {out / 'ais_slots.png'}")


if __name__ == "__main__":
    main()
