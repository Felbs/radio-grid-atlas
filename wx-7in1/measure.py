#!/usr/bin/env python3
"""Measure the grid of two consumer 7-in-1 weather stations from raw IQ.

Both are 2-FSK bursts from a solar/battery outdoor unit to an indoor
display, on fixed channels in different ISM bands:

  station A  433.92 MHz band   Emax W6 family (rtl_433 emax.c; sold under many
             names, this one answers as "IMETEO-X6"): preamble 0xAA..,
             sync 0xCACA54, 32-byte message, 8-bit additive checksum over
             the first 31 bytes; a new random 12-bit ID at every power-up
  station B  915.00 MHz band   Fine Offset WH65B family: preamble 0xAA..,
             sync 0x2DD4, 17-byte frame, CRC-8 (poly 0x31) + 8-bit sum

This script derives, from the IQ alone:
  1. the bit period of each station, from the run lengths of the FSK
     discriminator (no decoding needed), pooled and compared with nominal
  2. the deviation (tone spacing) and each burst's carrier offset
  3. verified frames: sync found, checksum/CRC pass, both polarities tried
  4. the report cadence: the interval between verified bursts of a station
  5. optional cross-check against rtl_433 run on the same file

PRIVACY: station IDs are masked to their last two hex digits and no
readings are printed. One of these stations is a neighbour's.

Usage:
  python measure.py --iq433 wx433.cs16 --iq915 wx915.cs16 --fs 250000 \
      [--ref433 rtl433_433.jsonl] [--ref915 rtl433_915.jsonl] [--png figures/wx.png]
"""
import argparse
import collections
import json
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

HERE = Path(__file__).resolve().parent

STATIONS = {
    "A-433": dict(sync=bytes([0xCA, 0xCA, 0x54]), nbytes=32, t_bit_nom=90e-6, dev_nom=37e3, interval_nom=30.0),
    "B-915": dict(sync=bytes([0x2D, 0xD4]), nbytes=17, t_bit_nom=58e-6, dev_nom=37e3, interval_nom=16.0),
}


# ---------------------------------------------------------------- integrity
def crc8(data, poly=0x31, init=0):
    c = init
    for b in data:
        c ^= b
        for _ in range(8):
            c = ((c << 1) ^ poly) & 0xFF if c & 0x80 else (c << 1) & 0xFF
    return c


def verify(kind, b):
    """Frame bytes (sync removed) -> dict or None."""
    if kind == "A-433":
        if len(b) < 32 or (sum(b[:31]) & 0xFF) != b[31]:
            return None
        return {"id": (b[2] << 4) | (b[3] >> 4), "ch": b[1] & 0x0F, "ctr": None}
    if kind == "B-915":
        if len(b) < 17 or b[0] != 0x24 or crc8(b[:15]) != b[15] or (sum(b[:16]) & 0xFF) != b[16]:
            return None
        return {"id": b[1], "ch": 0, "ctr": None}
    return None


# ---------------------------------------------------------------- synthetic TX
def synth_frame(kind, rng):
    if kind == "A-433":
        b = bytearray(rng.integers(0, 256, 32, dtype=np.uint8).tolist())
        b[0], b[1] = 0xAA, 0x04
        b[31] = sum(b[:31]) & 0xFF
    else:
        b = bytearray(rng.integers(0, 256, 17, dtype=np.uint8).tolist())
        b[0] = 0x24
        b[15] = crc8(b[:15])
        b[16] = sum(b[:16]) & 0xFF
    return bytes(b)


def fsk_iq(bits, fs, t_bit, dev, f_off, snr_db, rng):
    """NRZ bits -> 2-FSK IQ (1 = upper tone), with a carrier offset and AWGN."""
    n = int(len(bits) * t_bit * fs)
    idx = np.minimum((np.arange(n) / (t_bit * fs)).astype(int), len(bits) - 1)
    f = np.where(np.asarray(bits)[idx] == 1, dev, -dev) + f_off
    ph = 2 * np.pi * np.cumsum(f) / fs
    x = np.exp(1j * ph)
    pad = np.zeros(int(2e-3 * fs), complex)
    x = np.concatenate([pad, x, pad])
    noise = (rng.standard_normal(x.size) + 1j * rng.standard_normal(x.size)) / np.sqrt(2)
    return x + noise * 10 ** (-snr_db / 20)


# ---------------------------------------------------------------- receiver
def burst_scan(path, fs, win_s=100e-6, min_s=3e-3, thr_mult=6.0, blk=5_000_000):
    """Wideband envelope -> (start, duration) list in seconds; also returns the file length."""
    m = int(win_s * fs)
    mm = np.memmap(path, np.int16, mode="r")
    env = []
    for a in range(0, mm.size, 2 * blk):
        raw = np.asarray(mm[a:a + 2 * blk], dtype=np.float32)
        p = raw[0::2] ** 2 + raw[1::2] ** 2
        env.append(np.sqrt(p[:p.size // m * m].reshape(-1, m).mean(1)))
    e = np.concatenate(env)
    thr = np.median(e) * thr_mult
    idx = np.flatnonzero(e > thr)
    segs = np.split(idx, np.flatnonzero(np.diff(idx) > 5) + 1) if idx.size else []
    out = [(s[0] * win_s, (s[-1] - s[0] + 1) * win_s) for s in segs if (s[-1] - s[0] + 1) * win_s >= min_s]
    return out, mm.size // 2 / fs


def load(path, fs, t0, dur):
    raw = np.fromfile(path, dtype=np.int16, offset=max(int(t0 * fs), 0) * 4, count=int(dur * fs) * 2)
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0


def centroid(x, fs):
    """Power-weighted centre of the strongest 200 kHz region (FSK has no carrier line)."""
    n = min(1 << 14, x.size)
    F = np.abs(np.fft.fftshift(np.fft.fft(x[:n] * np.hanning(n), n))) ** 2
    fax = np.fft.fftshift(np.fft.fftfreq(n, 1 / fs))
    k = int(np.argmax(F))
    w = np.abs(fax - fax[k]) < 100e3
    return float(np.sum(fax[w] * F[w]) / np.sum(F[w]))


def demod(x, fs, f_off, bw=90e3):
    z = x * np.exp(-2j * np.pi * f_off / fs * np.arange(x.size))
    z = lfilter(firwin(101, bw / (fs / 2)), 1.0, z)
    d = np.angle(z[1:] * np.conj(z[:-1])) * fs / (2 * np.pi)
    on = np.abs(z[1:]) > 0.5 * np.percentile(np.abs(z[1:]), 90)
    return d, on


def bit_period(d, on, fs, t_nom):
    """Bit period from level runs of the discriminator: every run is an integer
    number of bits, so T = sum(run) / sum(round(run / T_nom)) over clean runs."""
    dd = d[on]
    if dd.size < 50:
        return None, None, None
    c = 0.5 * (np.percentile(dd, 10) + np.percentile(dd, 90))
    lv = (d > c).astype(np.int8)
    inside = np.flatnonzero(on)
    tr = np.flatnonzero(np.diff(lv[inside[0]:inside[-1] + 1]) != 0)
    runs = np.diff(tr) / fs
    k = np.round(runs / t_nom)
    good = (k >= 1) & (k <= 8) & (np.abs(runs / np.maximum(k, 1) - t_nom) < 0.2 * t_nom)
    if good.sum() < 20:
        return None, None, None
    T = runs[good].sum() / k[good].sum()
    return T, c, runs


def slice_bits(d, on, fs, T, c):
    best = None
    sps = T * fs
    for ph in np.arange(0, sps, sps / 8):
        idx = (ph + sps / 2 + np.arange(int((d.size - ph - sps) / sps)) * sps).astype(int)
        q = np.abs(d[idx] - c).mean()
        if best is None or q > best[0]:
            best = (q, (d[idx] > c).astype(np.uint8))
    return best[1]


def frames_in(bits, sync, nbytes):
    s = "".join(map(str, bits))
    pat = "".join(f"{v:08b}" for v in sync)
    out = []
    i = s.find(pat)
    while i >= 0:
        j = i + len(pat)
        if j + nbytes * 8 <= len(s):
            out.append((i, bytes(int(s[k:k + 8], 2) for k in range(j, j + nbytes * 8, 8))))
        i = s.find(pat, i + 1)
    return out


def process(x, fs, kind):
    st = STATIONS[kind]
    f_off = centroid(x, fs)
    d, on = demod(x, fs, f_off)
    T, c, runs = bit_period(d, on, fs, st["t_bit_nom"])
    if T is None:
        return None
    dd = d[on]
    dev = 0.5 * (np.percentile(dd, 90) - np.percentile(dd, 10))
    hit = None
    for inv in (0, 1):
        bits = slice_bits(d, on, fs, T, c) ^ inv
        for pos, b in frames_in(bits, st["sync"], st["nbytes"]):
            v = verify(kind, b)
            if v:
                hit = dict(v, inv=inv, pos=pos, bits=bits)
                break
        if hit:
            break
    return dict(f_off=f_off, T=T, dev=dev, runs=runs, d=d, c=c, frame=hit)


def mask(i):
    return f"..{i % 100:02d}"


# ================================================================= main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq433")
    ap.add_argument("--iq915")
    ap.add_argument("--fs", type=float, default=250e3)
    ap.add_argument("--ref433")
    ap.add_argument("--ref915")
    ap.add_argument("--png")
    ap.add_argument("--lora-log", help="JSONL from the LoRa-card receiver (github.com/Felbs/helium-upcycle, docs/08 + weather/): one verified frame per line with a 't' timestamp; cadence of station B is measured from it")
    a = ap.parse_args()
    fs = a.fs
    rng = np.random.default_rng(3)

    # ---- synthetic self-test: the receiver must verify its own frames
    ok = 0
    for kind in STATIONS:
        st = STATIONS[kind]
        for k in range(3):
            fr = synth_frame(kind, rng)
            raw = bytes([0xAA] * 6) + st["sync"] + fr
            bits = [int(c) for c in "".join(f"{v:08b}" for v in raw)]
            x = fsk_iq(bits, fs, st["t_bit_nom"], st["dev_nom"], rng.uniform(-30e3, 30e3), 20, rng)
            r = process(x, fs, kind)
            ok += bool(r and r["frame"] and verify(kind, raw[6 + len(st['sync']):]) == {k2: v for k2, v in r["frame"].items() if k2 in ("id", "ch", "ctr")})
    print(f"synthetic self-test: {ok}/6 frames verified bit-exact (3 per station, random carrier offsets, 20 dB SNR)")

    fig_data = {}
    for kind, path, ref in (("A-433", a.iq433, a.ref433), ("B-915", a.iq915, a.ref915)):
        if not path:
            continue
        st = STATIONS[kind]
        bursts, total_s = burst_scan(path, fs, min_s=(50e-3 if kind == "A-433" else 4e-3), thr_mult=8.0)
        print(f"\n== {kind}: {path}  ({total_s:.0f} s)")
        print(f"bursts >= {'50' if kind == 'A-433' else '4'} ms at 8x the floor: {len(bursts)}")
        results = []
        for t0, dur in bursts:
            x = load(path, fs, t0 - 1e-3, dur + 2e-3)
            r = process(x, fs, kind)
            if r:
                r["t0"] = t0
                r["dur"] = dur
                results.append(r)
        ver = [r for r in results if r["frame"]]
        byid = collections.defaultdict(list)
        for r in ver:
            byid[(r["frame"]["id"], r["frame"]["ch"])].append(r)
        print(f"verified frames: {len(ver)} / {len(results)} FSK-shaped bursts   polarity: "
              f"{collections.Counter(r['frame']['inv'] for r in ver)}  (inv=1: a logical 1 is the LOWER tone)")
        print(f"stations heard (ID masked): {len(byid)}")
        T_all = np.array([r["T"] for r in ver])
        dev_all = np.array([r["dev"] for r in ver])
        if ver:
            print(f"bit period: {np.median(T_all)*1e6:.2f} us median over verified bursts "
                  f"(nominal {st['t_bit_nom']*1e6:.0f} us -> {1/np.median(T_all):.0f} bit/s, "
                  f"{(np.median(T_all)/st['t_bit_nom']-1)*1e6:+.0f} ppm); spread {np.std(T_all)*1e6:.2f} us")
            print(f"deviation: ±{np.median(dev_all)/1e3:.1f} kHz median (10-90 % tone levels)  "
                  f"burst length: {np.median([r['dur'] for r in ver])*1e3:.1f} ms")
        for (sid, ch), rs in sorted(byid.items(), key=lambda kv: -len(kv[1])):
            ts = np.array(sorted(r["t0"] for r in rs))
            iv = np.diff(ts)
            base = st["interval_nom"]
            k = np.round(iv / base)
            cad = iv[k >= 1] / k[k >= 1] if iv.size else np.array([])
            offs = np.array([r["f_off"] for r in rs])
            line = (f"  station {mask(sid)} ch{ch}: {len(rs)} verified bursts, carrier {np.median(offs)/1e3:+.1f} kHz "
                    f"from dial (spread {np.std(offs)/1e3:.1f} kHz)")
            if cad.size:
                line += (f"; cadence {np.median(cad):.3f} s (nominal {base:.0f} s, rms jitter {np.std(cad)*1e3:.0f} ms, "
                         f"{int(np.sum(k==1))} consecutive / {int(np.sum(k>1))} with gaps)")
            print(line)
        # run-length evidence for the bit clock, pooled, decoding-free
        runs = np.concatenate([r["runs"] for r in results if r["runs"] is not None]) if results else np.array([])
        if runs.size:
            h, edges = np.histogram(runs * 1e6, bins=np.arange(0, 8.5 * st["t_bit_nom"] * 1e6, 2))
            peaks = edges[np.argsort(h)[::-1][:4]]
            print(f"run-length peaks (us, all FSK bursts, no decoding): {sorted(np.round(peaks, 0))}")
        if ref and Path(ref).exists():
            cnt = collections.Counter()
            ids = set()
            for line in open(ref):
                try:
                    j = json.loads(line)
                except ValueError:
                    continue
                cnt[j.get("model")] += 1
                if j.get("id") is not None:
                    ids.add(str(j["id"])[-2:])
            print(f"referee (rtl_433 on the same file): {sum(cnt.values())} frames {dict(cnt)}; IDs (last 2 digits) {sorted(ids)}")
        fig_data[kind] = dict(results=results, ver=ver, st=st)

    lora = None
    if a.lora_log and Path(a.lora_log).exists():
        import datetime as dt
        ts = np.array(sorted(dt.datetime.fromisoformat(json.loads(l)["t"]).timestamp() for l in open(a.lora_log) if l.strip()))
        iv = np.diff(ts)
        k = np.round(iv / 16.0)
        per = iv[k >= 1] / k[k >= 1]
        hours = (ts[-1] - ts[0]) / 3600
        lora = dict(ts=ts, per=per, k=k)
        print(f"\n== B-915 via the LoRa card ({a.lora_log}): {ts.size} CRC+sum-verified frames in {hours:.2f} h "
              f"= {ts.size/hours:.0f}/h of 225 possible ({100*ts.size/hours/225:.0f} %)")
        print(f"cadence from receive timestamps: {np.median(per):.3f} s per frame (nominal 16 s), rms {np.std(per)*1e3:.0f} ms "
              f"(includes USB/host timing); consecutive {int(np.sum(k==1))}, one missed {int(np.sum(k==2))}, two missed {int(np.sum(k==3))}")

    if a.png and fig_data:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        # left: one 433 burst's discriminator on its bit grid with the sync shaded
        kind = "A-433" if fig_data.get("A-433", {}).get("ver") else next(k for k in fig_data if fig_data[k]["ver"])
        r = max(fig_data[kind]["ver"], key=lambda r: r["dev"])
        fr = r["frame"]
        sps = r["T"] * fs
        t = np.arange(r["d"].size) / fs * 1e3
        ax[0].plot(t, r["d"] / 1e3, lw=0.6, color="C0")
        n_sync = len(fig_data[kind]["st"]["sync"]) * 8
        x0 = (fr["pos"]) * r["T"] * 1e3 + 0.0
        for kb in range(0, 70):
            ax[0].axvline(x0 + kb * r["T"] * 1e3, color="0.85", lw=0.5, zorder=0)
        ax[0].axvspan(x0, x0 + n_sync * r["T"] * 1e3, color="C1", alpha=0.25, label="sync word")
        ax[0].set_xlim(x0 - 1.0, x0 + 70 * r["T"] * 1e3)
        ax[0].set_ylim(-80, 80)
        ax[0].set_xlabel("ms")
        ax[0].set_ylabel("instantaneous frequency, kHz")
        ax[0].set_title(f"{kind}: discriminator on the {r['T']*1e6:.1f} µs bit grid")
        ax[0].legend(loc="upper right", fontsize=8)
        # middle: run-length histograms, both stations
        for kind2, col in (("A-433", "C0"), ("B-915", "C3")):
            if kind2 not in fig_data or not fig_data[kind2]["ver"]:
                continue
            fd = fig_data[kind2]
            runs = np.concatenate([x["runs"] for x in fd["ver"]]) * 1e6
            ax[1].hist(runs, bins=np.arange(0, 600, 4), color=col, alpha=0.6, label=f"{kind2}")
            Tm = np.median([x["T"] for x in fd["ver"]]) * 1e6
            for k in range(1, 6):
                ax[1].axvline(k * Tm, color=col, lw=0.6, ls="--")
        ax[1].set_xlabel("level run length, µs")
        ax[1].set_ylabel("count (verified bursts)")
        ax[1].set_title("runs land on integer bits")
        ax[1].legend(fontsize=8)
        # right: cadence
        if lora is not None:
            t = lora["ts"] - lora["ts"][0]
            ax[2].plot(t, t % 16.0, ".", color="C3", ms=3, label="B-915 (LoRa-card timestamps)")
        for kind2, col in (("A-433", "C0"),):
            if kind2 not in fig_data or not fig_data[kind2]["ver"]:
                continue
            fd = fig_data[kind2]
            byid = collections.defaultdict(list)
            for x in fd["ver"]:
                byid[(x["frame"]["id"], x["frame"]["ch"])].append(x["t0"])
            for sid, ts in byid.items():
                ts = np.array(sorted(ts))
                if ts.size > 2:
                    ax[2].plot(ts - ts[0], (ts - ts[0]) % fd["st"]["interval_nom"], ".", color=col, ms=4, label="A-433 (SDR)")
            ax[2].axhline(0, color="0.7", lw=0.5)
        ax[2].set_xlabel("time since a station's first burst, s")
        ax[2].set_ylabel("time modulo nominal interval, s")
        ax[2].set_title("cadence: 30 s (433) and 16 s (915); flat = on the grid")
        ax[2].legend(fontsize=8)
        fig.tight_layout()
        Path(a.png).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(a.png, dpi=110)
        print(f"\nfigure: {a.png}")


if __name__ == "__main__":
    main()
