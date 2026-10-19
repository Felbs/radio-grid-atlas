#!/usr/bin/env python3
"""Measure the ERT utility-meter grid off a 912.6 MHz capture.

Itron ERT (Encoder-Receiver-Transmitter) modules on gas, water and
electric meters "bubble up" a few-millisecond burst every tens of
seconds, hopping across the 910-920 MHz ISM band.  The grid:

  chip clock   32768 chips/s  (= a 32.768 kHz watch crystal)
  bit clock    16384 bit/s    (Manchester: two chips per bit, one of
                               them always a transition at bit centre)
  modulation   OOK - the carrier is simply keyed by the chips
  frames       SCM  96 bits  (21-bit preamble 0x1F2A60, CRC-16 0x6F63)
               SCM+ 128 bits (sync 0x16A3 + 0x1E, CRC-16/CCITT 0x1021 init 0x0971)
               IDM  720 bits (sync 0x16A3 + 0x1C, CRC 0x1021 init 0xD895)

This script derives, from raw IQ:
  1. CRC-verified frames (the gold standard of proof), per frame type
  2. the chip clock, from the unrestricted level-change spectrum of
     every burst - no decoding needed to see 32768 Hz
  3. the Manchester guarantee, counted: transitions at bit centres vs
     at bit boundaries over all verified frames
  4. which Manchester polarity the air uses (both are tried, the CRC
     decides)
  5. the hop channels: every verified burst's carrier offset, and a
     comb search for the spacing
  6. optional cross-check against rtl_433 run on the same file

PRIVACY: meter IDs are masked to their last two digits and consumption
values are never printed.  Neighbours' meters are not ours to publish.

Usage: python measure.py --iq ert_912p6.cs16 --fs 2400000 [--referee rtl433.jsonl]
"""
import argparse
import collections
import json
from pathlib import Path

import numpy as np
from scipy.signal import firwin, lfilter

HERE = Path(__file__).resolve().parent
CHIP = 32768.0                    # chips per second
BIT = CHIP / 2                    # Manchester: 2 chips per bit
SCM_PRE = 0x1F2A60                # 21 bits, sync bit + preamble
SCMP_SYNC = 0x16A31E              # 24 bits: 0x16A3 + protocol 0x1E
IDM_SYNC = 0x16A31C               # 24 bits: 0x16A3 + packet type 0x1C
FRAMES = {                        # name: (sync, sync bits, frame bits)
    "SCM": (SCM_PRE, 21, 96),
    "SCM+": (SCMP_SYNC, 24, 128),
    "IDM": (IDM_SYNC, 24, 720),
}
METER_TYPE = {4: "electric", 5: "electric", 7: "electric", 8: "electric",
              0: "gas", 1: "gas", 2: "gas", 9: "gas", 12: "gas",
              3: "water", 11: "water", 13: "water"}


# ---------------------------------------------------------------- CRC
def crc16(data, poly, init):
    """MSB-first CRC-16, the same routine rtl_433 uses."""
    crc = init
    for byte in data:
        crc ^= byte << 8
        for _ in range(8):
            crc = ((crc << 1) ^ poly) & 0xFFFF if crc & 0x8000 else (crc << 1) & 0xFFFF
    return crc


def to_bytes(bits):
    return bytes(int("".join(map(str, bits[i:i + 8])), 2)
                 for i in range(0, len(bits) - len(bits) % 8, 8))


def verify(kind, bits):
    """Frame bits (sync-aligned, MSB first) -> dict or None."""
    b = to_bytes(bits)
    if kind == "SCM":
        if len(b) < 12 or crc16(b[2:12], 0x6F63, 0) != 0:
            return None
        ert_id = ((b[2] & 0x06) << 23) | (b[7] << 16) | (b[8] << 8) | b[9]
        return {"kind": kind, "id": ert_id, "type": (b[3] >> 2) & 0x0F}
    if kind == "SCM+":
        if len(b) < 16 or crc16(b[2:14], 0x1021, 0x0971) != (b[14] << 8 | b[15]):
            return None
        return {"kind": kind, "id": int.from_bytes(b[4:8], "big"), "type": b[3] & 0x0F}
    if kind == "IDM":
        if len(b) < 90 or crc16(b[2:88], 0x1021, 0xD895) != (b[88] << 8 | b[89]):
            return None
        return {"kind": kind, "id": int.from_bytes(b[7:11], "big"), "type": b[6] & 0x0F}
    return None


# ---------------------------------------------------------------- synthetic TX
def scmplus_frame(rng, ert_id=None):
    """A valid SCM+ frame as 128 bits (random ID and consumption)."""
    ert_id = int(rng.integers(1, 2 ** 32)) if ert_id is None else ert_id
    body = bytes([0x1E, 0x9C]) + ert_id.to_bytes(4, "big") \
        + int(rng.integers(0, 2 ** 24)).to_bytes(4, "big") + bytes([0x00, 0x00])
    crc = crc16(body, 0x1021, 0x0971)
    frame = bytes([0x16, 0xA3]) + body + bytes([crc >> 8, crc & 0xFF])
    return [int(c) for c in "".join(f"{x:08b}" for x in frame)]


def scm_frame(rng):
    """A valid 96-bit SCM frame: 21-bit preamble, 59 payload bits, CRC."""
    while True:
        bits = [int(c) for c in f"{SCM_PRE:021b}"] + list(rng.integers(0, 2, 59))
        b = to_bytes(bits + [0] * 16)
        crc = crc16(b[2:10], 0x6F63, 0)
        full = bits + [int(c) for c in f"{crc:016b}"]
        if verify("SCM", full):
            return full


def manchester(bits, one=(1, 0)):
    zero = (1 - one[0], 1 - one[1])
    return [c for b in bits for c in (one if b else zero)]


def ook_burst(chips, fs, f_off, amp=1.0, lead_chips=20):
    """Chips -> OOK carrier at f_off (a 32768 Hz chip grid on an fs grid)."""
    chips = [0] * lead_chips + list(chips) + [0] * lead_chips
    n = int(len(chips) * fs / CHIP)
    idx = np.minimum((np.arange(n) * CHIP / fs).astype(int), len(chips) - 1)
    env = np.asarray(chips, float)[idx]
    return amp * env * np.exp(2j * np.pi * f_off / fs * np.arange(n))


# ---------------------------------------------------------------- receiver
def burst_scan(path, fs, blk_s=10.0, win_s=100e-6):
    """Wideband envelope at 100 us -> burst (start, duration) list in seconds."""
    m = int(win_s * fs)
    env = []
    with open(path, "rb") as fh:
        while True:
            raw = np.fromfile(fh, dtype=np.int16, count=int(blk_s * fs) * 2)
            if raw.size < 2 * m:
                break
            x = raw[0::2].astype(np.float32) ** 2 + raw[1::2].astype(np.float32) ** 2
            env.append(np.sqrt(x[:x.size // m * m].reshape(-1, m).mean(1)))
    e = np.concatenate(env)
    med = np.median(e)
    thr = med + 6 * 1.4826 * np.median(np.abs(e - med))
    idx = np.flatnonzero(e > thr)
    segs = np.split(idx, np.flatnonzero(np.diff(idx) > 3) + 1) if idx.size else []
    return [(s[0] * win_s, (s[-1] - s[0] + 1) * win_s) for s in segs
            if (s[-1] - s[0] + 1) * win_s >= 2.5e-3], e.size * win_s


def load(path, fs, t0, dur):
    raw = np.fromfile(path, dtype=np.int16, offset=max(int(t0 * fs), 0) * 4,
                      count=int(dur * fs) * 2)
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0


def carrier_offset(x, fs):
    """OOK keeps half its power in the carrier line: the FFT peak is it."""
    n = 1 << 15
    F = np.abs(np.fft.fft(x[:n] * np.hanning(min(n, x.size)), n))
    k = int(np.argmax(F))
    d = 0.5 * (F[k - 1] - F[k + 1]) / (F[k - 1] - 2 * F[k] + F[k + 1] + 1e-12)
    return ((k + d) / n * fs + fs / 2) % fs - fs / 2


def channelize(x, fs, f_off, taps, dec):
    z = x * np.exp(-2j * np.pi * f_off / fs * np.arange(x.size))
    return lfilter(taps, 1.0, z)[::dec]


def slice_chips(env, fsd):
    """Envelope -> chip levels, chip phase, transition times (chips).

    Returns None unless the envelope is two-level: OOK leaves the upper
    and lower quartiles an order of magnitude apart, while an FSK burst
    (the Neptune R900 population in this band) or a burst our carrier
    estimator mis-tuned has a Rayleigh-ish envelope with hi/lo ~ 2 and
    toggles on noise at ~100 kHz.  Those must not vote on the chip clock.
    """
    hi = np.median(env[env > np.percentile(env, 75)])
    lo = np.median(env[env < np.percentile(env, 25)])
    if hi < 4 * lo:
        return None
    lv = (env > 0.5 * (hi + lo)).astype(int)
    sps = fsd / CHIP
    on = np.flatnonzero(lv)
    if on.size < 8:
        return None
    # transitions inside the burst only: the silent margins toggle on noise
    tr = np.flatnonzero(np.diff(lv[on[0]:on[-1] + 1]) != 0).astype(float) + on[0] + 0.5
    if tr.size < 8:
        return None
    ph = np.angle(np.mean(np.exp(2j * np.pi * (tr % sps) / sps)))
    start = (ph / (2 * np.pi) * sps + sps / 2) % sps
    idx = np.arange(start, env.size, sps)
    chips = lv[idx.astype(int)]
    return chips, start, tr, sps


def hunt(chips):
    """Try both Manchester pairings x both polarities; the CRC decides."""
    out = []
    for pair in (0, 1):
        pairs = chips[pair:]
        pairs = pairs[:pairs.size // 2 * 2].reshape(-1, 2)
        valid = pairs[:, 0] != pairs[:, 1]             # Manchester: chips differ
        for pol, bits in (("10", pairs[:, 0]), ("01", pairs[:, 1])):
            bstr = "".join("x" if not v else str(b) for v, b in zip(valid, bits))
            for kind, (sync, nsync, nbits) in FRAMES.items():
                pat = f"{sync:0{nsync}b}"
                i = bstr.find(pat)
                while i != -1:
                    seg = bstr[i:i + nbits]
                    if len(seg) == nbits and "x" not in seg:
                        fr = verify(kind, [int(c) for c in seg])
                        if fr:
                            fr.update(pair=pair, pol=pol, bit0=i,
                                      bits=[int(c) for c in seg])
                            out.append(fr)
                    i = bstr.find(pat, i + 1)
    return out


def process_burst(x, fs, taps, dec):
    """One burst's IQ -> (frames, carrier offset, chips, transitions)."""
    fsd = fs / dec
    f_off = carrier_offset(x, fs)
    env = np.abs(channelize(x, fs, f_off, taps, dec))
    sl = slice_chips(env, fsd)
    if sl is None:
        return [], f_off, None
    chips, start, tr, sps = sl
    return hunt(chips), f_off, (chips, start, tr, sps, env)


# ================================================================= main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--referee", default=None,
                    help="rtl_433 JSON output for the same file")
    a = ap.parse_args()
    fs = a.fs
    dec = max(1, int(fs / 300e3))
    fsd = fs / dec
    taps = firwin(129, 60e3, fs=fs)
    rng = np.random.default_rng(7)

    # ---- SYNTHETIC TX FIRST: the receiver must pass its own gauntlet
    ok = 0
    for k in range(6):
        bits = scmplus_frame(rng) if k % 2 else scm_frame(rng)
        x = ook_burst(manchester(bits), fs, rng.uniform(-900e3, 900e3), amp=0.02)
        x = np.concatenate([np.zeros(2000), x, np.zeros(2000)])
        x = x + 0.003 * (rng.standard_normal(x.size) + 1j * rng.standard_normal(x.size))
        frames, _, _ = process_burst(x, fs, taps, dec)
        ok += any(f["bits"] == list(bits) for f in frames)
    print(f"synthetic self-test: {ok}/6 CRC-verified, bit-exact recoveries "
          f"(3 SCM + 3 SCM+, random offsets within +-900 kHz)")
    assert ok == 6, "receiver failed its own synthetic gauntlet"

    # ---- burst census
    bursts, total_s = burst_scan(a.iq, fs)
    print(f"bursts >= 2.5 ms: {len(bursts)} in {total_s:.0f} s")

    NF = 1 << 16
    clk = np.zeros(NF // 2 + 1)
    recs, durs = [], []
    centre_tr, centre_n, edge_tr, edge_n = 0, 0, 0, 0
    n_ook = 0
    fig_burst = None
    for t0, dur in bursts:
        x = load(a.iq, fs, t0 - 2e-3, dur + 4e-3)
        frames, f_off, sl = process_burst(x, fs, taps, dec)
        durs.append(dur)
        if sl is None:
            continue
        n_ook += 1
        chips, start, tr, sps, env = sl
        # unrestricted level-change spectrum (no decoding involved)
        imp = np.zeros(env.size)
        imp[np.minimum(tr.astype(int), env.size - 1)] = 1.0
        clk += np.abs(np.fft.rfft(imp, NF)) ** 2
        for fr in frames:
            # the Manchester guarantee, counted on this verified frame
            c0 = fr["pair"] + 2 * fr["bit0"]
            nb = len(fr["bits"])
            seg = chips[c0:c0 + 2 * nb]
            centre_tr += int(np.sum(seg[0::2] != seg[1::2]))
            centre_n += nb
            edge_tr += int(np.sum(seg[1:-1:2] != seg[2::2]))
            edge_n += nb - 1
            recs.append({"t": t0, "f": f_off, **{k: fr[k] for k in
                         ("kind", "id", "type", "pol", "pair")}})
            if fig_burst is None:
                fig_burst = (env, start, sps, fr, fsd)

    print(f"two-level (OOK) envelopes: {n_ook}; the other {len(bursts) - n_ook} "
          f"are FSK or too weak to tune (not ERT, or not ours to claim)")
    kinds = collections.Counter(r["kind"] for r in recs)
    print(f"CRC-verified frames: {len(recs)}  {dict(kinds)}")
    pols = collections.Counter(r["pol"] for r in recs)
    print(f"Manchester polarity that passed CRC: {dict(pols)}  "
          f"('10' = bit 1 is carrier-on then carrier-off)")
    print(f"transitions at bit centres: {centre_tr}/{centre_n}   "
          f"at bit boundaries: {edge_tr}/{edge_n} "
          f"({100 * edge_tr / max(edge_n, 1):.0f}% - data-dependent)")
    d = np.array(durs) * 1e3
    print(f"burst durations (ms): SCM-length 5-7: {np.sum((d > 5) & (d < 7))}, "
          f"SCM+-length 7-9: {np.sum((d > 7) & (d < 9))}, "
          f"IDM-length 43-48: {np.sum((d > 43) & (d < 48))}, other: "
          f"{np.sum((d <= 5) | ((d >= 9) & (d <= 43)) | (d >= 48))}")

    # ---- chip clock
    fx = np.fft.rfftfreq(NF, 1 / fsd)
    band = fx > 2000
    i = int(np.argmax(clk * band))
    dd = 0.5 * (clk[i - 1] - clk[i + 1]) / (clk[i - 1] - 2 * clk[i] + clk[i + 1])
    f_chip = fx[i] + dd * (fx[1] - fx[0])
    snr = 10 * np.log10(clk[i] / np.median(clk[band]))
    print(f"chip clock: {f_chip:.1f} Hz ({snr:+.1f} dB) - published 32768 "
          f"({1e6 * (f_chip / CHIP - 1):+.0f} ppm, all meters pooled)")

    # ---- census (masked)
    per = collections.defaultdict(list)
    for r in recs:
        per[(r["kind"], r["id"], r["type"])].append(r)
    print(f"meters heard (IDs masked to last 2 digits): {len(per)}")
    for (kind, mid, mt), rs in sorted(per.items(), key=lambda kv: -len(kv[1])):
        ts = np.sort([r["t"] for r in rs])
        iv = f", median interval {np.median(np.diff(ts)):.1f} s" if ts.size > 2 else ""
        fr = sorted(round(r["f"] / 1e3) for r in rs)
        print(f"  ...{mid % 100:02d} {kind:4s} {METER_TYPE.get(mt, 'unknown'):8s} "
              f"n={len(rs)}{iv}, carrier offsets kHz {fr}")

    # ---- hop comb
    f_all = np.array([r["f"] for r in recs])
    spac = np.arange(20e3, 600e3, 250.0)
    R = np.array([np.abs(np.mean(np.exp(2j * np.pi * f_all / s))) for s in spac])
    best = spac[int(np.argmax(R))]
    resid = ((f_all / best + 0.5) % 1 - 0.5) * best
    print(f"hop comb over {f_all.size} verified bursts: best spacing "
          f"{best / 1e3:.2f} kHz (R = {R.max():.3f}), residual {resid.std() / 1e3:.1f} kHz rms; "
          f"window is only {fs / 1e6:.1f} MHz of the 910-920 MHz hop band")
    for s in (50e3, 100e3, 200e3, 250e3, 400e3):
        k = int(np.argmin(np.abs(spac - s)))
        print(f"  R at {s / 1e3:.0f} kHz spacing: {R[k]:.3f}")
    chan = np.round(f_all / best).astype(int)
    print(f"  channels occupied in this window: {sorted(set(chan.tolist()))} x {best / 1e3:.0f} kHz "
          f"(counts {dict(sorted(collections.Counter(chan.tolist()).items()))})")

    # ---- referee
    if a.referee:
        ref = [json.loads(l) for l in open(a.referee)]
        ref = [r for r in ref if "id" in r]
        rk = collections.Counter(r["model"] for r in ref)
        ours = {r["id"] for r in recs}
        theirs = {r["id"] for r in ref if r["model"] in ("ERT-SCM", "SCMplus", "IDM")}
        print(f"referee (rtl_433 on the same file): {len(ref)} frames {dict(rk)}")
        print(f"  ERT IDs: ours {len(ours)}, rtl_433 {len(theirs)}, "
              f"common {len(ours & theirs)}, ours-only {len(ours - theirs)}, "
              f"rtl_433-only {len(theirs - ours)}")

    # ---- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(14, 4.2))

    env, start, sps, fr, fsd_ = fig_burst
    c0 = fr["pair"] + 2 * fr["bit0"]
    i0 = int(start + (c0 - 6) * sps)
    i1 = int(start + (c0 + 2 * 28) * sps)
    tt = (np.arange(i0, i1) - start) / fsd_ * 1e3
    a1.plot(tt, env[i0:i1] / env.max(), lw=.9, color="#0aa")
    for k in range(c0, c0 + 2 * 28, 2):
        a1.axvline((start + k * sps - start) / fsd_ * 1e3, color="#888", lw=.4)
    nsync = FRAMES[fr["kind"]][1]
    a1.axvspan(c0 / CHIP * 1e3, (c0 + 2 * nsync) / CHIP * 1e3, color="#dd3",
               alpha=.2, label=f"{fr['kind']} sync ({nsync} bits)")
    a1.legend(loc="lower right", fontsize=8)
    a1.set(title="one verified burst, envelope\n(grey: the 16384 Hz bit grid)",
           xlabel="ms from frame start", ylabel="envelope")

    mm = (fx > 31500) & (fx < 34000)
    a2.semilogy(fx[mm], clk[mm] + 1e-9, lw=.7, color="#0aa")
    a2.axvline(CHIP, color="#d33", ls=":", lw=1)
    a2.set(title=f"level-change spectrum, all bursts\nchip clock {f_chip:.1f} Hz "
                 f"(dotted: 32768)", xlabel="Hz")

    for (kind, mid, mt), rs in per.items():
        a3.plot([r["t"] for r in rs], [r["f"] / 1e3 for r in rs], ".", ms=5)
    for k in range(int(-fs / 2 / best) - 1, int(fs / 2 / best) + 2):
        a3.axhline(k * best / 1e3, color="#888", lw=.3)
    a3.set(title=f"every verified burst: time vs carrier offset\n"
                 f"(one colour per meter; grey comb {best / 1e3:.1f} kHz)",
           xlabel="capture time (s)", ylabel="kHz from 912.6 MHz dial")
    fig.tight_layout()
    fig.savefig(out / "ert.png", dpi=110)
    print(f"figure -> {out / 'ert.png'}")


if __name__ == "__main__":
    main()
