#!/usr/bin/env python3
"""Measure the Z-Wave R2 grid (908.4 MHz, 40 kbit/s FSK) and its symbol lock from raw IQ.

Z-Wave's 40 kbit/s mode is two-tone FSK with a ±20 kHz deviation: modulation index
h = 2·20/40 = 1. At h = 1 the two tones are exactly one bit-rate apart, so within one
bit the signal turns a whole extra cycle on one tone versus the other -- orthogonal
tones. The usual receiver (an FM discriminator: the phase step between samples) is
the wrong tool here: every tone change at h = 1 swings the phase through a full
half-turn inside a narrow channel filter and the discriminator spikes. This script:

  1. finds the bursts and measures the grid off the longest one (the FLiRS wake-up
     beam: a 1.1 s train of short fragments) -- tone pair, deviation, centre, and
     the bit clock from the transitions, with a wide search window;
  2. decodes every frame twice -- discriminator slicing and non-coherent tone-energy
     (|sum x e^-j2pi(c+20k)t|^2 vs |sum x e^-j2pi(c-20k)t|^2 per bit) -- and counts
     how many pass the Z-Wave XOR checksum each way;
  3. reads the frame structure, the beam, the bit polarity and the ACK turnaround.

Home IDs are masked in all output (they identify a network).

usage: python measure.py --iq cap.npy --fs 800e3 [--offset-khz 0] [--fig figures/zwave.png]
       --iq accepts .npy (complex) or interleaved int16 (.cs16); the Z-Wave channel
       should sit near --offset-khz from the capture centre.
"""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.signal import firwin, filtfilt

HERE = Path(__file__).resolve().parent
KIND = {1: "singlecast", 2: "multicast", 3: "ACK", 8: "routed"}


def load(path, fs, off_hz):
    p = Path(path)
    if p.suffix == ".npy":
        x = np.load(p, mmap_mode="r")
        x = np.asarray(x, np.complex64) if np.iscomplexobj(x) else (x[0::2] + 1j * x[1::2]).astype(np.complex64)
    else:
        r = np.fromfile(p, np.int16).astype(np.float32)
        x = (r[0::2] + 1j * r[1::2]).astype(np.complex64)
    if off_hz:
        x = x * np.exp(-2j * np.pi * off_hz / fs * np.arange(len(x))).astype(np.complex64)
    return x


def bursts(x, fs, thr_db=4.0, gap_s=0.5e-3, min_s=0.5e-3):
    k = max(int(fs * 1e-4), 1)
    env = np.convolve(np.abs(x) ** 2, np.ones(k) / k, "same")
    floor = np.median(env)
    on = np.where(env > floor * 10 ** (thr_db / 10))[0]
    segs = []
    if on.size:
        st = last = on[0]
        for i in on[1:]:
            if i - last > gap_s * fs:
                segs.append((st, last)); st = i
            last = i
        segs.append((st, last))
    return [(a, b, 10 * np.log10(env[a:b].mean() / floor)) for a, b in segs if (b - a) >= min_s * fs]


def discriminator(x, fs, bw=45e3):
    s = filtfilt(firwin(65, bw / (fs / 2)), 1, x)
    return np.angle(s[1:] * np.conj(s[:-1])) * fs / (2 * np.pi), s


def psd(x, fs, N=4096):
    m = x[: len(x) // N * N].reshape(-1, N)
    P = np.mean(np.abs(np.fft.fftshift(np.fft.fft(m * np.hanning(N), axis=1), axes=1)) ** 2, 0)
    f = np.fft.fftshift(np.fft.fftfreq(N, 1 / fs))
    return f, 10 * np.log10(P / np.median(P))


def measure_grid(x, bl, fs):
    """Tones from the data frames' spectrum, comb spacing from the wake-up beam,
    bit clock from the discriminator's transition spikes on the data frames."""
    fr = [(a, b) for a, b, _ in bl if 3e-3 * fs < (b - a) < 30e-3 * fs]
    beam = max(bl, key=lambda s: s[1] - s[0])
    f, Pd = psd(np.concatenate([x[a:b] for a, b in fr]), fs)
    sel = np.abs(f) < 150e3
    i1 = np.where(sel)[0][np.argmax(Pd[sel])]
    far = sel & (np.abs(f - f[i1]) > 10e3)
    i2 = np.where(far)[0][np.argmax(Pd[far])]
    lo, hi = sorted((f[i1], f[i2]))
    fb, Pb = psd(x[beam[0]: beam[1]], fs)
    selb = np.abs(fb) < 150e3
    lines = fb[selb][Pb[selb] > Pb[selb].max() - 15]
    lines = np.array([g.mean() for g in np.split(lines, np.where(np.diff(lines) > 2 * fs / 4096)[0] + 1)])
    comb = np.median(np.diff(lines)) if len(lines) > 2 else float("nan")
    T = None
    for a, b in fr:
        d, _ = discriminator(x[a:b], fs, bw=60e3)
        e = (np.abs(d - np.median(d)) > 45e3).astype(np.float32)
        E = np.abs(np.fft.rfft(e - e.mean(), 1 << 18)) ** 2
        T = E if T is None else T + E
    ft = np.fft.rfftfreq(1 << 18, 1 / fs)
    w = (ft > 5e3) & (ft < 200e3)
    k = np.where(w)[0][np.argmax(T[w])]
    return dict(f=fb, Pd=Pb, ff=f, Pf=Pd, lo=lo, hi=hi, centre=(lo + hi) / 2, dev=(hi - lo) / 2, comb=comb,
                clock=ft[k], clock_snr=10 * np.log10(T[k] / np.median(T[w])), nfr=len(fr),
                beam_ms=1e3 * (beam[1] - beam[0]) / fs)


def bits_tone(x, fs, c, dev, rate):
    sps = int(round(fs / rate))
    n = np.arange(len(x))
    up = x * np.exp(-2j * np.pi * (c + dev) / fs * n)
    dn = x * np.exp(-2j * np.pi * (c - dev) / fs * n)
    best = None
    for ph in range(sps):
        k = (len(x) - ph) // sps
        U = np.abs(up[ph: ph + k * sps].reshape(k, sps).sum(1)) ** 2
        D = np.abs(dn[ph: ph + k * sps].reshape(k, sps).sum(1)) ** 2
        q = np.mean(np.abs(U - D) / (U + D + 1e-12))
        if best is None or q > best[0]:
            best = (q, (U > D).astype(np.uint8), ph)
    return best


def bits_disc(x, fs, c, rate, bw=45e3, clip=None):
    d, _ = discriminator(x, fs, bw)
    d = d - c
    if clip:
        d = np.clip(d, -clip, clip)
    sps = fs / rate
    best = None
    for ph in range(int(sps)):
        k = int((len(d) - ph) / sps)
        I = np.array([d[int(ph + i * sps): int(ph + (i + 1) * sps)].mean() for i in range(k)])
        q = np.abs(I).mean()
        if best is None or q > best[0]:
            best = (q, (I > 0).astype(np.uint8), ph)
    return best


DISC_VARIANTS = [(45e3, None), (45e3, 40e3), (60e3, 40e3), (30e3, 30e3)]


def frames(bits):
    out = []
    for inv in (0, 1):
        s = "".join(map(str, bits ^ inv))
        i = 0
        while True:
            i = s.find("0101010111110000", i)
            if i < 0:
                break
            j = i + 16
            bb = (bits ^ inv)[j: j + (len(bits) - j) // 8 * 8]
            raw = np.packbits(bb) if len(bb) >= 8 else np.zeros(0, np.uint8)
            if len(raw) >= 3 and raw[0] == 0x55:
                out.append(dict(inv=inv, beam=True, node=int(raw[1]), bit=i)); i = j + 24; continue
            if len(raw) >= 10 and 10 <= raw[7] <= len(raw):
                L = int(raw[7]); fr = raw[:L]
                cs = 0xFF
                for v in fr[:-1]:
                    cs ^= int(v)
                if cs == fr[-1]:
                    fc = int(fr[5]) << 8 | int(fr[6])
                    out.append(dict(inv=inv, beam=False, home=fr[:4].tobytes().hex(), src=int(fr[4]),
                                    dst=int(fr[8]), kind=KIND.get(fc >> 8 & 0xF, "h%d" % (fc >> 8 & 0xF)),
                                    ackreq=bool(fc & 0x4000), seq=fc & 0xF, payload=fr[9:-1].tobytes().hex(),
                                    L=L, bit=i))
                    i = j + 8 * L; continue
            i = j
    return out


def mask(h):
    return "********"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", nargs="+", required=True)
    ap.add_argument("--fs", type=float, default=800e3)
    ap.add_argument("--offset-khz", type=float, default=0.0)
    ap.add_argument("--fig", default=str(HERE / "figures" / "zwave.png"))
    a = ap.parse_args()
    fs = a.fs
    allf = []
    figdata = {}
    for path in a.iq:
        x = load(path, fs, a.offset_khz * 1e3)
        bl = bursts(x, fs)
        print(f"\n== {Path(path).name}: {len(x)/fs:.1f} s, {len(bl)} bursts")
        g = measure_grid(x, bl, fs)
        print(f"  data frames ({g['nfr']} bursts 3-30 ms): tones {g['lo']/1e3:+.2f} / {g['hi']/1e3:+.2f} kHz"
              f" -> centre {g['centre']/1e3:+.2f} kHz, deviation ±{g['dev']/1e3:.2f} kHz")
        print(f"  wake-up beam ({g['beam_ms']:.1f} ms, mostly 0101 preamble): comb spacing {g['comb']/1e3:.3f} kHz"
              f" -> bit rate 2x = {2*g['comb']/1e3:.2f} kbit/s")
        print(f"  transition-spike line on the data frames (searched 5-200 kHz): {g['clock']/1e3:.3f} kHz (+{g['clock_snr']:.1f} dB)")
        rate = 40e3
        h = 2 * g["dev"] / rate
        print(f"  modulation index h = 2·dev/rate = {h:.3f}")
        figdata.setdefault("grid", g)
        n_tone = 0; n_disc = [0] * len(DISC_VARIANTS); maxbits = 0
        beams = Counter(); evs = []
        for b0, b1, snr in bl:
            seg = x[max(b0 - 40, 0): b1 + 40]
            dd = np.clip(discriminator(seg, fs)[0], -40e3, 40e3)
            c = (np.percentile(dd, 20) + np.percentile(dd, 80)) / 2
            qt, bt, pht = bits_tone(seg, fs, c, g["dev"], rate)
            ft = [f for f in frames(bt) if not f["beam"]]
            n_tone += len(ft)
            for vi, (bw, clip) in enumerate(DISC_VARIANTS):
                n_disc[vi] += len([f for f in frames(bits_disc(seg, fs, c, rate, bw, clip)[1]) if not f["beam"]])
            for f in frames(bt):
                if f["beam"]:
                    beams[(f["inv"], f["node"])] += 1
            for f in ft:
                t0 = (max(b0 - 40, 0) + pht + f["bit"] * fs / rate) / fs
                t1 = t0 + (16 + 8 * f["L"]) / rate
                maxbits = max(maxbits, 16 + 8 * f["L"])
                evs.append((t0, t1, f, snr, c, b0 / fs, b1 / fs))
        print(f"  frames passing the checksum: tone-energy {n_tone}; FM discriminator (filter, clip) "
              + ", ".join(f"{bw/1e3:.0f}k/{('±%.0fk' % (cl/1e3)) if cl else 'none'}: {n}" for (bw, cl), n in zip(DISC_VARIANTS, n_disc)))
        if maxbits:
            print(f"  checksum as a ruler: a {maxbits}-bit frame decodes on a FIXED 40000 bit/s clock with no tracking,"
                  f" so the true rate is within ±{40e3/(2*maxbits):.0f} bit/s (±{100/(2*maxbits):.2f} %)")
        if beams:
            (inv, node), cnt = beams.most_common(1)[0]
            print(f"  beam fragments: {sum(beams.values())} (most common: node {node}, {cnt}x)")
        pol = Counter(e[2]["inv"] for e in evs)
        print(f"  polarity: {dict(pol)}  (inv=1 means a logical 1 is the LOWER tone)")
        evs.sort(key=lambda e: e[0])
        for t0, t1, f, snr, c, *_ in evs:
            print(f"   {t0:8.4f}s +{snr:4.1f}dB c={c/1e3:+5.1f}k home={mask(f['home'])} {f['src']:>3}->{f['dst']:<3} "
                  f"{f['kind']:<10} seq={f['seq']:<2} ackreq={int(f['ackreq'])} payload={f['payload'][:12]}{'…' if len(f['payload'])>12 else ''}")
        # ACK turnaround: ackreq frame A->B, next frame is the ACK B->A.
        # air gap = ACK burst start - frame's last bit; then the ACK spends its preamble before SOF.
        ta = []
        for ea, eb in zip(evs, evs[1:]):
            fa, fb = ea[2], eb[2]
            if fa["ackreq"] and fb["kind"] == "ACK" and fb["src"] == fa["dst"] and fb["dst"] == fa["src"]:
                gap = (eb[5] - ea[1]) * 1e3
                pre = (eb[0] - eb[5]) * 1e3
                ta.append((fa["dst"], gap, pre))
        for node, gap, pre in ta:
            print(f"  ACK by node {node}: air gap {gap:.2f} ms after the last bit, then {pre:.2f} ms of preamble before its SOF")
        allf.append((Path(path).name, evs, ta, x, bl, g))
    make_fig(allf, fs, a.fig)


def make_fig(allf, fs, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # panels a-c from the capture with the longest wake-up beam, panel d from one with a nonce exchange
    name, evs, ta, x, bl, g = max(allf, key=lambda r: r[5]["beam_ms"])
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    # (a) beam PSD
    a0 = ax[0, 0]
    sel = np.abs(g["f"]) < 120e3
    a0.plot(g["f"][sel] / 1e3, g["Pd"][sel], lw=0.8, label="wake-up beam (0101…)")
    a0.plot(g["ff"][sel] / 1e3, g["Pf"][sel], lw=0.8, c="tab:orange", label="data frames")
    for v, lab in ((g["lo"], "logical 1"), (g["hi"], "logical 0")):
        a0.axvline(v / 1e3, c="r", ls="--", lw=0.8); a0.text(v / 1e3, max(g["Pd"][sel].max(), g["Pf"][sel].max()) + 1, lab, ha="center", fontsize=8, c="r")
    a0.legend(fontsize=8, loc="upper right")
    a0.set_title(f"tones {g['lo']/1e3:+.1f}/{g['hi']/1e3:+.1f} kHz (frames); beam comb {g['comb']/1e3:.2f} kHz = rate/2")
    a0.set_xlabel("kHz from 908.40 MHz"); a0.set_ylabel("dB over median"); a0.grid(alpha=.3)
    # (b) discriminator vs tone energy on one frame
    a1 = ax[0, 1]
    t0, t1, f, snr, c = next(e for e in evs if e[2]["kind"] == "singlecast")[:5]
    i0 = int(t0 * fs) + int(20 / 40e3 * fs); n = int(12 / 40e3 * fs)
    seg = x[i0 - 200: i0 + n + 200]
    d, _ = discriminator(seg, fs)
    d = d[200: 200 + n] - c
    tt = np.arange(n) / fs * 1e6
    a1.plot(tt, np.clip(d, -90e3, 90e3) / 1e3, lw=0.7, label="FM discriminator (kHz)")
    nn = np.arange(len(seg))
    up = np.abs(np.convolve(seg * np.exp(-2j * np.pi * (c + g["dev"]) / fs * nn), np.ones(20), "same"))[200: 200 + n]
    dn = np.abs(np.convolve(seg * np.exp(-2j * np.pi * (c - g["dev"]) / fs * nn), np.ones(20), "same"))[200: 200 + n]
    sc = 60 / max(up.max(), dn.max())
    a1.plot(tt, (up - dn) * sc, lw=1.2, label="tone energy: upper − lower (scaled)")
    for k in range(13):
        a1.axvline(k * 25, c="k", lw=0.3, alpha=.4)
    a1.set_title("12 bits of one frame: the discriminator spikes at h = 1, tone energy doesn't")
    a1.set_xlabel("µs"); a1.legend(fontsize=8, loc="lower right"); a1.set_ylim(-95, 95); a1.grid(alpha=.3)
    # (c) timing eye metric vs phase, one frame
    a2 = ax[1, 0]
    seg = x[int(t0 * fs) - 40: int(t1 * fs) + 40]
    nn = np.arange(len(seg)); sps = 20
    U0 = seg * np.exp(-2j * np.pi * (c + g["dev"]) / fs * nn); D0 = seg * np.exp(-2j * np.pi * (c - g["dev"]) / fs * nn)
    qs = []
    for ph in range(sps):
        k = (len(seg) - ph) // sps
        U = np.abs(U0[ph: ph + k * sps].reshape(k, sps).sum(1)) ** 2
        D = np.abs(D0[ph: ph + k * sps].reshape(k, sps).sum(1)) ** 2
        qs.append(np.mean(np.abs(U - D) / (U + D)))
    a2.plot(np.arange(sps) / fs * 1e6, qs, "o-")
    a2.set_title("symbol lock: tone contrast vs integration start (one 25 µs bit)")
    a2.set_xlabel("start offset in µs"); a2.set_ylabel("mean |U−D|/(U+D)"); a2.grid(alpha=.3)
    # (d) one device-initiated exchange (starts with an S2 nonce get, 9f01)
    a3 = ax[1, 1]
    pick = None
    for r in allf:
        for e in r[1]:
            if e[2]["kind"] == "singlecast" and e[2]["payload"].startswith("9f01"):
                pick = (r[1], e[0]); break
        if pick:
            break
    if pick:
        evs2, first = pick
        win = [e for e in evs2 if 0 <= e[0] - first < 0.1]
        nodes = sorted({e[2]["src"] for e in win})
        ymap = {nd: i for i, nd in enumerate(nodes)}
        for t0, t1, f, snr, c, bs, be in win:
            col = "tab:green" if f["kind"] == "ACK" else "tab:blue"
            a3.barh(ymap[f["src"]], (be - bs) * 1e3, left=(bs - first) * 1e3, color=col, height=0.45, alpha=.35)
            a3.barh(ymap[f["src"]], (t1 - t0) * 1e3, left=(t0 - first) * 1e3, color=col, height=0.45)
            lab = "ACK" if f["kind"] == "ACK" else f["payload"][:4]
            a3.text((t0 - first) * 1e3, ymap[f["src"]] + 0.3, f"{lab}→{f['dst']}", fontsize=8)
        a3.set_yticks(range(len(nodes))); a3.set_yticklabels([f"node {n}" for n in nodes])
        a3.set_ylim(-0.6, len(nodes) - 0.2)
    a3.set_xlabel("ms (pale = preamble, solid = SOF..checksum)")
    a3.set_title("device-initiated S2 exchange: nonce get, nonce report, encrypted report")
    a3.grid(alpha=.3)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=90)
    print(f"\nfigure: {out}")


if __name__ == "__main__":
    main()
