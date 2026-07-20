#!/usr/bin/env python3
"""Measure the HD Radio AM (NRSC-5 all-digital MA3) OFDM grid off the air.

The grid (all anchored to 744,187.5 / 16 = 46,511.71875 Hz — the same
elemental clock as FM IBOC, divided down; AM spacing = FM spacing / 2):
  FFT length      256 samples  -> subcarrier spacing 181.6864 Hz
  guard interval   14 samples  (raised-sine TAPERED, overlap-added)
  symbol          270 samples  -> 172.265625 symbols/s
  MA3 subcarriers: carrier at 0, BPSK reference at +-1, QAM64 primary
                   +-2..26, PIDS +-27, reduced-power QAM64 "enhanced"
                   secondary/tertiary +-28..52  (edge +-9.63 kHz)
  block = 32 symbols (5.3833 Hz), L1 frame = 8 blocks = 1.48607 s
  (constants cross-read from the nrsc5 source: FFT_AM/CP_AM/needle/psmi)

Measurements (each validated on a synthetic MA3 transmitter first):
 1. PSD: carrier line + a flat +-4.9 kHz OFDM pedestal, no analog
    sidebands, no hybrid (MA1) shelves at +-10..15 kHz.
 2. Tu from CP autocorrelation |sum x(t)x*(t+tau)| — peak at the useful
    symbol length (search bounded by physics only: 2.5..10 ms).
 3. Symbol rate from the c(t) = x(t)x*(t+Tu) line comb: c(t) is
    periodic at 1/Ts; harmonics regress to 172.2656 Hz to ~mHz.
    Sub-lines at multiples of 5.3833 Hz betray the 32-symbol block.
 4. Symbol-timing sweep: FFT windows of Tu, hopped by Ts; the +-1
    reference subcarriers' symbol-to-symbol products turn +-real at
    exactly one window offset — the grid appears.
 5. Reference-bit decode: the 32-symbol block-sync word from the nrsc5
    source, the block counter cycling 0..7, and the psmi field —
    the broadcast itself says "2" = MA3 all-digital.

Hard-won laws baked in (each cost a round of synthetic debugging):
  - Notch the carrier before the timing sweep: its leakage into the
    +-1 bins fakes BPSK coherence at every offset and steals the lock.
  - Squared-BPSK "ref lines" do NOT survive per-symbol phase restarts
    (the FM law in a new disguise) — c(t) is the honest clock.
  - The PSD comb is invisible at 256/270 = 94.8% occupancy; don't
    hunt for it.

Usage: python measure.py --iq capture.cs16 --fs 500000
(60 s captured at 500 kS/s centered on the AM channel, int16 I/Q.)
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NATIVE = 46511.71875            # 744187.5 / 16
SPACING = NATIVE / 256          # 181.6864013671875 Hz
SYMRATE = NATIVE / 270          # 172.265625 Hz
BLOCKRATE = SYMRATE / 32        # 5.38330078125 Hz
# 32-symbol block-sync pattern on the reference subcarriers
# (find_block_am in the nrsc5 source; -1 = variable bit)
NEEDLE = np.array([0, 1, 1, 0, 0, 1, 0, -1, -1, 1, -1, -1, -1, -1, 0, -1,
                   -1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1])


def load(path, fs, secs=None):
    count = -1 if secs is None else 2 * int(secs * fs)
    raw = np.fromfile(path, dtype=np.int16, count=count)
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0


def parabolic(vals, i, dx):
    if 0 < i < len(vals) - 1:
        y0, y1, y2 = vals[i - 1], vals[i], vals[i + 1]
        return dx * 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
    return 0.0


def bpsk_metric(d):
    """|sum d^2| / sum |d^2|: 1 for locked (+-real) products, ~0 noise;
    rotation-invariant, so residual CFO cannot wrap it."""
    d2 = d.astype(np.complex128) ** 2
    return float(np.abs(d2.sum()) / (np.abs(d2).sum() + 1e-12))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True)
    a = ap.parse_args()
    fs = a.fs
    assert abs(fs - 500000) < 1, "this script expects 500,000 S/s"
    x = load(a.iq, fs)
    print(f"{len(x)/fs:.1f} s loaded")

    # ---- 1. PSD: pedestal, edges, carrier ---------------------------
    N = 2 ** 19
    seg = x[: len(x) // N * N].reshape(-1, N)
    P = (np.abs(np.fft.fft(seg * np.hanning(N).astype(np.float32), axis=1)) ** 2).mean(0)
    fax = np.fft.fftfreq(N, 1 / fs)
    db = 10 * np.log10(P + 1e-12)
    ped = float(np.median(db[(np.abs(fax) > 500) & (np.abs(fax) < 4500)]))
    flo = float(np.median(db[(np.abs(fax) > 15500) & (np.abs(fax) < 19000)]))
    ma1 = float(np.median(db[(np.abs(fax) > 10500) & (np.abs(fax) < 14500)]))
    sm = 10 * np.log10(np.convolve(P, np.ones(33) / 33, "same") + 1e-12)
    edge = {}
    for sgn in (+1, -1):
        m = (sgn * fax > 3000) & (sgn * fax < 7000)
        ff, dd = sgn * fax[m], sm[m]
        o = np.argsort(ff)
        below = np.where(dd[o] < (ped + flo) / 2)[0]
        edge[sgn] = float(ff[o][below[0]]) if len(below) else float("nan")
    print(f"OFDM pedestal +{ped-flo:.1f} dB over floor; edges +{edge[+1]:.0f} / "
          f"-{edge[-1]:.0f} Hz (PIDS at 27 x {SPACING:.1f} = {27*SPACING:.0f})")
    print(f"MA1 primary region 10.5-14.5 kHz: {ma1-flo:+.1f} dB vs floor "
          f"-> {'HYBRID SHELVES' if ma1-flo > 3 else 'no hybrid shelves (not MA1)'}")

    # ---- 2. carrier line (CFO) --------------------------------------
    nc = min(2 ** 24, len(x))
    C = np.abs(np.fft.fft(x[:nc] * np.hanning(nc)))
    fc = np.fft.fftfreq(nc, 1 / fs)
    m = np.abs(fc) < 100
    i = int(np.argmax(C[m]))
    cfo = float(fc[m][i]) + parabolic(np.log(C[m] + 1e-30), i, fs / nc)
    print(f"carrier line: {cfo:+.2f} Hz from dial (residual analog carrier "
          f"is retained even in all-digital MA3)")

    # ---- 3. prep: derotate; filtered copies -------------------------
    n = len(x)
    x = x * np.exp(-2j * np.pi * cfo / fs * np.arange(n)).astype(np.complex64)
    X = np.fft.fft(x)
    f = np.fft.fftfreq(n, 1 / fs)
    Xs = X.copy()
    Xs[np.abs(f) < 40] = 0          # carrier notch — see law #1
    Xs[np.abs(f) > 6000] = 0        # channel filter (neighbors out)
    xs = np.fft.ifft(Xs).astype(np.complex64)
    del Xs
    X[np.abs(f) > 11000] = 0        # wide copy keeps the enhanced partitions
    xw = np.fft.ifft(X).astype(np.complex64)
    del X, x

    # ---- 4. Tu candidates via CP autocorrelation --------------------
    # Physics bound ONLY (spacing 100..400 Hz <-> Tu 2.5..10 ms). The
    # bound is honest, so junk peaks can enter: on this very capture a
    # co-channel analog skywave station's program-audio autocorrelation
    # put a bigger peak at 2.63 ms than the guard did at 5.50 ms. The
    # discriminator is step 5: only the TRUE Tu makes c(t) periodic at
    # a symbol rate — candidates that don't comb are thrown out.
    chunk = 2 ** 21
    lo, hi = int(2.5e-3 * fs), int(10e-3 * fs)
    acc = np.zeros(hi + 2, dtype=np.complex128)
    nch = 0
    for s in range(0, len(xs) - chunk + 1, chunk):
        F = np.fft.fft(xs[s:s + chunk], 2 * chunk)
        acc += np.fft.ifft(np.abs(F) ** 2)[: hi + 2]
        nch += 1
    ac = np.abs(acc[lo:hi]) / nch
    acmed = float(np.median(ac))
    cands = []
    tmp = ac.copy()
    for _ in range(6):
        i = int(np.argmax(tmp))
        if tmp[i] < 2 * acmed:
            break
        cands.append((lo + i) / fs + parabolic(ac, i, 1 / fs))
        tmp[max(0, i - int(3e-4 * fs)): i + int(3e-4 * fs)] = 0
    print("CP-autocorr Tu candidates (ms): "
          + ", ".join(f"{t*1e3:.4f} (x{ac[int(round(t*fs))-lo]/acmed:.1f} med)"
                      for t in cands))

    # ---- 5. the c(t) comb decides which candidate is the guard ------
    def comb_fit(tau):
        lagn = int(round(tau * fs))
        c = xs[:-lagn] * np.conj(xs[lagn:])
        ncc = min(2 ** 24, len(c))
        S = np.abs(np.fft.fft(c[:ncc] * np.hanning(ncc))) ** 2
        fS = np.fft.fftfreq(ncc, 1 / fs)
        m = (fS > 50) & (fS < 950)
        fSm, dbS = fS[m], 10 * np.log10(S[m] + 1e-30)
        dbS -= np.median(dbS)
        lines, i, dfa = [], 0, fs / ncc
        while i < len(fSm):
            if dbS[i] > 12:
                j = i
                while j < len(fSm) and dbS[j] > 12:
                    j += 1
                k = i + int(np.argmax(dbS[i:j]))
                lines.append((fSm[k] + parabolic(dbS, k, dfa), float(dbS[k])))
                i = j + int(5 / dfa)
            else:
                i += 1
        if not lines:
            return None
        pos = np.array([p for p, _ in lines])
        f1 = pos[np.argmax([s for _, s in lines])]
        fam, cand = None, None
        for div in (1, 2, 3):   # strongest line can be a harmonic
            fm = np.abs(pos - np.round(pos / (f1 / div)) * (f1 / div)) < 0.5
            if fam is None or fm.sum() > fam.sum():
                fam, cand = fm, f1 / div
        h = np.round(pos / cand)
        A = np.vstack([h[fam], np.ones(fam.sum())]).T
        SRf = float(np.linalg.lstsq(A, pos[fam], rcond=None)[0][0])
        score = float(sum(s for (p, s), ok in zip(lines, fam) if ok))
        return SRf, int(fam.sum()), score, lines, fam, (fSm, dbS)

    best = None
    for tau in cands:
        r = comb_fit(tau)
        n_fam = 0 if r is None else r[1]
        sc_ = 0 if r is None else r[2]
        print(f"  tau {tau*1e3:.4f} ms: {n_fam} comb lines"
              + (f", fundamental {r[0]:.4f} Hz, {sc_:.0f} dB total" if r else ""))
        if r and (best is None or (r[1], r[2]) > (best[1][1], best[1][2])):
            best = (tau, r)
    tau0, (SR, nfam, _, lines, fam, (fSm, dbS)) = best
    pos = np.array([p for p, _ in lines])
    print(f"Tu = {tau0*1e3:.4f} ms -> 1/Tu = {1/tau0:.2f} Hz "
          f"(published {SPACING:.4f})")
    print(f"symbol-rate comb: {SR:.4f} Hz over {nfam} harmonics "
          f"(published {SYMRATE:.6f}, {(SR/SYMRATE-1)*1e6:+.1f} ppm)")
    for p, s in lines:
        if not fam[np.where(pos == p)[0][0]]:
            k = p / BLOCKRATE
            if abs(k - round(k)) < 0.02:
                print(f"  sub-line {p:.3f} Hz = {k:.3f} x block rate "
                      f"{BLOCKRATE:.4f} Hz (the 32-symbol block, seen blind)")
            else:
                print(f"  unrelated line at {p:.3f} Hz (+{s:.0f} dB) -- "
                      f"not on the symbol or block comb")
    print(f"Tu/Ts = {SR*tau0:.4f} (256/270 = {256/270:.5f}) -> grid-snapped "
          f"spacing = SR*270/256 = {SR*270/256:.4f} Hz (published {SPACING:.6f})")

    # ---- 6. symbol-timing sweep on the reference comb ---------------
    # windows of exactly Tu (2752 samples at 500k), hopped by Ts
    # (2902.336 samples, rounded per symbol so no drift accumulates):
    # the +-1 ref bins' symbol-to-symbol products go +-real at ONE offset.
    W = int(round(1 / SPACING * fs))            # 2752
    TS = fs / SYMRATE                           # 2902.336...
    n_sym = 128
    kk = np.round(np.arange(n_sym) * TS).astype(np.int64)
    twiddle = np.exp(-2j * np.pi * np.arange(W) / W).astype(np.complex64)

    def refs_at(base):
        segs = xs[base + kk[:, None] + np.arange(W)[None, :]]
        r_up = segs @ twiddle          # bin +1... (sign of exponent): bin 1
        r_dn = segs @ np.conj(twiddle)  # bin -1
        return r_up, r_dn

    step = 4
    offs = np.arange(0, int(round(TS)), step)
    scores = np.empty(len(offs))
    base0 = int(2 * fs)
    for i, off in enumerate(offs):
        u, d = refs_at(base0 + off)
        scores[i] = bpsk_metric(np.concatenate([u[1:] * np.conj(u[:-1]),
                                                d[1:] * np.conj(d[:-1])]))
    bi = int(np.argmax(scores))
    best, mbest = int(offs[bi]), scores[bi]
    for off in range(max(0, best - step), best + step + 1):
        u, d = refs_at(base0 + off)
        mm = bpsk_metric(np.concatenate([u[1:] * np.conj(u[:-1]),
                                         d[1:] * np.conj(d[:-1])]))
        if mm > mbest:
            best, mbest = off, mm
    print(f"timing sweep: lock at offset {best}/{int(round(TS))} "
          f"(metric {mbest:.3f}; sweep median {np.median(scores):.3f})")

    # ---- 7. reference bits: block sync, counter, psmi ---------------
    n_all = int((len(xw) - base0 - best - W) // TS) - 2
    kk2 = np.round(np.arange(n_all) * TS).astype(np.int64)
    segs = xw[base0 + best + kk2[:, None] + np.arange(W)[None, :]]
    car = segs.mean(1)                       # bin 0 = carrier (phase ref)
    u = (segs @ twiddle) * np.exp(-1j * np.angle(car))
    bits = (u.imag > 0).astype(np.int8)
    fixed = np.where(NEEDLE >= 0)[0]
    bestblk = (0, 0, -1.0)
    for pol in (0, 1):
        b = bits ^ pol
        for ph in range(32):
            nb = (len(b) - ph) // 32
            blocks = b[ph:ph + nb * 32].reshape(nb, 32)
            hits = float((blocks[:, fixed] == NEEDLE[fixed]).all(1).mean())
            if hits > bestblk[2]:
                bestblk = (pol, ph, hits, blocks)
    pol, ph, hits, blocks = bestblk
    bc = (blocks[:, 17] << 2) | (blocks[:, 18] << 1) | blocks[:, 19]
    psmi = [(b[26] << 4) | (b[27] << 3) | (b[28] << 2) | (b[29] << 1) | b[30]
            for b in blocks if ((b[17] << 2) | (b[18] << 1) | b[19]) == 0]
    vals, cnts = np.unique(np.array(psmi), return_counts=True)
    print(f"block sync word: hit rate {hits:.2f} over {len(blocks)} blocks "
          f"(32 symbols = {32/SYMRATE*1e3:.2f} ms each)")
    print(f"block counter: {list(map(int, bc[:16]))} ... cycles 0..7 -> "
          f"L1 frame = {256/SYMRATE:.4f} s")
    print(f"psmi: {dict(zip(map(int, vals), map(int, cnts)))} "
          f"(1 = MA1 hybrid, 2 = MA3 all-digital)")

    # ---- 8. per-subcarrier power profile ----------------------------
    nprof = min(n_all, 6000)
    prof = np.zeros(4 * W)
    for s in range(0, nprof, 500):
        F = np.fft.fft(segs[s:s + 500], n=4 * W, axis=1)
        prof += (np.abs(F) ** 2).sum(0)
    prof /= nprof
    fprof = np.fft.fftfreq(4 * W, 1 / fs)
    o = np.argsort(fprof)
    fprof, prof = fprof[o], 10 * np.log10(prof[o] + 1e-12)

    def at(i):
        return float(prof[np.argmin(np.abs(fprof - i * SPACING))])
    # ambient reference: IN-band (inside the 11 kHz filter), outside the
    # MA3 edge, dodging the +-10 kHz neighbor carriers at index +-55.04.
    # (A stopband reference flatters every number — been there.)
    nse = np.median([at(s * (i + 0.5)) for i in range(56, 60) for s in (1, -1)])
    prim = np.mean([at(s * i) for i in range(2, 27) for s in (1, -1)]) - nse
    pids = (at(27) + at(-27)) / 2 - nse
    enh = np.mean([at(s * i) for i in range(28, 53) for s in (1, -1)]) - nse
    out_ = np.mean([at(s * i) for i in (54, 56, 57, 58, 59) for s in (1, -1)]) - nse
    print(f"partitions vs in-band ambient: primary +-2..26: {prim:+.1f} dB, "
          f"PIDS +-27: {pids:+.1f} dB, enhanced +-28..52: {enh:+.1f} dB, "
          f"beyond +-53: {out_:+.1f} dB")
    if enh < 2:
        print("  (enhanced partitions are at/below the nighttime ambient "
              "here -> present per spec at -15 dB, but NOT claimed measured)")

    # ---- figures -----------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    outdir = HERE / "figures"
    outdir.mkdir(exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    mz = np.abs(fax) < 16000
    o = np.argsort(fax[mz])
    ax1.plot(fax[mz][o] / 1e3, db[mz][o], lw=.5, color="#0aa")
    for s in (+1, -1):
        ax1.axvspan(s * .18, s * 4.82, color="#d33", alpha=.10)
        ax1.axvspan(s * 5.0, s * 9.63, color="#d90", alpha=.10)
    ax1.annotate("carrier", (0, db[mz].max() - 1), ha="center", fontsize=9)
    ax1.annotate("primary\n(QAM64)", (2.6, ped + 6), ha="center", fontsize=8, color="#d33")
    ax1.annotate("enhanced\n(-7 dB)", (7.3, ped - 2), ha="center", fontsize=8, color="#b70")
    ax1.set_xlabel("kHz from 820 kHz dial")
    ax1.set_ylabel("PSD (dB)")
    ax1.set_title("all-digital AM: no analog sidebands, just grid")
    ax2.plot(offs, scores, lw=.6, color="#0aa")
    ax2.axvline(best, color="#d33", ls=":")
    ax2.annotate(f"the grid appears:\noffset {best}", (best, mbest),
                 fontsize=9, color="#d33")
    ax2.set_xlabel(f"FFT window offset (samples; one symbol = {TS:.1f})")
    ax2.set_ylabel("reference-comb BPSK metric")
    ax2.set_title("symbol-timing sweep: one true offset per symbol")
    fig.tight_layout()
    fig.savefig(outdir / "hdam_grid.png", dpi=110)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax1.plot(fSm, dbS, lw=.5, color="#0aa")
    for hh in range(1, 6):
        ax1.axvline(hh * SYMRATE, color="#d33", ls=":", lw=.8)
    for p, s in lines:
        if not fam[np.where(pos == p)[0][0]]:
            ax1.axvline(p, color="#b70", ls=":", lw=.8)
    ax1.set_xlabel("Hz")
    ax1.set_ylabel("dB over median")
    ax1.set_title("c(t)=x(t)x*(t+Tu) spectrum: 172.266 Hz symbol comb\n"
                  "(orange: sub-lines on the 5.383 Hz block-rate grid)")
    idx = fprof / SPACING
    mi = np.abs(idx) < 60
    ax2.plot(idx[mi], prof[mi] - nse, lw=.6, color="#0aa")
    for a_, b_, c_ in [(2, 26, "#d33"), (28, 52, "#b70")]:
        for s in (1, -1):
            ax2.axvspan(s * (a_ - .5), s * (b_ + .5), color=c_, alpha=.10)
    for s in (1, -1):
        ax2.annotate("neighbor\ncarrier", (s * 55.04, at(s * 55) - nse),
                     fontsize=7, ha="center", color="#666")
    ax2.set_xlabel("subcarrier index (x 181.686 Hz)")
    ax2.set_ylabel("dB over in-band ambient")
    ax2.set_title("locked per-band power: MA3 core stands, enhanced drowns")
    fig.tight_layout()
    fig.savefig(outdir / "hdam_structure.png", dpi=110)
    print(f"figures -> {outdir}")


if __name__ == "__main__":
    main()
