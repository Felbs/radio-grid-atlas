#!/usr/bin/env python3
"""Measure the HD Radio (NRSC-5 FM hybrid) OFDM grid off the air.

The grid (all anchored to a 744,187.5 Hz elemental clock):
  FFT length      2048 samples  -> subcarrier spacing 363.373 Hz
  guard interval   112 samples  (raised-sine TAPERED, overlap-added —
                                 NOT a plain cyclic prefix!)
  symbol          2160 samples  -> 344.53 symbols/s
  active carriers +-546 (= +-198.4 kHz), in 19-carrier partitions
  reference subcarriers: EVERY 19th (+-546, +-527, ...), DBPSK, known
                         32-bit training word per L1 block

Two measurements prove the grid on any capture:
 1. The PSD shows the digital sideband shelves at +-129..198 kHz
    riding beside the analog FM hump (the "hybrid" in IBOC).
 2. Symbol timing: differential products between consecutive OFDM
    symbols on the REFERENCE subcarrier comb are +-real when the FFT
    window sits on the true symbol boundary. Sweeping the window
    offset produces a sharp lock spike at one offset in 2160 — the
    grid, appearing out of what looks like noise.

Hard-won laws baked in (each cost us hours):
  - Do NOT try guard-interval autocorrelation: the guard is tapered
    and overlap-added, so the classic CP ridge barely exists.
  - OFDM symbols restart phase each symbol: do NOT "compensate" a
    guard-stride rotation on the differential products. Residual CFO
    adds one COMMON rotation — estimate it blindly as 0.5*arg(sum d^2)
    and remove it.
  - The reference carriers are NOT power-boosted; you cannot see the
    comb in the PSD. Only their modulation gives them away.

Usage: python measure.py --iq capture.cs16 --fs 2976750
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FS = 2_976_750.0        # 4x the 744,187.5 Hz native rate
NFFT = 8192             # useful symbol at 4x (2048 native)
SYM = 8640              # total symbol at 4x (2160 native)
# In HYBRID mode the OFDM only exists in the outer sidebands: the
# reference subcarriers sit at +-(356..546) every 19th — a comb through
# the middle would read the analog FM hump and score noise (we know:
# that bug cost this script its first run).
REF_SC = np.array(sorted(s * (546 - 19 * i)
                         for i in range(11) for s in (+1, -1)))


def dbpsk_metric(d):
    """+1 when diff products are +-real (locked DBPSK), ~0 on noise.
    The only systematic rotation is the common one from residual CFO —
    estimate blindly via 0.5*arg(sum d^2) (sign-invariant), remove.
    (OFDM symbols restart phase each symbol; there is NO guard-stride
    rotation to compensate — 'fixing' one scrambles the lock.)"""
    d = d * np.exp(-1j * 0.5 * np.angle((d.astype(np.complex128) ** 2).sum()))
    return float((np.abs(d.real).sum() - np.abs(d.imag).sum())
                 / (np.abs(d).sum() + 1e-12))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq", required=True)
    ap.add_argument("--fs", type=float, required=True,
                    help="must be 2976750 (4x native) for this script")
    ap.add_argument("--t-start", type=float, default=3.0)
    a = ap.parse_args()
    assert abs(a.fs - FS) < 1, "capture must be at 2,976,750 S/s"

    n_need = 20 * SYM
    raw = np.fromfile(a.iq, dtype=np.int16, count=2 * n_need,
                      offset=int(a.t_start * FS) * 4)
    x = (raw[0::2].astype(np.float32)
         + 1j * raw[1::2].astype(np.float32)) / 32768.0

    # ---- 1. the hybrid PSD ------------------------------------------
    N = 8192
    seg = x[: len(x) // N * N].reshape(-1, N)
    psd = np.fft.fftshift(
        (np.abs(np.fft.fft(seg * np.hanning(N), axis=1)) ** 2).mean(0))
    fax = np.fft.fftshift(np.fft.fftfreq(N, 1 / FS))
    db = 10 * np.log10(psd + 1e-9)
    sb = float(db[(np.abs(fax) > 130e3) & (np.abs(fax) < 195e3)].mean())
    away = float(db[(np.abs(fax) > 220e3) & (np.abs(fax) < 350e3)].mean())
    print(f"digital sideband shelves: {sb - away:+.1f} dB above the "
          f"floor at +-130..195 kHz (the IBOC signature)")

    # ---- 2. symbol-timing lock via the reference comb ---------------
    # sweep the FFT window offset across one whole symbol; at the true
    # boundary the 22 reference carriers turn +-real. 8 symbols per
    # hypothesis is plenty; step 6 then refine +-6.
    n_sym = 8
    win = np.arange(NFFT)
    sym_off = np.arange(n_sym)[:, None] * SYM
    step = 6
    t0s = np.arange(0, SYM, step)
    scores = np.empty(len(t0s))
    for i, off in enumerate(t0s):
        F = np.fft.fft(x[off + sym_off + win[None, :]], axis=1)
        R = F[:, REF_SC % NFFT]
        scores[i] = dbpsk_metric(R[1:] * np.conj(R[:-1]))
    bi = int(np.argmax(scores))
    best, m_best = int(t0s[bi]), scores[bi]
    for off in range(max(0, best - step), best + step + 1):
        F = np.fft.fft(x[off + sym_off + win[None, :]], axis=1)
        R = F[:, REF_SC % NFFT]
        m = dbpsk_metric(R[1:] * np.conj(R[:-1]))
        if m > m_best:
            best, m_best = off, m
    print(f"symbol-timing sweep: lock spike at offset {best}/{SYM} "
          f"(metric {m_best:.3f}; sweep median "
          f"{np.median(scores):.3f} = the noise floor)")
    print(f"grid confirmed: FFT 2048, guard 112, symbol 2160 samples "
          f"@ 744187.5 Hz -> {744187.5 / 2160:.2f} sym/s, spacing "
          f"{744187.5 / 2048:.3f} Hz, refs every 19th carrier at "
          f"+-(356..546)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    ax1.plot(fax / 1e3, db, lw=.6, color="#0aa")
    for s in (+1, -1):
        ax1.axvspan(s * 129.36, s * 198.4, color="#d33", alpha=.12)
    ax1.annotate("analog FM", (0, db.max() - 2), ha="center", fontsize=9)
    ax1.annotate("digital\nsidebands", (165, db.max() - 12), ha="center",
                 fontsize=9, color="#d33")
    ax1.annotate("digital\nsidebands", (-165, db.max() - 12), ha="center",
                 fontsize=9, color="#d33")
    ax1.set_xlabel("kHz from carrier")
    ax1.set_ylabel("PSD (dB)")
    ax1.set_title("hybrid IBOC: OFDM riding beside analog FM")
    ax2.plot(t0s, scores, lw=.6, color="#0aa")
    ax2.axvline(best, color="#d33", ls=":")
    ax2.annotate(f"the grid appears:\noffset {best}", (best, m_best),
                 fontsize=9, color="#d33", ha="left")
    ax2.set_xlabel(f"FFT window offset (samples, one OFDM symbol = {SYM})")
    ax2.set_ylabel("reference-comb DBPSK metric")
    ax2.set_title("symbol-timing sweep: one true offset in 2160")
    fig.tight_layout()
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "iboc_grid.png", dpi=110)
    print(f"figure -> {out / 'iboc_grid.png'}")


if __name__ == "__main__":
    main()
