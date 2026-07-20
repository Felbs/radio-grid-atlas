#!/usr/bin/env python3
"""Measure the VOR grid: 30 Hz AM + a 9960 Hz FM subcarrier (+-480 Hz at
30 Hz), whose phase difference IS the receiver's magnetic bearing from
the station.

The two 30 Hz signals travel very different roads — one as plain AM, one
as FM on a subcarrier — precisely so their phase difference survives the
receiver chain. This script measures every grid element coherently, which
matters because a rooftop-antenna VOR carrier can sit ~19 dB below the
noise in its own channel: envelope detection dies (noise x noise), but a
carrier-locked coherent chain with 60 s of integration recovers the whole
grid. Run --selftest to see the method prove itself on a synthetic
transmitter with a known bearing before you trust it on air.

PRIVACY: the measured absolute bearing locates the receiver (two stations
= a position fix), so by default only the bearing STABILITY is printed.
Pass --reveal-bearing if you want the number on your own screen.

Usage:  python measure.py --iq your_capture.cs16 --fs 250000
        python measure.py --selftest
Output: printed grid measurements + figures/ plots (only with --figures).
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.signal import firwin, fftconvolve

HERE = Path(__file__).resolve().parent

MORSE = {
    ".-": "A", "-...": "B", "-.-.": "C", "-..": "D", ".": "E", "..-.": "F",
    "--.": "G", "....": "H", "..": "I", ".---": "J", "-.-": "K", ".-..": "L",
    "--": "M", "-.": "N", "---": "O", ".--.": "P", "--.-": "Q", ".-.": "R",
    "...": "S", "-": "T", "..-": "U", "...-": "V", ".--": "W", "-..-": "X",
    "-.--": "Y", "--..": "Z",
}


# ---------- primitives ----------------------------------------------------

def find_carrier(x, fs, search_hz=2000.0):
    """Strongest line within +-search_hz of the dial (bounded by the 50 kHz
    VOR channel plan + SDR ppm physics — not by the answer). Windowed, or
    leakage from out-of-band birdies steals the pick. Returns f0, C/N0."""
    n = min(1 << 20, len(x))
    seg = x[: n * (len(x) // n)].reshape(-1, n)
    S = (np.abs(np.fft.fft(seg * np.hanning(n), axis=1)) ** 2).mean(axis=0)
    f = np.fft.fftfreq(n, 1 / fs)
    m = np.abs(f) < search_hz
    f0 = f[m][np.argmax(S[m])]
    binw = fs / n
    C = S[np.abs(f - f0) < 1.5].sum()
    N0 = np.median(S[(np.abs(f - f0) > 3000) & (np.abs(f - f0) < 8000)]) / binw
    return f0, 10 * np.log10(C / N0)


def carrier_lock(x, fs, f0, lock_bw, out_bw=15000.0, dec=8):
    """Derotate to baseband, then divide out the narrowband carrier phase.
    The carrier lands on the positive real axis; Re() is coherent AM."""
    n = np.arange(len(x))
    y = x * np.exp(-2j * np.pi * f0 * n / fs)
    y = fftconvolve(y, firwin(801, out_bw, fs=fs), mode="same")[::dec]
    fsd = fs / dec
    c = fftconvolve(y, firwin(20001, lock_bw, fs=fsd), mode="same")
    u = c / (np.abs(c) + 1e-12)
    return y * np.conj(u), fsd


def precise_tone(sig, fs, f_lo, f_hi):
    """Line frequency: coarse FFT peak, then coherent fine scan. (Parabolic
    interpolation alone left a 3 mHz bias = 73 deg of phase over 60 s.)"""
    S = np.abs(np.fft.rfft(sig - sig.mean()))
    fax = np.fft.rfftfreq(len(sig), 1 / fs)
    m = (fax >= f_lo) & (fax <= f_hi)
    f0 = fax[m][np.argmax(S[m])]
    snr = 20 * np.log10(S[m].max() / np.median(S[m]))
    t = np.arange(len(sig)) / fs
    d = sig - sig.mean()
    span = 2.0 * fs / len(sig)
    for _ in range(4):
        cand = f0 + np.linspace(-span, span, 21)
        f0 = cand[int(np.argmax([np.abs(np.sum(d * np.exp(-2j * np.pi * f * t)))
                                 for f in cand]))]
        span /= 5.0
    return f0, snr


def tone_amp(sig, fs, f):
    t = np.arange(len(sig)) / fs
    return 2 * np.abs(np.sum((sig - sig.mean()) * np.exp(-2j * np.pi * f * t))) / len(sig)


def band_complex(sig, fs, fc, cutoff, taps=1001):
    z = sig * np.exp(-2j * np.pi * fc * np.arange(len(sig)) / fs)
    return fftconvolve(z, firwin(taps, cutoff, fs=fs), mode="same")


def circ_stats(deg):
    v = np.exp(1j * np.radians(deg))
    mean = np.degrees(np.angle(v.mean())) % 360
    dev = np.degrees(np.angle(v * np.conj(v.mean() / np.abs(v.mean()))))
    return mean, float(np.sqrt((dev ** 2).mean()))


# ---------- the reference subcarrier: matched filter ----------------------

def matched_subcarrier(s, fs, f30, beta=16.0):
    """Find the FM subcarrier by correlating against its known waveform
    class e^{j beta sin(2 pi f30 t + psi)} — optimal at low SNR. f_sub is
    scanned over the published +-1% tolerance (9860..10060), psi over all
    360 deg; nothing narrower than physics."""
    t = np.arange(len(s)) / fs
    fc = 9960.0
    z = band_complex(s, fs, fc, 900.0)
    nfft = 1 << int(np.ceil(np.log2(len(s))))
    fax = np.fft.fftfreq(nfft, 1 / fs)
    sel = np.abs(fax) <= 100.0
    best = (-1.0, None, None)
    for psi in np.arange(0.0, 360.0, 4.0):
        q = z * np.exp(-1j * beta * np.sin(2 * np.pi * f30 * t + np.radians(psi)))
        Q = np.abs(np.fft.fft(q, nfft))
        i = np.argmax(Q[sel])
        if Q[sel][i] > best[0]:
            best = (Q[sel][i], fc + fax[sel][i], psi)
    return best[1], best[2]


def refine_fit(s, fs, f30, fsub, psi0, beta0=16.0):
    """Local max of |corr| over (psi, beta, delta-f). Returns psi, beta,
    fsub, amplitude."""
    t = np.arange(len(s)) / fs
    z = band_complex(s, fs, fsub, 900.0)

    def corr(psi, beta, df):
        m = np.exp(1j * (2 * np.pi * df * t
                         + beta * np.sin(2 * np.pi * f30 * t + np.radians(psi))))
        return np.abs(np.sum(z * np.conj(m)))

    psi, beta, df = psi0, beta0, 0.0
    for span in (8.0, 2.0, 0.5):
        g = psi + np.linspace(-span, span, 17)
        psi = g[int(np.argmax([corr(p, beta, df) for p in g]))]
        g = beta * (1 + np.linspace(-0.12, 0.12, 13) * span / 8)
        beta = g[int(np.argmax([corr(psi, b, df) for b in g]))]
        g = df + np.linspace(-0.05, 0.05, 11) * span / 8
        df = g[int(np.argmax([corr(psi, beta, d) for d in g]))]
    m = np.exp(1j * (2 * np.pi * df * t
                     + beta * np.sin(2 * np.pi * f30 * t + np.radians(psi))))
    amp = 2 * np.abs(np.sum(z * np.conj(m))) / len(z)
    return psi % 360, beta, fsub + df, amp


def disambiguate(s, fs, f30, fsub0, psi0):
    """A matched filter one Bessel line off still correlates (sum J_k J_k+n);
    the TRUE center maximizes the fitted amplitude. Try the neighbors."""
    fits = [refine_fit(s, fs, f30, fsub0 + k * f30, psi0) for k in (-2, -1, 0, 1, 2)]
    return max(fits, key=lambda r: r[3])


# ---------- bearing --------------------------------------------------------

def block_phases(s, fs, f30, fsub, beta, block_s):
    """Per-block phase of variable (30 Hz AM) and reference (subcarrier FM).
    Same f30 for both, so its small error cancels in the difference.
    bearing = psi_ref - phi_var  (convention PROVEN by --selftest)."""
    t = np.arange(len(s)) / fs
    z = band_complex(s, fs, fsub, 900.0)
    nb = int(block_s * fs)
    brg = []
    for k in range(len(s) // nb):
        sl = slice(k * nb, (k + 1) * nb)
        c = np.sum((s[sl] - s[sl].mean()) * np.exp(-2j * np.pi * f30 * t[sl]))
        phi_var = np.degrees(np.angle(c))
        g = np.arange(0, 360, 2.0)
        v = [np.abs(np.sum(z[sl] * np.conj(np.exp(
            1j * beta * np.sin(2 * np.pi * f30 * t[sl] + np.radians(p)))))) for p in g]
        p0 = g[int(np.argmax(v))]
        g = p0 + np.linspace(-2, 2, 21)
        v = [np.abs(np.sum(z[sl] * np.conj(np.exp(
            1j * beta * np.sin(2 * np.pi * f30 * t[sl] + np.radians(p)))))) for p in g]
        psi_ref = g[int(np.argmax(v))]
        brg.append((psi_ref - phi_var) % 360)
    return np.array(brg)


# ---------- main measurement ----------------------------------------------

def measure(x, fs, reveal=False, block_s=None, figures=False, tag="vor"):
    f0, cn0 = find_carrier(x, fs)
    print(f"carrier offset:   {f0:+9.1f} Hz from dial   C/N0 = {cn0:5.1f} dB-Hz")
    if cn0 < 12:
        print("carrier too weak to continue (C/N0 < 12 dB-Hz) - honest negative")
        return
    s = None
    yc, fsd = carrier_lock(x, fs, f0, lock_bw=3.0 if cn0 >= 20 else 0.8)
    s = yc.real
    dc = s.mean()

    f30, snr30 = precise_tone(s, fsd, 25.0, 35.0)
    am = tone_amp(s, fsd, f30) / dc
    print(f"variable 30 Hz:   {f30:9.4f} Hz ({snr30:.0f} dB)   AM depth {am:.3f} "
          f"(published 0.30)")
    if snr30 < 10:
        print("no 30 Hz AM - not a VOR (or below this method's floor)")
        return

    fsub0, psi0 = matched_subcarrier(s, fsd, f30)
    psi, beta, fsub, amp = disambiguate(s, fsd, f30, fsub0, psi0)
    print(f"subcarrier:       {fsub:9.2f} Hz   depth {amp / dc:.3f} (published 9960, 0.30)")
    print(f"  ratio fsub/f30: {fsub / f30:9.3f}      (published 332 - nested clocks)")
    print(f"  FM deviation:   {beta * f30:9.1f} Hz   (published +-480, index 16)")

    if block_s is None:
        block_s = 10.0 if cn0 >= 20 else 15.0
    brg = block_phases(s, fsd, f30, fsub, beta, block_s)
    mean, rms = circ_stats(brg)
    print(f"bearing:          recovered; stability {rms:.1f} deg rms over "
          f"{len(brg)} x {block_s:.0f} s blocks")
    if reveal:
        print(f"  ABSOLUTE BEARING (do not publish): {mean:.1f} deg magnetic FROM station")
    else:
        print("  absolute value withheld - it would locate the receiver "
              "(--reveal-bearing to print)")

    if figures:
        make_figures(s, fsd, f30, fsub, beta, brg, block_s, tag)
    return dict(f0=f0, cn0=cn0, f30=f30, am=am, fsub=fsub, beta=beta,
                sub=amp / dc, brg=brg, mean=mean, rms=rms)


def fold_waveforms(s, fs, f30, fsub, nbins=120):
    """Privacy-safe folded waveforms, each with its OWN arbitrary origin:
    (a) the 30 Hz AM cycle, (b) the subcarrier instantaneous frequency."""
    t = np.arange(len(s)) / fs
    ph = (t * f30) % 1.0
    idx = (ph * nbins).astype(int) % nbins
    lo = fftconvolve(s - s.mean(), firwin(1001, 60.0, fs=fs), mode="same")
    am_fold = np.bincount(idx, lo, nbins) / np.bincount(idx, None, nbins)
    z = band_complex(s, fs, fsub, 900.0)
    w = np.zeros(nbins, complex)
    np.add.at(w, idx, z)
    w /= np.bincount(idx, None, nbins)
    phase = np.unwrap(np.angle(w))
    finst = np.gradient(phase) * nbins * f30 / (2 * np.pi) + fsub
    # display smoothing only (circular, 5 bins) — measurements happen elsewhere
    ker = np.ones(5) / 5
    finst = np.real(np.fft.ifft(np.fft.fft(finst) * np.fft.fft(ker, len(finst))))
    am_fold = np.roll(am_fold, -int(np.argmax(am_fold)))
    finst = np.roll(finst, -int(np.argmax(finst)))
    return am_fold, finst


def make_figures(s, fs, f30, fsub, beta, brg, block_s, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = HERE / "figures"
    out.mkdir(exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8.6))

    S = np.abs(np.fft.rfft((s - s.mean()) * np.hanning(len(s)))) ** 2
    fax = np.fft.rfftfreq(len(s), 1 / fs)
    m = fax < 11000
    Sm, fm = S[m], fax[m]
    k = 20  # max-pool so narrow lines survive the plot decimation
    nk = len(Sm) // k
    Sp = Sm[: nk * k].reshape(nk, k).max(axis=1)
    fp = fm[: nk * k].reshape(nk, k).mean(axis=1)
    a = ax[0, 0]
    a.plot(fp / 1e3, 10 * np.log10(Sp + 1e-12), lw=.35, color="#0aa")
    a.set_title("coherent AM spectrum (carrier locked out)")
    a.set_xlabel("kHz")
    a.set_ylabel("dB")
    a.annotate("30 Hz variable", (0.03, .96), xycoords="axes fraction", fontsize=8)
    a.axvline(9.96, color="#c33", lw=.6, alpha=.6)
    a.annotate("9960 subcarrier", (9.96, .9), xycoords=("data", "axes fraction"),
               fontsize=8, ha="right", color="#c33")

    am_fold, finst = fold_waveforms(s, fs, f30, fsub)
    cyc = np.linspace(0, 1, len(am_fold), endpoint=False)
    a = ax[0, 1]
    a.plot(cyc, am_fold / np.abs(am_fold).max(), color="#0aa")
    a.set_title("variable signal: one 30 Hz AM cycle, folded 60 s\n(origin arbitrary)")
    a.set_xlabel("fraction of 30 Hz cycle")

    a = ax[1, 0]
    a.plot(cyc, finst, color="#c33")
    for y, lab in ((9480, "9960-480"), (9960, "9960"), (10440, "9960+480")):
        a.axhline(y, color="#888", lw=.5)
        a.annotate(lab, (1.0, y), fontsize=7, va="bottom", ha="right")
    a.set_title("reference signal: subcarrier instantaneous frequency,\n"
                "folded 60 s (origin arbitrary, display-smoothed)")
    a.set_xlabel("fraction of 30 Hz cycle")
    a.set_ylabel("Hz")

    a = ax[1, 1]
    mean, rms = circ_stats(brg)
    dev = (brg - mean + 180) % 360 - 180
    tb = (np.arange(len(brg)) + 0.5) * block_s
    a.plot(tb, dev, "o-", color="#0aa")
    a.axhline(0, color="#888", lw=.5)
    a.set_ylim(-max(10, 1.2 * np.abs(dev).max()), max(10, 1.2 * np.abs(dev).max()))
    a.set_title(f"bearing minus its mean: {rms:.1f} deg rms\n"
                "(absolute value withheld — it locates the receiver)")
    a.set_xlabel("s")
    a.set_ylabel("deg")
    fig.tight_layout()
    fig.savefig(out / f"{tag}.png", dpi=110)
    print(f"figure -> {out / f'{tag}.png'}")


# ---------- synthetic proof ------------------------------------------------

def selftest():
    """SYNTHETIC TX FIRST: known bearing 237 deg, C/N0 27 dB-Hz (as weak as
    our weakest usable capture). The pipeline must recover it."""
    fs, secs = 250_000.0, 60.0
    t = np.arange(int(fs * secs)) / fs
    rng = np.random.default_rng(11)
    az, phi0, f30t, fsubt = 237.0, np.radians(51.0), 30.004, 9962.0
    env = (1.0
           + 0.30 * np.cos(2 * np.pi * f30t * t - np.radians(az) + phi0)
           + 0.30 * np.cos(2 * np.pi * fsubt * t
                           + 16.0 * np.sin(2 * np.pi * f30t * t + phi0)))
    x = env * np.exp(2j * np.pi * -479.2 * t)
    x += np.sqrt(10 ** -2.7 * fs / 2) * (rng.standard_normal(len(t))
                                         + 1j * rng.standard_normal(len(t)))
    r = measure(x, fs, reveal=True)
    err = abs((r["mean"] - az + 180) % 360 - 180)
    assert err < 2.0, f"bearing error {err:.2f} deg"
    assert abs(r["f30"] - f30t) < 0.001 and abs(r["fsub"] - fsubt) < 0.5
    print(f"SELFTEST PASS: true bearing {az}, recovered {r['mean']:.2f} "
          f"(+-{r['rms']:.2f} rms) at C/N0 27 dB-Hz")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iq")
    ap.add_argument("--fs", type=float, default=250000.0)
    ap.add_argument("--secs", type=float, default=60.0)
    ap.add_argument("--block", type=float, default=None,
                    help="bearing block length, s (default: auto by C/N0)")
    ap.add_argument("--reveal-bearing", action="store_true")
    ap.add_argument("--figures", action="store_true")
    ap.add_argument("--tag", default="vor", help="figure filename stem")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    if not a.iq:
        ap.error("--iq required (or --selftest)")
    raw = np.fromfile(a.iq, dtype=np.int16, count=int(2 * a.secs * a.fs))
    x = raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)
    measure(x, a.fs, reveal=a.reveal_bearing, block_s=a.block,
            figures=a.figures, tag=a.tag)


if __name__ == "__main__":
    main()
