#!/usr/bin/env python3
"""hi_line.py - a backyard 21 cm neutral-hydrogen spectrometer.

Neutral hydrogen across the galaxy emits at 1420.405751 MHz (the spin-flip "21 cm
line"). Point an antenna at the Milky Way, integrate long enough, and the line
rises out of the noise; its Doppler shift is the velocity of the gas, which is
how amateurs measure the galactic rotation curve (the evidence for dark matter)
and our own motion through the galaxy.

The enabling trick for a modest antenna is FREQUENCY SWITCHING: take a spectrum
with the tuner ON the line and another OFF it, then divide. The receiver's own
bandpass shape (which dwarfs the line) is identical in both and cancels, leaving
the sky. Everything else is INTEGRATION: the line SNR grows as sqrt(time), so
this is built to run for hours - the perfect leave-it-overnight experiment.

  python hi_line.py selftest                 # synthetic HI line, no hardware
  python hi_line.py observe --minutes 480    # overnight run (needs SDR + 1420 LNA)
  python hi_line.py show run.npz             # re-plot a finished run

Privacy: observer latitude/longitude only refine the small Earth-rotation term
and default to 0 - never hardcode or commit real coordinates (satellite/site
geometry can dox). Target is given in galactic (l, b); results stay local.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
F_HI = 1420.405751e6            # neutral-hydrogen rest frequency, Hz
C_KMS = 299_792.458            # speed of light, km/s


# ── core spectral estimation ────────────────────────────────────────────
def welch_psd(iq, fs, nfft=4096):
    """Averaged power spectrum (Welch, Hann, 50% overlap). Returns (freqs, psd)
    with freqs relative to the tuner center."""
    win = np.hanning(nfft)
    step = nfft // 2
    n = (len(iq) - nfft) // step + 1
    if n < 1:
        raise ValueError("segment shorter than one FFT")
    acc = np.zeros(nfft)
    for i in range(n):
        seg = iq[i * step:i * step + nfft] * win
        acc += np.abs(np.fft.fftshift(np.fft.fft(seg))) ** 2
    psd = acc / n
    freqs = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs))
    return freqs, psd


def freq_switch(psd_on, f_on, psd_off, f_off, freqs):
    """Frequency-switched spectrum: divide ON by OFF (bandpass cancels), then
    place on an absolute-frequency axis using the two tuner centers. Returns
    (abs_freqs_hz, ratio) on the ON grid."""
    ratio = psd_on / (psd_off + 1e-30)
    return f_on + freqs, ratio


# ── the line: find it, convert to velocity ──────────────────────────────
def doppler_velocity(f_obs, f_rest=F_HI):
    """Radio-convention radial velocity (km/s). Positive = receding."""
    return C_KMS * (f_rest - f_obs) / f_rest


def lsr_correction(l_deg, b_deg):
    """Projection of the Sun's peculiar motion onto the sightline, km/s - the
    'our own motion' term that must be removed to get the gas velocity relative
    to the Local Standard of Rest. Uses the standard solar motion (Schoenrich
    2010: U,V,W = 11.1, 12.24, 7.25 km/s in galactic X,Y,Z). Earth's orbital
    (+/-30 km/s annual) and rotation (~0.5 km/s) terms need the observation date
    and are added by observe() when astropy is available; this is the dominant,
    date-independent part."""
    l, b = np.radians(l_deg), np.radians(b_deg)
    U, V, W = 11.1, 12.24, 7.25
    # unit vector to target in galactic Cartesian, dotted with solar motion
    return U * np.cos(l) * np.cos(b) + V * np.sin(l) * np.cos(b) + W * np.sin(b)


def find_line(abs_freqs, spectrum, search_khz=500.0):
    """Locate the HI line near 1420.405 MHz as the largest deviation from 1.0
    (the frequency-switched baseline) within +/- search window. Returns
    (f_peak_hz, amplitude, snr)."""
    m = np.abs(abs_freqs - F_HI) < search_khz * 1e3
    if not m.any():
        return None, 0.0, 0.0
    band = spectrum[m] - 1.0
    fb = abs_freqs[m]
    k = int(np.argmax(band))
    # SNR vs the RMS of the off-line part of the window
    off = np.abs(fb - fb[k]) > 60e3
    noise = np.std(band[off]) if off.any() else np.std(band)
    return float(fb[k]), float(band[k]), float(band[k] / (noise + 1e-9))


# ── synthetic sky for the self-test ─────────────────────────────────────
def _synthetic_capture(rng, fs, nsamp, center, line_amp, line_v_kms=40.0):
    """IQ that mimics a receiver: a broad bandpass hump PLUS standing-wave RIPPLE
    (the real enemy - sharp features a raw line-finder mistakes for a line) + a
    faint HI line at its Doppler-shifted frequency + noise. The hump and ripple
    are functions of RELATIVE frequency, so they are identical ON and OFF and
    cancel under frequency switching - exactly why the method works. line_amp is
    the line's fraction (small - buried in a single dump, only integration and
    bandpass removal reveal it)."""
    freqs = np.fft.fftshift(np.fft.fftfreq(nsamp, 1 / fs))
    hump = 2.0 * np.exp(-((freqs - 0.15e6) / (fs * 0.35)) ** 2)
    ripple = 0.4 * np.cos(2 * np.pi * freqs / 0.25e6)     # standing-wave ripple
    shape = 1.0 + hump + ripple
    # HI line: a narrow gaussian, Doppler-shifted by line_v_kms, only if in band
    f_line = F_HI * (1 - line_v_kms / C_KMS) - center
    if line_amp > 0 and abs(f_line) < fs / 2:
        shape += line_amp * np.exp(-((freqs - f_line) / 30e3) ** 2)
    shape = np.clip(shape, 0.05, None)
    # build noise with that power spectral shape
    spec = (rng.standard_normal(nsamp) + 1j * rng.standard_normal(nsamp)) * np.sqrt(shape)
    iq = np.fft.ifft(np.fft.ifftshift(spec))
    return iq.astype(np.complex64)


def selftest():
    print("=" * 68)
    print("HYDROGEN-LINE SPECTROMETER selftest")
    print("  goal: pull a faint HI line out from under a huge receiver bandpass")
    print("  via frequency switching + integration, and read its velocity.")
    print("=" * 68)
    rng = np.random.default_rng(21)
    fs = 2.4e6
    nfft = 4096
    cap = nfft * 40          # ~per-dump samples
    true_v = 40.0            # km/s the injected line is Doppler-shifted by
    line_amp = 0.5           # buried in one dump (< ripple 0.4); integration wins

    # ON captures centered on the line; OFF captures 1.2 MHz away.
    f_on = F_HI
    f_off = F_HI + 1.2e6
    acc_on = np.zeros(nfft)
    acc_off = np.zeros(nfft)
    n_dumps = 60             # simulate integration over many dumps
    for _ in range(n_dumps):
        on = _synthetic_capture(rng, fs, cap, f_on, line_amp, true_v)
        off = _synthetic_capture(rng, fs, cap, f_off, 0.0, true_v)  # no line off-window
        _, p_on = welch_psd(on, fs, nfft)
        fr, p_off = welch_psd(off, fs, nfft)
        acc_on += p_on
        acc_off += p_off

    # a) WITHOUT frequency switching: can we see the line in raw ON spectrum?
    raw = acc_on / acc_on.max()
    raw_abs = f_on + fr
    _, _, raw_snr = find_line(raw_abs, raw / np.median(raw))

    # b) WITH frequency switching: divide out the bandpass
    abs_f, ratio = freq_switch(acc_on, f_on, acc_off, f_off, fr)
    f_peak, amp, snr = find_line(abs_f, ratio)
    v = doppler_velocity(f_peak)

    print(f"  raw ON spectrum:       line SNR {raw_snr:5.1f}  "
          f"({'buried under bandpass' if raw_snr < 5 else 'visible'})")
    print(f"  frequency-switched:    line SNR {snr:5.1f}  at {f_peak/1e6:.4f} MHz")
    print(f"  recovered velocity:    {v:+.1f} km/s   (injected {true_v:+.1f})")
    lsr = lsr_correction(120.0, 0.0)
    print(f"  LSR solar-motion term @ (l=120,b=0): {lsr:+.1f} km/s "
          f"(the 'our own motion' correction)")
    ok = snr > 6 and abs(v - true_v) < 8 and snr > raw_snr
    print("=" * 68)
    print(f"RESULT: freq-switching lifted the line from SNR {raw_snr:.1f} -> {snr:.1f}, "
          f"velocity within {abs(v-true_v):.1f} km/s. "
          f"{'PASS' if ok else 'NEEDS WORK'}")
    return 0 if ok else 1


# ── live overnight run ───────────────────────────────────────────────────
def observe(minutes, l_deg, b_deg, fs=2.4e6, antenna="Antenna B"):
    """Frequency-switched integration for `minutes`, saving cumulative spectra.
    Uses the bias-T LNA at 1420 MHz. Writes hi_run.npz periodically so a crash or
    a morning Ctrl-C still leaves the accumulated result."""
    import time
    sys.path.insert(0, r"Z:\src\hamTuna\tools")
    import cw
    from SoapySDR import SOAPY_SDR_RX
    nfft = 8192
    f_on, f_off = F_HI, F_HI + 3.0e6
    sdr, st = cw._open_sdr(antenna, fs)
    for key in ("biasT_ctrl", "biasT", "bias_tee"):
        try:
            sdr.writeSetting(key, "true")     # power the 1420 MHz LNA
        except Exception:
            pass
    acc_on = np.zeros(nfft)
    acc_off = np.zeros(nfft)
    fr = None
    out = HERE / "hi_run.npz"
    t_end = None
    dumps = 0
    print(f"[hi] integrating {minutes} min on {antenna}, bias-T on, "
          f"freq-switch {f_on/1e6:.3f}/{f_off/1e6:.3f} MHz")
    try:
        import datetime
        start = datetime.datetime.now()
        while (datetime.datetime.now() - start).total_seconds() < minutes * 60:
            for f, acc in ((f_on, acc_on), (f_off, acc_off)):
                sdr.setFrequency(SOAPY_SDR_RX, 0, f)
                time.sleep(0.05)
                iq = cw._grab(sdr, st, 0.3, fs)
                fr, p = welch_psd(iq, fs, nfft)
                acc += p
            dumps += 1
            if dumps % 20 == 0:
                abs_f, ratio = freq_switch(acc_on, f_on, acc_off, f_off, fr)
                fp, amp, snr = find_line(abs_f, ratio)
                v = doppler_velocity(fp) if fp else 0.0
                print(f"[hi] {dumps} dumps  line SNR {snr:4.1f}  v {v:+.0f} km/s", flush=True)
                np.savez(out, abs_f=abs_f, ratio=ratio, l=l_deg, b=b_deg, dumps=dumps)
    finally:
        try:
            sdr.deactivateStream(st); sdr.closeStream(st)
        except Exception:
            pass
        if fr is not None:
            abs_f, ratio = freq_switch(acc_on, f_on, acc_off, f_off, fr)
            np.savez(out, abs_f=abs_f, ratio=ratio, l=l_deg, b=b_deg, dumps=dumps)
            fp, amp, snr = find_line(abs_f, ratio)
            v = doppler_velocity(fp) if fp else 0.0
            vlsr = v + lsr_correction(l_deg, b_deg)
            print(f"[hi] DONE {dumps} dumps -> {out.name}")
            print(f"[hi] line SNR {snr:.1f} at velocity {v:+.1f} km/s "
                  f"(LSR-corrected {vlsr:+.1f} km/s)")


def show(path):
    d = np.load(path)
    abs_f, ratio = d["abs_f"], d["ratio"]
    fp, amp, snr = find_line(abs_f, ratio)
    v = doppler_velocity(fp) if fp else 0.0
    print(f"{Path(path).name}: {int(d['dumps'])} dumps, line SNR {snr:.1f} "
          f"at {fp/1e6:.4f} MHz, v {v:+.1f} km/s "
          f"(LSR {v + lsr_correction(float(d['l']), float(d['b'])):+.1f})")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("selftest")
    o = sub.add_parser("observe")
    o.add_argument("--minutes", type=float, default=480)
    o.add_argument("--l", type=float, default=120.0, help="target galactic longitude")
    o.add_argument("--b", type=float, default=0.0, help="target galactic latitude")
    o.add_argument("--antenna", default="Antenna B")
    s = sub.add_parser("show")
    s.add_argument("path")
    a = ap.parse_args()
    if a.cmd == "observe":
        observe(a.minutes, a.l, a.b, antenna=a.antenna)
    elif a.cmd == "show":
        show(a.path)
    else:
        sys.exit(selftest())


if __name__ == "__main__":
    main()
