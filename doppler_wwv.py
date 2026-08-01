#!/usr/bin/env python3
"""doppler_wwv.py - the ionosphere breathing, in milli-Hertz.

WWV's 10 MHz carrier is a rock-stable atomic reference. But the ionosphere it
bounces off is not stable - as reflecting layers rise and fall (sunrise/sunset),
as travelling ionospheric disturbances (atmospheric gravity waves) ripple
through, and as solar flares suddenly ionize the D-layer, the PATH LENGTH
changes, Doppler-shifting the received carrier by milli-Hertz. Track that shift
finely enough and you are watching near-space weather from a wire in the yard.
The dlayer-diary sequel: that measured absorption; this measures MOTION.

Method: FFT-locate the carrier, mix to ~DC, then estimate instantaneous
frequency from the phase slope of the complex carrier over 1 s blocks (a
carrier this strong gives sub-milli-Hz precision). Detrend the receiver's own
offset; what remains is ionospheric Doppler.

  python doppler_wwv.py captures/wwv10_night.cs16
  python doppler_wwv.py --selftest
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FS = 250_000.0  # rate-ok: offline analyzer of archival 250k WWV captures - no
#                 SDR open here; NEW captures for this tool must be made at
#                 >=2.048M (law 8/01) and decimated before analysis


def carrier_offset(iq, fs, search=2000.0):
    """Coarse carrier offset from the averaged spectrum (Hz from center)."""
    N = 1 << 16
    m = len(iq) // N * N
    seg = iq[:m].reshape(-1, N) * np.hanning(N).astype(np.float32)
    P = (np.abs(np.fft.fftshift(np.fft.fft(seg, axis=1), axes=1)) ** 2).mean(0)
    c = N // 2
    k = int(search / (fs / N))
    band = P[c - k:c + k]
    return (int(np.argmax(band)) - k) * fs / N


def doppler_series(iq, fs, off, block_s=1.0):
    """Instantaneous carrier frequency (Hz, relative to off) per block, from
    the phase slope of the mixed-to-DC carrier. Narrow LP first so only the
    carrier (not the WWV tones/voice) drives the phase."""
    from scipy.signal import firwin, lfilter
    n = np.arange(len(iq), dtype=np.float64)
    x = iq * np.exp(-2j * np.pi * off / fs * n)
    lp = firwin(1023, 20.0 / (fs / 2)).astype(np.float32)     # +/-20 Hz around carrier
    xf = lfilter(lp, 1.0, x)
    nb = int(block_s * fs)
    amp_all = np.abs(xf)
    strong = 0.5 * np.median(amp_all)            # reject only genuine fades (< half median)
    freqs, times = [], []
    for b in range(0, len(xf) - nb, nb):
        seg = xf[b:b + nb]
        a = np.abs(seg)
        # GATE: skip blocks with deep fades or AM-tick amplitude collapse
        if np.median(a) < strong or (a.min() / (a.max() + 1e-12)) < 0.25:
            freqs.append(np.nan); times.append(b / fs); continue
        ph = np.unwrap(np.angle(seg))
        t = np.arange(len(ph))
        c = np.polyfit(t, ph, 1)
        resid = np.std(ph - np.polyval(c, t))    # fading/modulation corruption
        if resid > 0.5:                          # rad - too corrupt to trust
            freqs.append(np.nan); times.append(b / fs); continue
        freqs.append(c[0] * fs / (2 * np.pi))
        times.append(b / fs)
    f = np.array(freqs)
    return np.array(times), f


def selftest():
    print("doppler_wwv selftest: synth 10 MHz carrier + 8 mHz ionospheric ripple")
    fs = 250_000.0
    dur = 60.0
    t = np.arange(int(fs * dur)) / fs
    off = 137.0
    # a slow 8 mHz sinusoidal Doppler (a gravity-wave-like ripple) + noise
    dopp = 0.008 * np.sin(2 * np.pi * t / 25.0)
    phase = 2 * np.pi * (off * t + np.cumsum(dopp) / fs)
    iq = (np.exp(1j * phase) + 0.05 * (np.random.randn(len(t)) + 1j * np.random.randn(len(t)))).astype(np.complex64)
    o = carrier_offset(iq, fs)
    _, f = doppler_series(iq, fs, o)
    f = f[np.isfinite(f)]
    f = f - np.median(f)
    amp = (f.max() - f.min()) / 2
    print(f"  recovered carrier offset {o:.1f} Hz (true 137)")
    print(f"  Doppler ripple amplitude {amp*1e3:.1f} mHz (injected 8 mHz)  "
          f"{'PASS' if 4e-3 < amp < 16e-3 else 'FAIL'}")
    return 0 if 4e-3 < amp < 16e-3 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("iq", nargs="?")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest or not a.iq:
        sys.exit(selftest())
    p = Path(a.iq)
    raw = np.fromfile(p, np.int16).astype(np.float32) / 32768.0
    iq = (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)
    off = carrier_offset(iq, FS)
    t, f = doppler_series(iq, FS, off)
    good = np.isfinite(f)
    fg = f[good] - np.median(f[good])            # remove receiver offset
    print(f"{p.name}: carrier {off:+.1f} Hz, {good.sum()}/{len(f)} clean 1-s blocks "
          f"({100*good.mean():.0f}% - rest faded/modulation-gated)")
    if good.sum() < 5:
        print("  too few clean blocks (weak/fading signal - 10 MHz is a night band)")
        return
    print(f"  Doppler std {np.std(fg)*1e3:.1f} mHz, peak-peak {(fg.max()-fg.min())*1e3:.1f} mHz")
    print(f"  {'STABLE ionosphere' if np.std(fg) < 0.05 else 'MOVING (layer motion / TID / disturbance)'}")
    f = f - np.nanmedian(f)
    out = HERE / "captures" / (p.stem + "_doppler.json")
    out.write_text(json.dumps({"carrier_off_hz": round(off, 2),
                               "t_s": [round(x, 1) for x in t],
                               "doppler_mhz": [round(x * 1e3, 2) for x in f]}, indent=1))
    print(f"  -> {out.name}")


if __name__ == "__main__":
    main()
