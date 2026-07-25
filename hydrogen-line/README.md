# 21 cm hydrogen-line observatory

Measure neutral hydrogen across the Milky Way at **1420.405751 MHz** — the same
signal professional radio astronomers use to map the galaxy. With enough
integration you can see the galactic rotation curve (the evidence for dark
matter) and our own motion through the galaxy.

## How it works
- **Frequency switching** — take a spectrum ON the line and one OFF it, divide.
  The receiver's own bandpass + standing-wave ripple (which dwarf the faint line)
  are identical in both and cancel, leaving only the sky.
- **Integration** — the line SNR grows as √time, so this is built to run for
  hours. It's the ideal leave-it-overnight experiment.

## Run it
```bash
python hi_line.py selftest                 # validate the method (no hardware)
python hi_line.py observe --minutes 480    # ~8 h overnight run
python hi_line.py show hi_run.npz          # re-read a finished run
```
`observe` saves `hi_run.npz` every ~20 dumps, so a morning Ctrl-C (or a crash)
still keeps everything integrated so far. Point at the Milky Way — the default
target is galactic (l=120°, b=0°), on the plane.

## Hardware honesty
- The **bias-T LNA at 1420 MHz** is essential and you have it — `observe` powers
  it automatically.
- A **directional antenna is the real limit.** HI is faint and broad; a horn or
  dish gives the gain to see it. With a wideband/omni antenna this is a genuine
  long shot — but because detection is integration-limited, a long overnight run
  is exactly how you'd give it the best chance. Worst case is an honest null that
  says "the horn is the missing piece"; best case is a faint line at the plane.

## Selftest result
Frequency switching lifts a synthetic line from SNR 3.6 (buried under bandpass
ripple) to 9.4, recovering its velocity within ~2 km/s of truth. That validates
the DSP; the sky is the remaining variable.

*Privacy: observer coordinates only refine a ~0.5 km/s term and default to 0 —
they are never hardcoded or committed. Targets are given in galactic (l, b).*
