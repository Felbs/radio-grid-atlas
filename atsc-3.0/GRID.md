# ATSC 3.0 (NextGen TV) — the lighthouse in front of an unknown harbour

## The grid

| parameter | value | why |
|---|---|---|
| Bootstrap sample rate | **6.144 MHz**, always | fixed by A/321 no matter what the rest of the frame does — this is the whole trick |
| Bootstrap bandwidth | **4.5 MHz** | fits inside the 6, 7 and 8 MHz channel rasters used worldwide, so one detector works everywhere |
| FFT | **2048** → **3000 Hz** subcarrier spacing | 6.144 MHz / 2048 |
| Symbol | **3072 samples = 500 µs**, laid out **C(520) │ A(2048) │ B(504)** | C is a copy of the *last* 520 samples of A — a cyclic prefix by another name |
| Bootstrap | **4 symbols = 2.000 ms**, opening every frame | symbol 0 acquires; symbols 1–3 carry signalling as cyclic shifts |
| Sequence | Zadoff-Chu, **root q = 137**, length **1499** | 1499 is the largest prime that fits 4.5 MHz at 3 kHz spacing; primality is what makes every cyclic shift orthogonal to every other |
| Scrambling | 16-bit PN, p(x) = x¹⁶ + x¹⁵ + x¹⁴ + x + 1 | flattens the ZC's autocorrelation sidelobes so a shift can't be mistaken for its neighbour |

## The problem this grid exists to solve

ATSC 3.0 is deliberately, aggressively configurable. FFT size, guard
interval, pilot pattern, constellation, code rate, the number of PLPs and
whether they're layered — a broadcaster picks all of it, and can change it.
A receiver arriving cold knows none of it.

So the standard nails exactly one thing down. Every frame opens with a
bootstrap at a fixed 6.144 MHz and a fixed 4.5 MHz width, *whatever follows*.
It is a known lighthouse in front of an unknown harbour: find it, and it
tells you the sample rate of everything after it, how long until the next
frame, and which version of the standard you're looking at. Everything else
in the receiver is downstream of that one fixed handhold.

This is the opposite design choice from [ATSC 1.0](../atsc-8vsb/GRID.md),
whose grid is rigid all the way down — 8-VSB has exactly one symbol rate,
inherited from NTSC's sound carrier, and a receiver can assume it. The two
entries next to each other are a nice illustration of what you buy and what
you pay for with configurability.

## What we measured

Off our own antenna: an ATSC 3.0 UHF channel, yagi, RSPdx at 8 MS/s, 3 s of IQ.

```
Stage 2, Zadoff-Chu matched filter
  ZC root 137 (A/321 6.1)            peak/median   27.74
  CONTROL root 136                   peak/median    5.00
  CONTROL root 138                   peak/median    4.95
  -> root 137 beats its neighbours by 5.5x

  bootstraps in 3.00 s                12
  frame period                        247.111 ms
  spread                                0.000 ms over 11 gaps
  frame rate                            4.047 frames/s
```

Zero spread across eleven consecutive frame gaps. The transmitter's frame
clock and our receiver's resampled clock agree to below our measurement
resolution over three seconds.

## The trick worth stealing — and the measurement that nearly fooled us

You can find the bootstrap knowing nothing about Zadoff-Chu at all. Part C
is a copy of the last 520 samples of part A, so the waveform repeats itself
at a lag of exactly 2048. Correlate the capture against itself at that lag,
comb four symbols together, and the bootstrap stands up:

```
single lag 2048, no comb           peak/median   25.91
CONTROL wrong lag 2000             peak/median    4.94   <- fails, as it must
selectivity in lag                       3.82x   PASS

symbol peaks at   [0, 3078, 6150, 9222] samples
measured symbol period  3074.0 samples = 500.33 us   (grid says 3072 = 500)
error 2.0 samples                                     PASS
```

Two samples out of 3072 — 0.065 % — with no sequence, no table, no version
number. The four peaks 500 µs apart in the left-hand figure *are* the grid.

**The part that nearly fooled us.** Our first instinct for proving the
spacing was to sweep the comb and score it. That sweep peaks 180 samples
away from 3072 and is 720 samples wide at half height, and a control set 96
samples off scores 16.29 against the true 18.85 — barely a difference. Read
naively that says the detector is weak.

It doesn't. Part C is 520 samples long, so the correlation sits on a plateau
~520 wide, and a four-tap comb offset by 96 per tap has drifted only 288 by
the last tap — still on the plateau. **The control was inside the
measurement's own resolution, where nothing can fail.** The blunt instrument
was the *score*, not the detector: the peak *positions* were carrying the
answer to two samples the whole time.

The tempting fix is to widen the control until it fails and report a pass.
That is tuning a gate to produce the answer you wanted. The honest fix is to
measure the resolution, say what the detector can and cannot distinguish,
and then find a statistic that actually resolves it. `measure.py` prints all
three: the sweep, its resolution, and the peak positions.

**What geometry still can't do:** it finds a cyclic-prefixed OFDM structure
and measures its period, but it cannot tell you *which standard* you're
looking at — every OFDM payload symbol has a cyclic prefix too. Only the
Zadoff-Chu root identifies the signal as ATSC 3.0 major version 0. Geometry
finds it; the sequence names it.

## An open honesty note on the PN

A/321 Figure 5.2 leaves the shift-register wiring open to more than one
honest reading, so `measure.py` tries four and reports which locks:

```
PN wiring spec            27.74
PN wiring recip           27.74
PN wiring spec_revseed     4.51
PN wiring recip_revseed    4.51
```

The tap set turns out not to matter — `spec` and `recip` score identically —
while the **seed bit order** does. We have not resolved whether that
degeneracy is a property of this polynomial or a coincidence of the metric,
and we would rather say so than present one reading as settled.

## Reproduce it

```
python measure.py --iq your_capture.cs16 --fs 8e6
```

Interleaved int16 IQ, centred on any ATSC 3.0 channel, anything from about
6.5 MS/s up. The script resamples to 6.144 MHz itself, runs both stages with
their controls, prints the frame period, and draws `figures/`.

You do not need our capture — you need *a* capture. That's the point of a
grid that's fixed by the standard.
