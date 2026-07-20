# 📡 The Radio Grid Atlas

**Every digital radio signal is built on a grid** — a known-in-advance
skeleton of sync words, pilots, symbol clocks, and subcarriers that the
transmitter promises and the receiver counts on. If you can *measure*
the grid, you can lock to it, equalize against it, and turn its error
into a live quality dial.

This atlas is a collection of grids we have **personally measured off
the air** with a hobby SDR (SDRplay RSPdx, consumer antennas, a
suburban rooftop). Every entry contains three things:

1. **`GRID.md`** — the grid's numbers (symbol rates, sync sequences,
   pilot positions, FEC structure) *plus the physics of why it's built
   that way* and the story of how we measured it.
2. **`measure.py`** — a self-contained Python script (numpy/scipy)
   that re-derives the key numbers from a raw IQ capture. Not a
   decoder — a *proof*. Run it on your own capture and watch the grid
   appear.
3. **`figures/`** — the plots those scripts drew from **our actual
   captures**: annotated spectra, correlation spikes, constellations.
   Nothing simulated, nothing copied from a datasheet.

## The method (works on any signal)

Everything here was found the same way:

1. **Find the known structure.** Standards docs, open-source decoders
   (a "referee"), or plain spectrum staring: something in the signal
   is transmitted the same way every time.
2. **Hypothesize → experiment → conclude.** Guess a parameter, write
   ten lines of numpy that would reveal it if true, run it on a
   capture, keep what survives. Wrong guesses are data too.
3. **Build the synthetic transmitter first.** Before trusting a
   measurement on a real capture, generate a perfect fake signal and
   confirm the measurement recovers what you put in. (A 90° phase
   convention bug will read as "no stereo" on real FM forever; the
   synthetic TX finds it in one run.)
4. **Cross-check against a referee.** When an open decoder exists,
   your derived numbers must match its behavior (our CFO estimate vs
   nrsc5's "Frequency offset:" line, to the tenth of a Hz).
5. **Turn the grid into a dial.** Once locked, the residual error
   against the known structure is a live signal-quality meter — the
   foundation for adaptive receivers.

## The grids (measured so far)

| entry | signal | the grid, in one line | status |
|---|---|---|---|
| [atsc-8vsb](atsc-8vsb/) | ATSC 1.0 digital TV | 10.762 Msym/s 8-VSB, 832-symbol segments, PN511 field sync every 313 segments, pilot tone | ✅ measured |
| [fm-stereo](fm-stereo/) | Broadcast FM composite | 19 kHz pilot, (L−R) DSB-SC on 38 kHz = 2× pilot, RDS on 57 kHz = 3× pilot | ✅ measured |
| [rds](rds/) | Radio Data System | 1187.5 bps BPSK on the 57 kHz subcarrier (= 19 kHz × 3; bitrate = 57 k/48) | ✅ measured |
| [nrsc5-fm-hybrid](nrsc5-fm-hybrid/) | HD Radio (FM IBOC) | OFDM: 2048-FFT @ 744.1875 kHz, 363.373 Hz spacing, tapered 112-sample guard, reference subcarriers every 19th carrying 32-bit DBPSK | ✅ measured |
| [rs41-radiosonde](rs41-radiosonde/) | Vaisala RS41 weather balloon | 4800 Bd GFSK, 64-bit whitened sync, LFSR whitening, RS(255,231) over GF(256)/0x11D | ✅ measured |
| [adsb-1090](adsb-1090/) | ADS-B / Mode S | 1090 MHz PPM, 8 µs preamble (pulses at 0, 1, 3.5, 4.5 µs), 112-bit extended squitter, CRC-24 | ✅ measured — 31 CRC-verified frames on rabbit ears |
| [meteor-lrpt](meteor-lrpt/) | Meteor-M LRPT weather satellite | 72 kSym/s QPSK, K=7 r=1/2 Viterbi, CCSDS frames behind 0x1ACFFC1D ASM, RS(255,223) | ✅ measured |
| [wwv-time](wwv-time/) | WWV / NIST time (10 MHz) | 1000 Hz ticks, 100 Hz BCD time code, 500/600 Hz minute tones — every element cesium-derived | ✅ measured (tone grid; per-second fold honest-inconclusive under evening QSB) |
| [ft8-hf](ft8-hf/) | FT8 (40 m amateur) | 15 s slots aligned to ABSOLUTE UTC, 79 × 0.16 s symbols, 8-GFSK with spacing = baud = 6.25 Hz | ✅ measured — the band folds at the wall clock, 6.8 dB |
| [wwvb-60khz](wwvb-60khz/) | WWVB / NIST time (60 kHz) | pulse-width AM at 1 bit/s: carrier drops 0.2/0.5/0.8 s for 0/1/marker, markers every 10 s | ✅ measured — three pulse islands, marker comb R = 1.000 |
| [flex-pager](flex-pager/) | FLEX paging (929 MHz) | 1.875 s frame comb (240 s / 128), 4-FSK at ±1600/±4800 Hz, 3200 sym/s | ✅ measured — went in hunting POCSAG, the band said FLEX |
| [gps-l1-ca](gps-l1-ca/) | GPS L1 C/A | 1023-chip Gold codes @ 1.023 Mchip/s (1 ms epochs), 50 bps nav bits, code/carrier locked 1:1540 | ✅ measured — 4 SVs pulled from 20 dB below the noise floor, on rabbit ears |
| [ais-vhf](ais-vhf/) | AIS marine (161.975/162.025 MHz) | 9600 bps GMSK, NRZI + HDLC 0x7E flags, SOTDMA 26.67 ms slot comb, CRC-16/X.25 | ✅ measured — 179 CRC-verified frames; slot comb aligned across 3 transmitters to 68 µs |
| [hdam-iboc](hdam-iboc/) | HD Radio (AM all-digital, 820 kHz) | OFDM: 181.6864 Hz spacing = 172.265625 sym/s × 270/256, BPSK refs at ±181.7 Hz, 1.48607 s L1 frames | ✅ measured — block sync 90%, psmi=2 says all-digital in the broadcast's own bits |

## Why publish this

Signal identification wikis tell you what a signal *looks like*. Open
decoders tell you what to *run*. This atlas is the missing middle: the
**measured numbers, the capture they came from, and the script that
proves them** — so anyone with an SDR can reproduce the measurement,
and anyone building a receiver can start from ground truth.

If you measure a grid we haven't, PRs are welcome — bring the numbers,
the script, and the figure.

## Family

Built by the tuna fleet: [Software-TV-Tuner](https://github.com/Felbs/Software-TV-Tuner)
(the ATSC grid put to work), [albacore](https://github.com/Felbs/albacore)
(instrumented HD Radio), [wxTuna](https://github.com/Felbs/wxTuna)
(LRPT + radiosondes), [aeroTuna](https://github.com/Felbs/aeroTuna)
(ADS-B), [gr-radiotuna](https://github.com/Felbs/gr-radiotuna)
(FM/RDS + the method's science deck).

*License: MIT (our text, scripts, and figures). The grids themselves
are physics — nobody owns them.*
