# Hydrogen-line first light (2026-07-25)

**Setup:** RSPdx @ 1420.405751 MHz, frequency-switched 1420.406 / 1423.406 MHz,
bias-T on, **Antenna B** (Old Faithful, the bias-T port), 2-minute integration.

## What worked (the chain)
The observatory runs end-to-end on real hardware: the SDR tunes to the hydrogen
frequency, the bias-T powers, frequency switching and integration proceed, and
`hi_run.npz` is written with a full ratio spectrum. 200 dumps in 2 minutes. The
plumbing is validated on-sky, not just in the selftest.

## What it did NOT find (honestly)
**No real 21 cm line** — and the tool's naive `find_line` correctly *looked like*
it found one, which is the lesson. It reported "SNR ~10 at v = -83 km/s", but:
- the velocity **flipped between -83 and +84 km/s** across dumps — a real line
  sits at ONE fixed velocity; a sign-flipping pair does not;
- the dominant structure is a **dip at band center** (the frequency-switch
  self-residual), and the "peak" is a **narrow RSPdx spur (birdie)** near
  +/-393 kHz, not sky.

Both are receiver artifacts. The reason there's no real signal is expected:
**Antenna B is a UHF TV antenna, effectively deaf at 1420 MHz (L-band).** With no
sky signal, everything in the ratio is instrument. This is exactly the
"genuine long shot without a horn" the README warned about.

## The honest call on the finder
First light tempted an "improvement" to `find_line` (detrend + spur rejection).
It was reverted: it broke the selftest (over-subtracted the injected line) and
still flagged the birdie. **Tuning artifact-rejection against a deaf antenna
overfits noise.** The right sequence is: get a 1420-capable antenna (the ~$50
horn on the shopping list, or a dish) so there's a REAL line to calibrate
against, THEN build spur/baseline rejection using that ground truth.

## Verdict
Chain: **validated.** Detection: **pending the horn.** The observatory is ready;
the antenna is the missing piece, precisely as designed. Next real step is an
overnight integration once a 1420 antenna is on the bias-T port.
