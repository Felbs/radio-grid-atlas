# relativity.py — checkpoint 2026-07-25 morning (80% built)

WORKS: acquire/track PRN15 (C/N0 40.3), 5899 nav bits, preamble grid found
(18 subframes, mod-300 phase), IS-GPS-200 parity verified correct (61 clean
words), subframes 1/3/4 decoded w/ consistent TOWs 24408-24418 + WN. We ARE
reading real GPS nav data from our own capture.

BLOCKER: only ~25% of stream decodes cleanly — carrier cleanup insufficient.
Tried: (a) global stitch (p² block-freq + unwrap ref) = 26% words; (b)
transition decode + preamble anchor = 17% (WORSE — adjacent-bit coherence
broken); (c) per-subframe cubic polyfit on p² phase = 8% (worse still).
Model falsified: failure is NOT slow drift or π-seams — likely genuine phase
dynamics (possibly receiver clock/TCXO wander at seconds scale) needing a
REAL Costas PLL (per-ms closed loop, bandwidth ~10-20 Hz) instead of open-
loop batch correction. Subframe 2 (e, sqrtA — the physics fields) always
falls in bad stretches; never clean.

NEXT (fresh context): implement classic Costas on the 1ms prompts:
phase err = atan(Q/I) per prompt, 2nd-order loop (Bn~15Hz), run over full
118s → bits from I sign → reuse harvest (union). Physics code + wobble plot
+ RELATIVITY.md writer are DONE and waiting — they fire the moment e/sqrtA
decode. Also: capture t0 known (meta json) → could cross-check decoded
ephemeris against public broadcast ephemeris for that date (honesty check).

## Update (morning, cont.): root cause FOUND — it's uniform ~4% BER
Falsified: phase (3 fixes incl. Costas PLL: no change), code alignment (power
flat 98-101dB all 120s, resid 0.59 samp), grid slip (127us, negligible).
The tell: 26% word success == (0.96)^30 — a UNIFORM ~4% bit error rate.
THE FIX: ephemeris repeats identically every 30 s -> we hold 3-4 copies of
each subframe. Assign sfid cycle to the 18 grid starts (anchor from any clean
HOW; sfid advances +1 mod 5 per 300 bits), BIT-WISE MAJORITY VOTE across
same-sfid repeats, then parity-harvest the voted bits. 4% -> ~0.5% BER ->
~86% words -> subframe 2 (e, sqrtA) should decode. Physics stage then fires.
