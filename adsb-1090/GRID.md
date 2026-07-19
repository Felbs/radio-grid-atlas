# ADS-B / Mode S (1090 MHz) — the fastest grid most people ever receive

Every airliner overhead shouts its position twice a second at
1090 MHz. The grid is pure pulse timing — no carrier tracking, no
equalizer, just a microsecond ruler.

## The grid

| parameter | value | why |
|---|---|---|
| Carrier | 1090 MHz | worldwide allocation |
| Modulation | PPM (pulse-position), **1 µs per bit**, 0.5 µs pulses | amplitude-only: decodable with the crudest receiver imaginable |
| Preamble | pulses at **0, 1.0, 3.5, 4.5 µs**, then 3.5 µs of silence | an aperiodic pattern that can't false-trigger on data |
| Frame | 56 or **112 bits** (extended squitter = ADS-B) | 120 µs total on air |
| Bit encoding | pulse in 1st half-µs = 1, in 2nd = 0 | self-clocking |
| Integrity | **CRC-24** (poly 0xFFF409) | doubles as the addressing mechanism |
| Position | Compact Position Reporting (CPR), even/odd frame pairs | 17-bit lat/lon halves that interlock |

## Status: grid documented, atlas verification capture pending

Our ADS-B receiver ([aeroTuna](https://github.com/Felbs/aeroTuna))
decodes live traffic daily — including through a confidence-guided
"rescue" pass that recovers 2.7× more frames than hard slicing at
range. It decodes live rather than from files, so tonight the SDR is
busy with better things (an HD station that just decoded for the
first time ever) and this entry ships without its `figures/`.
`measure.py` is ready: point it at a 30 s capture at 2.4 MS/s centered
on 1090 MHz and it will find preambles by matched correlation, verify
CRC-24 on the hits, and plot the pulse timeline of a real squitter.

The measurement to expect: preamble pulse spacing 0/1.0/3.5/4.5 µs to
the sample, and a CRC pass rate that tells you honestly how many of
your "detections" are real (at 1090 MHz, impulse noise fakes preambles
constantly — the CRC is the referee).

## Reproduce it

```
python measure.py --iq capture.cs16 --fs 2400000
```
Any antenna works surprisingly well at 1090 MHz — we log airliners on
*rabbit ears*. Range scales with antenna quality and altitude.
