# A retired Helium miner as a weather receiver — the LoRa card's spare FSK modem

The Semtech SX1302 inside every Helium IoT hotspot is a LoRa
concentrator, but it also carries one **plain FSK demodulator** (LoRaWAN
has an FSK data rate). With the vendor HAL patched to hand over raw
bits — fixed-length blocks, no CRC, no whitening, no length byte, the
same patch the [zwave-908](../zwave-908/) entry used — that modem can be
pointed at anything in 902–928 MHz that is 2-FSK with a known bit rate
and sync word. Here it reads a neighbourhood weather station, all day,
on hardware that was mining a cryptocurrency until last year.

## The recipe

| element | value | why |
|---|---|---|
| Card | RAK5146 (SX1302 + SX1250), USB, inside a FreedomFi FFG-HL-4-64 hotspot | the same card Helium gateways use |
| RF chain | 915.2 MHz, FSK IF channel at −200 kHz → **915.000 MHz** | the station sits at 915.03 MHz |
| Modem | FSK, bandwidth 125 kHz, **17,241 bit/s** (58 µs) | rtl_433's figure for the Fine Offset WH65B |
| Sync word | **`AA 2D D4`** (3 bytes) | the last preamble byte plus the station's sync, so the block starts at the family byte `0x24` |
| Block | fixed **40 bytes** | the 17-byte frame plus whatever follows; the host finds the frame |
| Integrity | software: **CRC-8 poly 0x31** over bytes 0–14, **8-bit sum** over bytes 0–15 | the modem cannot check it; the host does |
| Transmit | never | listen-only |

## What it hears (2.9 h, attic antenna, −108 … −118 dBm)

```
LoRa-card log: 572 CRC+sum-verified frames in 2.88 h = 199/h of 225 possible (88 %)
cadence from receive timestamps: 16.000 s per frame; consecutive 508, one missed 57, two missed 5
rssi dBm: min/median/max -118 / -117 / -108
```

The timestamps are whole seconds: the 16 s comb slips one second after
2.1 h, so the station's true period is about 15.998 s (130 ppm fast),
and the 249 ms rms is mostly the one-second quantization plus USB latency.

The same station on an SDR with a rooftop discone, 70 s, high gain:
bursts 7–8 dB above the floor exactly on this 16 s comb, **no frame
decodes**. The attic antenna's 8 dBi and a modem that integrates whole
bits beat a general-purpose receiver's wider view. The number to carry
away: **88 % of frames verified at −117 dBm.**

## What the card cannot do, found the hard way

- **No sample access.** The modem returns bits or nothing; deviation,
  spectrum and timing finer than a bit are invisible. The weather
  entry's deviation for this station is blank for that reason.
- **One setting at a time.** Bit rate, sync word and frequency are
  fixed per listen. Searching for an unknown station means guessing:
  the same card heard nothing from a 433 MHz station on 915 or 868 MHz
  at four bit rates and two polarities, because it was on the wrong band
  (the front-end filter and the antenna are built for 902–928 MHz).
- **No on-off keying.** The Itron ERT and Neptune R900 meters in the same
  band are OOK; the FSK modem has no mode for them.
- **Fixed-length blocks swallow neighbours.** After a sync the modem
  reads the full block (40 bytes = 18.6 ms here, 64 bytes = 12.8 ms in
  the Z-Wave work) before it can re-arm; a second frame inside that
  window is lost. Z-Wave's thermostat replies, 10 ms behind the
  controller's beam, were lost exactly this way.
- **Polarity.** The SX1302 marks a logical 1 on the upper tone; Z-Wave
  puts it on the lower one, this station on the upper. `wxrx1302`
  accepts `inv` to try the other.

## Files

- `wxrx1302.c` — the receiver: configure the card, print every block as
  hex with its RSSI. Build against the patched `sx1302_hal` from the
  zwave-908 entry (`lgw_fsk_raw`).
- `decode_wh65.py` — reads that output, applies the CRC and the sum,
  prints the fields with the ID masked.

```
./wxrx1302 3600 | python3 decode_wh65.py
```

**Privacy:** the station is a neighbour's. Its ID is masked everywhere
here and is not in this repository; the readings are not published.
