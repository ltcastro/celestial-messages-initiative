# Message Encoding Methods

Goal: a signal that is (1) obviously artificial, (2) robust to noise and
dispersion, and (3) decodable from mathematics / physics without sharing a
human language first.

## Principles

- Primer first: primes, counting, hydrogen line, simple arithmetic.
- Redundancy: repeats + forward error correction (Reed-Solomon, LDPC).
- Multi-modal payload after the primer: text, bitmap, audio spectrogram.
- Radio for a wide hailing beam; laser for a high-rate directed payload.

## Radio

Bands: the water hole near 1.42-1.66 GHz is a traditional quiet window
(21 cm H I to OH), not a legal requirement.

Modulation (simple to denser):

- OOK / ASK — carrier on/off
- FSK — two tones
- BPSK / QPSK / MSK — phase or continuous-phase, used on real deep-space links

Message structure:

- Arecibo-style 2-D binary bitmap whose dimensions are primes (e.g. 23x73)
- Lone Signal / Busch-style binary language with a defined hailing lexicon
- Zaitsev-style multi-section spectral languages (constant / continuous / discrete)

Deep-space practice concatenates an inner convolutional or LDPC code with an
outer Reed-Solomon block (Voyager, CCSDS heritage).

## Laser / optical

- OOK and pulse-position modulation (PPM)
- Serially concatenated PPM (SCPPM): convolutional outer code + accumulator + PPM
- Phase / polarization / wavelength as extra dimensions
- Adaptive optics on the uplink to fight atmosphere

Lasers buy beam gain (small diffraction angle) and bandwidth. They demand
precise ephemeris and a cooperative (or very large) collector at the far end.

## CMI encoding profiles by tier

| Tier | Profile |
|------|---------|
| 1 | Binary OOK or FSK, short RS, text only |
| 2 | BPSK + RS(shortened), 2-D prime bitmap |
| 3 | Arecibo-style grid + Reed-Solomon + laser PPM fallback |
| 4 | Busch binary + multi-modal (image + audio) + SCPPM laser + full RS(255,223)-class FEC |

Each outbound record stores an encoding_profile so a future decoder knows
which stack to try first.
