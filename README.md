# Celestial Messages Initiative (CMI)

Fictional teaching package for directed radio/laser messages to celestial targets, membership operations, equipment tracking, and Reed-Solomon decoding (Galois fields, Berlekamp-Massey, Chien search, Forney).

This does **not** represent a real transmitter network, FCC/ITU filing, or an actual METI campaign.

## Run the decoder demo

```bash
cd code
python3 demo.py
```

Requires only the Python 3 standard library.

Expected result: Galois arithmetic checks, the GF(2) Berlekamp-Massey example, a constructed two-root Chien search, and a successful RS(15,9) correction of two injected symbol errors.

## Layout

```
README.md
docs/
  01_System_Overview.md     CMI concept, tiers, laws, ocean site
  02_Encoding_Methods.md    radio + laser encodings
  03_Error_Correction.md    RS, BM, GF arithmetic, Chien, Forney
data/
  csv/                      spreadsheet sheets as CSV
code/
  galois_field.py           GF(2^4) and GF(2^8)
  berlekamp_massey.py       error-locator from syndromes
  chien_search.py           root scan of the locator
  reed_solomon.py           RS(15,9) encode / decode over GF(16)
  demo.py                   end-to-end runnable demo
```

The original Excel workbook (`CMI_Database.xlsx`) lives in the project folder that produced this repo; CSV exports of every sheet are under `data/csv/` so the tables stay diffable on GitHub.

## Disclaimer

Distance, power, and cost figures are order-of-magnitude teaching values, not an engineering budget. Interstellar one-way light time is years to millennia; no reply is assumed. High-power transmissions on Earth require national licenses and ITU coordination.
