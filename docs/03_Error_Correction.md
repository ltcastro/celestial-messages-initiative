# Error Correction Stack

Implementation lives in `code/`.

## Reed-Solomon

Teaching codec: RS(15, 9) over GF(16), t = 3. Classic deep-space form is RS(255, 223) over GF(256), t = 16.

Generator: g(x) = product from i=1 to 2t of (x - alpha^i).

## Galois field

Addition is XOR. Multiplication uses log/antilog tables. See `code/galois_field.py`.

- GF(2^4), primitive polynomial x^4 + x + 1 (0b10011)
- GF(2^8), primitive polynomial 0x11d

## Pipeline

1. Syndromes S_j = r(alpha^j) for j = 1..2t
2. Berlekamp-Massey builds the error-locator Lambda(x)
3. Chien search finds roots; a root at alpha^{-i} is an error in symbol i
4. Forney computes error values e_i = Omega(X_i^{-1}) / Lambda'(X_i^{-1})

Run `python3 code/demo.py` to encode, inject two symbol errors, and recover the block.
