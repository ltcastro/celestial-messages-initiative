"""
Small Reed-Solomon codec over GF(2^m). RS(n, k) with 2t = n-k parity symbols.

Polynomials are stored low-degree-first: p[i] is the coefficient of x^i.

Pipeline:
  encode -> add symbol errors -> syndromes -> Berlekamp-Massey
        -> Chien search -> Forney error values -> correct

Teaching implementation (not a production CCSDS encoder).
"""

from __future__ import annotations

from galois_field import GaloisField, GF16
from berlekamp_massey import berlekamp_massey


def _poly_mul(a, b, field):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] = field.add(out[i + j], field.mul(ai, bj))
    return out


def _poly_divmod(dividend, divisor, field):
    d = dividend[:]
    while d and d[-1] == 0:
        d.pop()
    divisor = divisor[:]
    while divisor and divisor[-1] == 0:
        divisor.pop()
    if not divisor:
        raise ZeroDivisionError("zero polynomial")
    quot = [0] * max(1, len(d) - len(divisor) + 1)
    inv_lead = field.inv(divisor[-1])
    while len(d) >= len(divisor) and d:
        shift = len(d) - len(divisor)
        coef = field.mul(d[-1], inv_lead)
        quot[shift] = coef
        for i, c in enumerate(divisor):
            d[shift + i] = field.add(d[shift + i], field.mul(coef, c))
        while d and d[-1] == 0:
            d.pop()
    return quot, d


class ReedSolomon:
    def __init__(self, field, n, k):
        if n > field.n_nonzero:
            raise ValueError("n cannot exceed 2^m - 1")
        if k >= n or k < 1:
            raise ValueError("need 1 <= k < n")
        self.field = field
        self.n = n
        self.k = k
        self.t = (n - k) // 2
        self.two_t = n - k
        self.g = [1]
        for i in range(1, self.two_t + 1):
            self.g = _poly_mul(self.g, [field.exp[i], 1], field)

    def encode(self, message):
        if len(message) != self.k:
            raise ValueError(f"message must have {self.k} symbols")
        shifted = [0] * self.two_t + list(message)
        _, rem = _poly_divmod(shifted, self.g, self.field)
        rem = rem + [0] * (self.two_t - len(rem))
        return rem[: self.two_t] + list(message)

    def syndromes(self, received):
        return [
            self.field.evaluate_poly(received, self.field.exp[j])
            for j in range(1, self.two_t + 1)
        ]

    def forney(self, locator, syndromes, positions):
        gf = self.field
        omega = [0] * self.two_t
        for i, si in enumerate(syndromes):
            for j, lj in enumerate(locator):
                if i + j < self.two_t:
                    omega[i + j] = gf.add(omega[i + j], gf.mul(si, lj))
        values = []
        for pos in positions:
            Xi = gf.exp[pos]
            Xi_inv = gf.inv(Xi)
            num = gf.evaluate_poly(omega, Xi_inv)
            deriv = 0
            xp = 1
            for idx, c in enumerate(locator):
                if idx % 2 == 1:
                    deriv = gf.add(deriv, gf.mul(c, xp))
                if idx >= 1:
                    xp = gf.mul(xp, Xi_inv)
            if deriv == 0:
                raise ValueError("Forney derivative vanished")
            values.append(gf.div(num, deriv))
        return values

    def decode(self, received):
        if len(received) != self.n:
            raise ValueError(f"received must have {self.n} symbols")
        s = self.syndromes(received)
        if all(v == 0 for v in s):
            return received[:], []
        locator, L = berlekamp_massey(s, self.field)
        positions = []
        gf = self.field
        for i in range(self.n):
            x_inv = 1 if i == 0 else gf.exp[(self.n - i) % self.n]
            if gf.evaluate_poly(locator, x_inv) == 0:
                positions.append(i)
        if len(positions) != L:
            raise ValueError(f"Chien found {len(positions)} roots but BM degree is {L}")
        errors = self.forney(locator, s, positions)
        corrected = received[:]
        for pos, e in zip(positions, errors):
            corrected[pos] = gf.add(corrected[pos], e)
        if any(v != 0 for v in self.syndromes(corrected)):
            raise ValueError(f"syndromes still nonzero after correction (locator={locator}, pos={positions})")
        return corrected, positions


if __name__ == "__main__":
    gf = GF16()
    rs = ReedSolomon(gf, n=15, k=9)
    msg = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    cw = rs.encode(msg)
    print("g(x):", rs.g)
    print("message :", msg)
    print("codeword:", cw)
    print("cw syndromes:", rs.syndromes(cw))
    received = cw[:]
    received[2] = gf.add(received[2], 7)
    received[11] = gf.add(received[11], 4)
    print("received:", received)
    decoded, pos = rs.decode(received)
    print("decoded :", decoded)
    print("positions:", pos)
    print("success:", decoded == cw)
