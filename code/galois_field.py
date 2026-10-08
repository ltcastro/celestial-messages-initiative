"""
Galois Field arithmetic for Reed-Solomon / Berlekamp-Massey / Chien search.

Provides:
  - GF16  : GF(2^4) with primitive polynomial x^4 + x + 1  (0b10011)
  - GF256 : GF(2^8) with primitive polynomial x^8 + x^4 + x^3 + x^2 + 1 (0x11d)
            — the field used by classic RS(255,223) / CCSDS / Voyager-style codes.

All operations are integer-valued field elements in [0, 2^m - 1].
Addition is XOR. Multiplication uses log/antilog tables.
"""

from __future__ import annotations


class GaloisField:
    """Generic GF(2^m) with log/antilog tables."""

    def __init__(self, m: int, primitive_poly: int):
        if m < 2 or m > 16:
            raise ValueError("Supported range is GF(2^2) .. GF(2^16)")
        self.m = m
        self.order = 1 << m          # 2^m
        self.n_nonzero = self.order - 1
        self.primitive_poly = primitive_poly
        self.exp = [0] * (self.n_nonzero * 2)  # doubled for wrap-free multiply
        self.log = [0] * self.order
        self._build_tables()

    def _build_tables(self) -> None:
        x = 1
        high_bit = 1 << self.m
        for i in range(self.n_nonzero):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & high_bit:
                x ^= self.primitive_poly
        # duplicate so exp[i + n] = exp[i] without modulo
        for i in range(self.n_nonzero):
            self.exp[i + self.n_nonzero] = self.exp[i]
        self.log[0] = -1  # undefined

    def add(self, a: int, b: int) -> int:
        return a ^ b

    def sub(self, a: int, b: int) -> int:
        return a ^ b  # characteristic 2

    def mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def div(self, a: int, b: int) -> int:
        if b == 0:
            raise ZeroDivisionError("division by zero in Galois field")
        if a == 0:
            return 0
        return self.exp[self.log[a] - self.log[b] + self.n_nonzero]

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("zero has no inverse")
        return self.exp[self.n_nonzero - self.log[a]]

    def pow(self, a: int, e: int) -> int:
        if a == 0:
            return 0 if e > 0 else 1
        e %= self.n_nonzero
        return self.exp[(self.log[a] * e) % self.n_nonzero]

    def evaluate_poly(self, coeffs: list[int], x: int) -> int:
        """Horner evaluation of sum coeffs[i] * x^i."""
        result = 0
        xp = 1
        for c in coeffs:
            result = self.add(result, self.mul(c, xp))
            xp = self.mul(xp, x)
        return result

    def __repr__(self) -> str:
        return f"GF(2^{self.m}) prim=0x{self.primitive_poly:x}"


def GF16() -> GaloisField:
    """GF(2^4), primitive polynomial x^4 + x + 1."""
    return GaloisField(4, 0b10011)


def GF256() -> GaloisField:
    """GF(2^8), primitive polynomial 0x11d (CCSDS / AES-adjacent)."""
    return GaloisField(8, 0x11D)


if __name__ == "__main__":
    gf = GF16()
    print(gf)
    print("3 + 5 =", gf.add(3, 5))
    print("3 * 5 =", gf.mul(3, 5))
    print("inv(3) =", gf.inv(3), "  check 3*inv =", gf.mul(3, gf.inv(3)))
    gf8 = GF256()
    print(gf8)
    print("7 * 11 =", gf8.mul(7, 11), "  inv(7)*7 =", gf8.mul(gf8.inv(7), 7))
