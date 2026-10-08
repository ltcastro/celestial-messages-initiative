"""
Berlekamp-Massey algorithm over a Galois field.

Finds the shortest linear recurrence (error-locator polynomial) that
explains a syndrome sequence. Used as the core of Reed-Solomon decoding.

Returns coefficients [Λ0, Λ1, ..., ΛL] with Λ0 = 1 in a well-formed run.
"""

from __future__ import annotations

from galois_field import GaloisField, GF16


def _poly_scale(p: list[int], s: int, field: GaloisField) -> list[int]:
    return [field.mul(c, s) for c in p]


def _poly_add(a: list[int], b: list[int], field: GaloisField) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i, c in enumerate(a):
        out[i] = field.add(out[i], c)
    for i, c in enumerate(b):
        out[i] = field.add(out[i], c)
    return out


def berlekamp_massey(syndromes: list[int], field: GaloisField) -> tuple[list[int], int]:
    """
    Discrepancy-form BM.

    syndromes[0] is S1 (first syndrome), syndromes[k] is S_{k+1}.
    Returns locator Λ(x) = 1 + Λ1 x + ... + ΛL x^L  as [1, Λ1, ..., ΛL].
    """
    err_loc = [1]       # current Λ
    old_loc = [1]       # previous Λ
    for i, syn in enumerate(syndromes):
        # Δ = S_{i+1} + Σ_{j=1..deg} Λ_j S_{i+1-j}
        delta = syn
        for j in range(1, len(err_loc)):
            delta = field.add(delta, field.mul(err_loc[j], syndromes[i - j]))
        # shift old locator (multiply by x)
        old_loc = [0] + old_loc
        if delta != 0:
            if len(old_loc) > len(err_loc):
                new_loc = _poly_scale(old_loc, delta, field)
                old_loc = _poly_scale(err_loc, field.inv(delta), field)
                err_loc = new_loc
            err_loc = _poly_add(err_loc, _poly_scale(old_loc, delta, field), field)

    # trim leading zeros
    while len(err_loc) > 1 and err_loc[-1] == 0:
        err_loc.pop()
    return err_loc, len(err_loc) - 1


def berlekamp_massey_gf2(syndromes: list[int]) -> tuple[list[int], int]:
    """Binary (GF(2)) convenience wrapper used in the chat examples."""

    class _GF2:
        def add(self, a, b):
            return a ^ b

        def mul(self, a, b):
            return a & b

        def inv(self, a):
            if a == 0:
                raise ZeroDivisionError
            return 1

    return berlekamp_massey(syndromes, _GF2())  # type: ignore[arg-type]


if __name__ == "__main__":
    print("=== GF(2) example from the chat ===")
    s = [1, 0, 1, 1, 0, 1]
    poly, L = berlekamp_massey_gf2(s)
    print("syndromes:", s)
    print("Λ(x) coeffs:", poly, "degree", L)
    print(GF16(), "ready")
