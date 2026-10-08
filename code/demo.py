"""
CMI end-to-end demo:
  1. Galois field arithmetic
  2. Berlekamp-Massey (GF(2) chat example + field version)
  3. Chien search on a constructed locator
  4. Full RS(15,9) encode / corrupt / decode
"""

from galois_field import GF16, GF256
from berlekamp_massey import berlekamp_massey, berlekamp_massey_gf2
from chien_search import chien_search
from reed_solomon import ReedSolomon


def section(title: str) -> None:
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)


def main() -> None:
    section("1. Galois field arithmetic")
    gf16 = GF16()
    print(gf16)
    print("  3 + 5 =", gf16.add(3, 5), "(XOR)")
    print("  3 * 5 =", gf16.mul(3, 5))
    print("  inv(3) =", gf16.inv(3), "  3*inv(3) =", gf16.mul(3, gf16.inv(3)))
    gf256 = GF256()
    print(gf256)
    print("  7 * 11 =", gf256.mul(7, 11), "  7*inv(7) =", gf256.mul(7, gf256.inv(7)))

    section("2. Berlekamp-Massey — GF(2) chat example")
    s = [1, 0, 1, 1, 0, 1]
    poly, L = berlekamp_massey_gf2(s)
    print("  syndromes:", s)
    print("  Lambda(x) coeffs [const .. high]:", poly, "  degree L =", L)

    section("3. Chien search on constructed two-error locator")
    a3, a7 = gf16.exp[3], gf16.exp[7]
    locator = [1, gf16.add(a3, a7), gf16.mul(a3, a7)]
    print("  Lambda(x) = (1 + alpha^3 x)(1 + alpha^7 x)  coeffs:", locator)
    roots = chien_search(locator, gf16, verbose=True)
    expected = sorted([(15 - 3) % 15, (15 - 7) % 15])
    print("  found roots i:", roots)
    print("  expected alpha^{-3}, alpha^{-7} indices:", expected)

    section("4. Full RS(15,9) encode -> 2 errors -> decode")
    rs = ReedSolomon(gf16, n=15, k=9)
    msg = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    cw = rs.encode(msg)
    received = cw[:]
    received[2] = gf16.add(received[2], 7)
    received[11] = gf16.add(received[11], 4)
    print("  message :", msg)
    print("  codeword:", cw)
    print("  received:", received)
    try:
        decoded, pos = rs.decode(received)
        print("  decoded :", decoded)
        print("  positions:", pos)
        print("  SUCCESS" if decoded == cw else "  MISMATCH")
    except ValueError as exc:
        print("  decode failed:", exc)
        print("  (Forney/syndrome convention may need a tweak; BM+Chien still stand.)")


if __name__ == "__main__":
    main()
