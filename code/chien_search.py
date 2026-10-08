"""
Chien search: evaluate the error-locator polynomial at every field element
α^i and report the indices i where Λ(α^i) = 0.

Convention used here (narrow-sense RS over GF(2^m) of length n = 2^m - 1):
  A root at α^i means an error in symbol position i
  (0-based index into a length-n codeword whose positions are α^0 .. α^{n-1}).
"""

from __future__ import annotations

from galois_field import GaloisField, GF16


def chien_search(
    locator_coeffs: list[int],
    field: GaloisField,
    n: int | None = None,
    verbose: bool = False,
) -> list[int]:
    """
    Parameters
    ----------
    locator_coeffs : coefficients [Λ0, Λ1, ..., Λv]
    field          : GaloisField used to build the locator
    n              : code length; defaults to field.n_nonzero (2^m - 1)
    verbose        : print the evaluation table

    Returns
    -------
    list of error position indices i where Λ(α^i) = 0
    """
    if n is None:
        n = field.n_nonzero
    roots: list[int] = []
    if verbose:
        print("i   α^i   Λ(α^i)   root?")
        print("-" * 32)
    for i in range(n):
        alpha_i = field.exp[i]
        val = field.evaluate_poly(locator_coeffs, alpha_i)
        is_root = val == 0
        if verbose:
            print(f"{i:3d}  {alpha_i:4d}    {val:4d}    {'YES' if is_root else 'no'}")
        if is_root:
            roots.append(i)
    return roots


if __name__ == "__main__":
    gf = GF16()
    # Construct Λ(x) = (1 + α^3 x)(1 + α^7 x)
    a3 = gf.exp[3]
    a7 = gf.exp[7]
    lam0 = 1
    lam1 = gf.add(a3, a7)
    lam2 = gf.mul(a3, a7)
    locator = [lam0, lam1, lam2]
    print("Locator coeffs:", locator)
    roots = chien_search(locator, gf, verbose=True)
    print("Error positions i:", roots)
    # Roots of (1 + α^k x) are x = α^{-k} = α^{15-k}
    expected = sorted([(15 - 3) % 15, (15 - 7) % 15])
    print("Expected (α^{-3}, α^{-7} indices):", expected)
    print("Match:", sorted(roots) == expected)
