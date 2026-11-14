"""Reproduce every number in front-matter chapter f2 of Error Correction Coding, "Setting Up the Lab Bench".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/chf2_ecc_lab_bench.py
    python3 python/chf2_ecc_lab_bench.py --bench

The chapter installs the book's tools and then teaches the Python the listings
use, so its numbers are of three kinds. The first is integer and bit arithmetic:
binary forms, shifts, masks, floor division and remainders with negative
operands. The second is floating point: the exact value stored for 0.1, the
spacing of doubles, the sum 0.1 + 0.2, rounding half to even, and the half-unit
test that every companion script's check() helper applies. The third is
version numbers, compared as tuples, and the compatible-release rule that pip
reads in "cocotb~=2.1". A short NumPy-style preview closes the chapter with a
mod-2 matrix product and one product in the field of 256 elements, both done
here in plain Python. Each theorem the chapter proves is also checked by brute
force over a range of inputs.

Every value computed here is printed beside the value the chapter states. The
script stops with exit status 1 at the first pair that disagrees by more than
half a unit in the last digit the chapter prints, and at the first count, list
or string that differs at all.

With --bench the script also inspects the machine it runs on. It asks each
tool of the lab bench for its version, compares the answer with the oldest
version the book's listings need where an official source states one, and
repeats the NumPy and galois checks with those libraries when they are
installed. That part exits 1 when anything the companion scripts and the first
hardware lab need is missing or too old: Python outside 3.10 to 3.14, git,
Verilator, a C++ compiler, make, Yosys, Surfer, cocotb, NumPy or galois. Icarus
Verilog, SymbiYosys and Z3 are reported but never fail the bench, since a later
Part is the first to ask for them.
"""

import argparse
import re
import shutil
import subprocess
import sys
from decimal import Decimal
from fractions import Fraction
from itertools import product
from math import inf, nextafter, ulp


def half_unit(stated: str) -> float:
    text = stated.replace(",", "").rstrip("%")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def agrees(computed, stated: str) -> bool:
    """The test check() applies: within half a unit of the last digit stated."""
    value = float(stated.replace(",", "").rstrip("%"))
    if stated.endswith("%"):
        computed = 100 * float(computed)
    return abs(float(computed) - value) <= half_unit(stated) * (1 + 1e-9)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", "").rstrip("%"))
    if stated.endswith("%"):
        computed = 100 * float(computed)
    print(f"  {label:<64} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<64} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


# Integers and bits ------------------------------------------------------------

def check_integers_and_bits() -> None:
    print("Integers and bits")
    same("bin(100)", bin(100), "0b1100100")
    same("format(100, '08b')", format(100, "08b"), "01100100")
    same("int('01100100', 2)", int("01100100", 2), 100)
    same("100 ^ 0b10000000", 100 ^ 0b10000000, 228)
    same("100 >> 2", 100 >> 2, 25)
    same("100 << 1", 100 << 1, 200)
    same("100 & 1, the lowest bit", 100 & 1, 0)
    same("(100 >> 2) & 1, bit 2", (100 >> 2) & 1, 1)
    same("100 | 1", 100 | 1, 101)
    same("ones in 100, bin(100).count('1')", bin(100).count("1"), 3)
    same("0x64, the same number in hexadecimal", 0x64, 100)
    same("2 ** 100, exactly", 2 ** 100, 1267650600228229401496703205376)
    same("digits in 2 ** 100", len(str(2 ** 100)), 31)
    same("100 ^ 29", 100 ^ 29, 121)
    same("121 ^ 29, undoing the XOR", 121 ^ 29, 100)


def check_shift_theorem() -> None:
    """n >> k is floor(n / 2^k) and n << k is n * 2^k, for every integer n and k >= 0."""
    print("Shifts are floor division and multiplication by powers of 2")
    cases = 0
    for n in range(-1024, 1025):
        for k in range(12):
            assert n >> k == n // 2 ** k, (n, k)
            assert n << k == n * 2 ** k, (n, k)
            cases += 1
    same("cases checked, n from -1024 to 1024, k from 0 to 11", cases, 24588)
    same("-7 >> 1", -7 >> 1, -4)
    same("13 >> 2", 13 >> 2, 3)


def check_floor_division() -> None:
    print("Floor division and remainder")
    same("7 // 2", 7 // 2, 3)
    same("7 % 2", 7 % 2, 1)
    same("-7 // 2", -7 // 2, -4)
    same("-7 % 2", -7 % 2, 1)
    same("7 // -2", 7 // -2, -4)
    same("7 % -2", 7 % -2, -1)
    same("int(-7 / 2), truncation toward zero", int(-7 / 2), -3)
    same("divmod(17, 5)", divmod(17, 5), (3, 2))
    same("-17 % 5", -17 % 5, 3)
    same("1000 % 7", 1000 % 7, 6)
    pairs = 0
    for a in range(-200, 201):
        for b in range(1, 21):
            q, r = a // b, a % b
            assert a == b * q + r and 0 <= r < b, (a, b)
            # The quotient and remainder are the only pair with that property.
            assert sum(1 for qq in range(q - 3, q + 4) if 0 <= a - b * qq < b) == 1
            pairs += 1
    same("pairs checked, a from -200 to 200, b from 1 to 20", pairs, 8020)


# Floating point ---------------------------------------------------------------

def binary_digits(x: Fraction, n: int) -> str:
    """The first n bits after the binary point of a fraction in [0, 1)."""
    out = []
    for _ in range(n):
        x *= 2
        bit = int(x >= 1)
        out.append(str(bit))
        x -= bit
    return "".join(out)


def terminates(x: Fraction, limit: int = 64) -> bool:
    """Does repeated doubling clear the fraction within the limit?"""
    x -= int(x)
    for _ in range(limit):
        if x == 0:
            return True
        x *= 2
        x -= int(x)
    return x == 0


def check_finite_expansion_theorem() -> None:
    """p/q in lowest terms has a finite binary expansion exactly when q is a power of 2."""
    print("Finite binary expansions")
    tried = 0
    for q in range(1, 257):
        power = q & (q - 1) == 0
        for p in range(q):
            x = Fraction(p, q)
            if x.denominator != q:
                continue
            assert terminates(x) == power, (p, q)
            tried += 1
    same("fractions p/q in lowest terms checked, q up to 256", tried, 19948)
    same("first 24 bits of 1/10", binary_digits(Fraction(1, 10), 24), "000110011001100110011001")
    same("first 8 bits of 1/10, as a fraction", Fraction(int(binary_digits(Fraction(1, 10), 8), 2), 256), Fraction(25, 256))
    check("25/256", 25 / 256, "0.09765625")
    check("error of the 8-bit truncation", 0.1 - 25 / 256, "0.00234375")
    same("does 3/8 terminate", terminates(Fraction(3, 8)), True)
    same("3/8 in binary", binary_digits(Fraction(3, 8), 3), "011")
    same("does 1/3 terminate", terminates(Fraction(1, 3)), False)
    same("first 8 bits of 1/3", binary_digits(Fraction(1, 3), 8), "01010101")


def check_point_one() -> None:
    print("The double nearest 0.1")
    stored = Fraction(0.1)
    same("Fraction(0.1)", stored, Fraction(3602879701896397, 36028797018963968))
    same("its denominator is 2 ** 55", stored.denominator == 2 ** 55, True)
    same("2 ** 55", 2 ** 55, 36028797018963968)
    same("2 ** 56 / 10, as a fraction", Fraction(2 ** 56, 10), Fraction(36028797018963968, 5))
    check("2 ** 56 / 10", Fraction(2 ** 56, 10), "7205759403792793.6")
    same("nearest integer to 2 ** 56 / 10", round(Fraction(2 ** 56, 10)), 7205759403792794)
    same("7205759403792794 / 2 ** 56 in lowest terms", Fraction(7205759403792794, 2 ** 56), stored)
    same("2 ** 52 <= 7205759403792794 < 2 ** 53", 2 ** 52 <= 7205759403792794 < 2 ** 53, True)
    same("stored minus 1/10, exactly", stored - Fraction(1, 10), Fraction(1, 5 * 2 ** 55))
    check("stored minus 1/10", stored - Fraction(1, 10), "5.551e-18")
    same("Decimal(0.1)", str(Decimal(0.1)), "0.1000000000000000055511151231257827021181583404541015625")
    same("repr(0.1)", repr(0.1), "0.1")
    same("0.1 + 0.2", repr(0.1 + 0.2), "0.30000000000000004")
    same("0.1 + 0.2 == 0.3", 0.1 + 0.2 == 0.3, False)
    same("Decimal(0.1 + 0.2)", str(Decimal(0.1 + 0.2)), "0.3000000000000000444089209850062616169452667236328125")
    same("Decimal(0.3)", str(Decimal(0.3)), "0.299999999999999988897769753748434595763683319091796875")
    same("(0.1 + 0.2) - 0.3 is 2 ** -54", (0.1 + 0.2) - 0.3 == 2.0 ** -54, True)
    check("2 ** -54", 2.0 ** -54, "5.551e-17")
    same("0.1 + 0.2 is the double just above 0.3", nextafter(0.3, inf) == 0.1 + 0.2, True)
    same("abs((0.1 + 0.2) - 0.3) < 1e-9", abs((0.1 + 0.2) - 0.3) < 1e-9, True)


def check_spacing_theorem() -> None:
    """Consecutive doubles in [2^e, 2^(e+1)) are 2^(e-52) apart."""
    print("The spacing of doubles")
    for e in range(-60, 61):
        x = 2.0 ** e
        gap = 2.0 ** (e - 52)
        assert ulp(x) == gap, e
        assert nextafter(x, inf) - x == gap, e
        y = nextafter(2.0 ** (e + 1), 0.0)
        assert y - nextafter(y, 0.0) == gap, e
    same("intervals [2**e, 2**(e+1)) checked, e from -60 to 60", len(range(-60, 61)), 121)
    same("ulp(1.0) is 2 ** -52", ulp(1.0) == 2.0 ** -52, True)
    check("ulp(1.0)", ulp(1.0), "2.220e-16")
    same("ulp(0.1) is 2 ** -56", ulp(0.1) == 2.0 ** -56, True)
    check("ulp(0.1)", ulp(0.1), "1.388e-17")
    same("ulp(1000.0) is 2 ** -43", ulp(1000.0) == 2.0 ** -43, True)
    check("ulp(1000.0)", ulp(1000.0), "1.137e-13")
    same("ulp(2.0 ** 53) is 2", ulp(2.0 ** 53), 2.0)
    same("2 ** 53", 2 ** 53, 9007199254740992)
    same("float(2 ** 53) + 1.0 == float(2 ** 53)", float(2 ** 53) + 1.0 == float(2 ** 53), True)
    same("float(2 ** 53 + 1) == float(2 ** 53)", float(2 ** 53 + 1) == float(2 ** 53), True)
    same("float(2 ** 53 + 2) - float(2 ** 53)", float(2 ** 53 + 2) - float(2 ** 53), 2.0)
    same("doubles in [1, 2)", 2 ** 52, 4503599627370496)


def check_rounding() -> None:
    print("Rounding half to even")
    same("round(2.5)", round(2.5), 2)
    same("round(3.5)", round(3.5), 4)
    same("round(-2.5)", round(-2.5), -2)
    same("round(0.5)", round(0.5), 0)
    same("round(1.5)", round(1.5), 2)
    same("round(0.125, 2)", round(0.125, 2), 0.12)
    same("round(0.375, 2)", round(0.375, 2), 0.38)
    same("0.125 is stored exactly", Fraction(0.125), Fraction(1, 8))
    same("round(2.675, 2)", round(2.675, 2), 2.67)
    same("Decimal(2.675)", str(Decimal(2.675)), "2.67499999999999982236431605997495353221893310546875")
    same("f'{2.675:.2f}'", f"{2.675:.2f}", "2.67")
    same("f'{0.96583769633507:.5f}'", f"{0.96583769633507:.5f}", "0.96584")
    # Ties to even leave no drift on average: the four ties 0.5, 1.5, 2.5, 3.5
    # round to 0, 2, 2, 4, which sum to 8, the same as the ties themselves.
    same("sum of round(t) for t in 0.5, 1.5, 2.5, 3.5", sum(round(t) for t in (0.5, 1.5, 2.5, 3.5)), 8)
    same("sum of 0.5, 1.5, 2.5, 3.5", 0.5 + 1.5 + 2.5 + 3.5, 8.0)
    same("sum when every tie rounds up", sum(int(t + 0.5) for t in (0.5, 1.5, 2.5, 3.5)), 10)


# The half-unit test -----------------------------------------------------------

def check_half_unit_theorem() -> None:
    """A d-place decimal s is a correct rounding of x exactly when |x - s| <= 10^-d / 2."""
    print("The half-unit test")
    same("half_unit('0.96584')", half_unit("0.96584"), 0.5 * 10.0 ** -5)
    check("half_unit('0.96584')", half_unit("0.96584"), "0.000005")
    check("half_unit('7,379')", half_unit("7,379"), "0.5")
    check("half_unit('2.220e-16')", half_unit("2.220e-16"), "5e-20")
    check("half_unit('96.58%')", half_unit("96.58%"), "0.005")
    check("7379 / 7640", 7379 / 7640, "0.96584")
    same("7379 / 7640 agrees with '0.96584'", agrees(7379 / 7640, "0.96584"), True)
    same("7379 / 7640 agrees with '0.9658'", agrees(7379 / 7640, "0.9658"), True)
    same("7379 / 7640 agrees with '0.9659'", agrees(7379 / 7640, "0.9659"), False)
    same("7379 / 7640 agrees with '96.58%'", agrees(7379 / 7640, "96.58%"), True)
    check("distance from 0.9659", abs(7379 / 7640 - 0.9659), "0.0000623")
    # Brute force: every x on a fine grid, every d from 0 to 3. A d-place decimal
    # passes the test exactly when it is a nearest d-place decimal to x.
    tried = 0
    for d in range(4):
        unit = Fraction(1, 10 ** d)
        for i in range(-400, 401):
            x = Fraction(i, 160)
            nearest = min(abs(x - unit * j) for j in range(-30 * 10 ** d, 30 * 10 ** d + 1))
            low = (x / unit).__floor__()
            for j in (low - 1, low, low + 1, low + 2):
                s = unit * j
                passes = abs(x - s) <= unit / 2
                assert passes == (abs(x - s) == nearest), (x, d, s)
            tried += 1
    same("grid points checked, 801 points times 4 precisions", tried, 3204)


# Version numbers --------------------------------------------------------------

def version(text: str) -> tuple:
    """The first dotted number in a tool's output, padded to three parts."""
    match = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text)
    if not match:
        raise ValueError(f"no version number in {text!r}")
    return tuple(int(g or 0) for g in match.groups())


def compatible(v: tuple, spec: str) -> bool:
    """PEP 440's ~= operator: at least spec, and the same release series."""
    parts = tuple(int(x) for x in spec.split("."))
    prefix = parts[:-1]
    padded = v + (0,) * (len(parts) - len(v))
    return padded >= parts and padded[:len(prefix)] == prefix


def check_versions() -> None:
    print("Version numbers")
    same("version('Verilator 5.054 2026-11-01')", version("Verilator 5.054 2026-11-01 rev vUNKNOWN"), (5, 54, 0))
    same("version('Yosys 0.71+post (git sha1 ...)')", version("Yosys 0.71+post (git sha1 5628c13)"), (0, 71, 0))
    same("version('SBY v0.71')", version("SBY v0.71"), (0, 71, 0))
    same("version('surfer 0.7.0 (git: v0.7.0)')", version("surfer 0.7.0 (git: v0.7.0)"), (0, 7, 0))
    same("version('Icarus Verilog version 13.0 (stable)')", version("Icarus Verilog version 13.0 (stable) (v13_0)"), (13, 0, 0))
    same("version('git version 2.56.1')", version("git version 2.56.1"), (2, 56, 1))
    same("Verilator 5.054 meets the floor 5.036", version("5.054") >= version("5.036"), True)
    same("Verilator 5.020 meets the floor 5.036", version("5.020") >= version("5.036"), False)
    same("'3.10' < '3.9' as strings", "3.10" < "3.9", True)
    same("3.10 < 3.9 as decimals", 3.10 < 3.9, True)
    same("(3, 10) < (3, 9) as tuples", (3, 10) < (3, 9), False)
    same("'0.71' < '0.9' as strings", "0.71" < "0.9", True)
    same("(0, 71) > (0, 9) as tuples", (0, 71) > (0, 9), True)
    same("(3, 10) < (3, 10, 0), Python's prefix rule", (3, 10) < (3, 10, 0), True)
    same("padded, (3, 10, 0) == (3, 10, 0)", version("3.10") == version("3.10.0"), True)
    same("sys.version_info compares as a tuple", sys.version_info >= (3, 10), True)
    for v, spec, want in [((2, 0, 0), "2.1", False), ((2, 1, 0), "2.1", True), ((2, 1, 5), "2.1", True),
                          ((2, 9, 3), "2.1", True), ((3, 0, 0), "2.1", False),
                          ((2, 1, 0), "2.1.0", True), ((2, 1, 7), "2.1.0", True), ((2, 2, 0), "2.1.0", False)]:
        same(f"{'.'.join(map(str, v))} satisfies ~={spec}", compatible(v, spec), want)
    # Lexicographic order on padded tuples is total and transitive.
    triples = list(product(range(3), repeat=3))
    for a in triples:
        for b in triples:
            assert (a < b) + (a == b) + (a > b) == 1, (a, b)
            for c in triples:
                if a < b and b < c:
                    assert a < c, (a, b, c)
    same("triples checked for totality and transitivity", len(triples), 27)


# A preview of NumPy and galois, in plain Python --------------------------------

G = [[1, 0, 0, 0, 1, 1, 0],
     [0, 1, 0, 0, 1, 0, 1],
     [0, 0, 1, 0, 0, 1, 1],
     [0, 0, 0, 1, 1, 1, 1]]


def encode(m: list) -> list:
    """The row vector m times G, reduced mod 2: what (m @ G) % 2 computes in NumPy."""
    return [sum(m[i] * G[i][j] for i in range(len(G))) % 2 for j in range(len(G[0]))]


def clmul(a: int, b: int) -> int:
    """Carry-less product: polynomial multiplication with coefficients mod 2."""
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


def reduce(x: int, poly: int = 0x11D) -> int:
    """The remainder of x after division by the degree-8 polynomial poly, mod 2."""
    while x.bit_length() > 8:
        x ^= poly << (x.bit_length() - 9)
    return x


def poly_add(p: list, q: list) -> list:
    """Coefficient-wise sum mod 2 of two coefficient lists, lowest power first."""
    return [(a + b) % 2 for a, b in zip(p, q)]


def check_preview() -> None:
    print("A preview of NumPy and galois")
    same("m = [1, 0, 1, 1] times G, mod 2", encode([1, 0, 1, 1]), [1, 0, 1, 1, 0, 1, 0])
    same("the sum before reduction", [sum(r[j] for r, mi in zip(G, [1, 0, 1, 1]) if mi) for j in range(7)],
         [1, 0, 1, 1, 2, 3, 2])
    same("[1, 1, 0, 0] times G, mod 2", encode([1, 1, 0, 0]), [1, 1, 0, 0, 0, 1, 1])
    codewords = {tuple(encode(list(m))) for m in product((0, 1), repeat=4)}
    same("distinct codewords from 16 messages", len(codewords), 16)
    same("smallest weight of a nonzero codeword", min(sum(c) for c in codewords if any(c)), 3)
    same("0x11D in binary", bin(0x11D), "0b100011101")
    same("100 in binary", bin(100), "0b1100100")
    same("29 in binary", bin(29), "0b11101")
    same("carry-less 100 times 29", clmul(100, 29), 1172)
    same("1172 in binary", bin(1172), "0b10010010100")
    same("ordinary 100 times 29", 100 * 29, 2900)
    same("0x11D shifted left by 2", 0x11D << 2, 1140)
    same("1172 ^ 1140", 1172 ^ 1140, 224)
    same("1172 reduced by 0x11D", reduce(1172), 224)
    same("100 plus 29 in GF(2^8), an XOR", 100 ^ 29, 121)
    # Polynomial addition mod 2 is XOR of the integers that hold the coefficients.
    pairs = 0
    for a in range(256):
        for b in range(256):
            pa = [(a >> i) & 1 for i in range(8)]
            pb = [(b >> i) & 1 for i in range(8)]
            assert poly_add(pa, pb) == [((a ^ b) >> i) & 1 for i in range(8)], (a, b)
            pairs += 1
    same("pairs of degree-7 polynomials checked", pairs, 65536)


# Worked examples and exercises --------------------------------------------------

def check_examples_and_exercises() -> None:
    print("The truncation error of 1/10, worked examples and exercises")
    tenth = Fraction(1, 10)
    for k in range(1, 25):
        e_k = Fraction(2 ** k % 10, 10 * 2 ** k)
        assert e_k == tenth - Fraction(2 ** k // 10, 2 ** k), k
        assert Fraction(1, 5 * 2 ** k) <= e_k <= Fraction(8, 10 * 2 ** k) < Fraction(1, 2 ** k), k
    same("e_k formula checked, k from 1 to 24", len(range(1, 25)), 24)
    same("last digit of 2**k, k = 1 to 8", [2 ** k % 10 for k in range(1, 9)], [2, 4, 8, 6, 2, 4, 8, 6])
    def exact(x: Fraction) -> Decimal:
        """A terminating fraction as an exact decimal, so it prints the way the chapter writes it."""
        return Decimal(x.numerator) / Decimal(x.denominator)

    same("e_8", exact(Fraction(2 ** 8 % 10, 10 * 2 ** 8)), Decimal("0.00234375"))
    same("e_4", exact(Fraction(2 ** 4 % 10, 10 * 2 ** 4)), Decimal("0.0375"))
    same("e_4 as 6/160", exact(Fraction(6, 160)), Decimal("0.0375"))
    same("1/10 cut to four places", exact(Fraction(2 ** 4 // 10, 2 ** 4)), Decimal("0.0625"))
    same("2**56 / 10", exact(Fraction(2 ** 56, 10)), Decimal("7205759403792793.6"))
    same("2**52 <= 2**56/10 < 2**53", 2 ** 52 <= Fraction(2 ** 56, 10) < 2 ** 53, True)
    same("2**-4 <= 0.1 < 2**-3", 2.0 ** -4 <= 0.1 < 2.0 ** -3, True)
    same("ulp(0.3) is 2 ** -54", ulp(0.3) == 2.0 ** -54, True)

    ties = [0.5, 1.5, 2.5, 3.5]
    up = [int(t + 0.5) for t in ties]
    same("the ties rounded half up", up, [1, 2, 3, 4])
    same("their sum", sum(up), 10)

    print("Example: the byte that holds 29")
    same("29 == 16 + 8 + 4 + 1", 29 == 16 + 8 + 4 + 1, True)
    same("format(29, '08b')", format(29, "08b"), "00011101")
    same("29 >> 3", 29 >> 3, 3)
    same("29 // 8", 29 // 8, 3)
    same("29 & 1", 29 & 1, 1)

    print("Example: dividing -17 by 5")
    same("-17 // 5", -17 // 5, -4)
    same("-17 % 5", -17 % 5, 3)
    same("(-4) * 5 + 3", (-4) * 5 + 3, -17)
    same("int(-17 / 5)", int(-17 / 5), -3)
    same("-17 - 5 * int(-17 / 5)", -17 - 5 * int(-17 / 5), -2)

    print("Example: which fractions terminate in binary")
    same("3/8 terminates", terminates(Fraction(3, 8)), True)
    same("5/12 terminates", terminates(Fraction(5, 12)), False)
    same("7/64 terminates", terminates(Fraction(7, 64)), True)
    same("1/5 terminates", terminates(Fraction(1, 5)), False)
    same("6/12 terminates", terminates(Fraction(6, 12)), True)
    same("6/12 in lowest terms", Fraction(6, 12), Fraction(1, 2))

    print("Example: the gap at 3.0 and at 1024.0")
    same("ulp(3.0) is 2 ** -51", ulp(3.0) == 2.0 ** -51, True)
    same("ulp(1024.0) is 2 ** -42", ulp(1024.0) == 2.0 ** -42, True)
    below = 1024.0 - nextafter(1024.0, 0.0)
    same("gap just below 1024.0 is 2 ** -43", below == 2.0 ** -43, True)
    same("gap just below 1024.0 equals ulp(1000.0)", below == ulp(1000.0), True)

    print("Example: reading the pin cocotb~=2.1")
    for v, verdict in [((2, 0, 5), False), ((2, 1, 3), True), ((2, 4, 0), True), ((3, 1, 0), False)]:
        same(f"{'.'.join(map(str, v))} accepted by ~=2.1", compatible(v, "2.1"), verdict)

    print("Exercises")
    same("[(100 >> i) & 1 for i in range(8)]", [(100 >> i) & 1 for i in range(8)], [0, 0, 1, 0, 0, 1, 1, 0])
    same("round(4.5)", round(4.5), 4)
    same("round(5.5)", round(5.5), 6)
    same("round(-0.5) == 0", round(-0.5) == 0, True)
    same("7379/7640 agrees with '0.966'", agrees(7379 / 7640, "0.966"), True)
    check("distance from 0.966", abs(7379 / 7640 - 0.966), "0.000162")
    same("half unit of '0.966'", half_unit("0.966"), 0.5 * 10.0 ** -3)
    sums = [sum(G[i][j] for i in (0, 3)) for j in range(7)]
    same("column sums of rows 0 and 3", sums, [1, 0, 0, 1, 2, 2, 1])
    same("encode([1, 0, 0, 1])", encode([1, 0, 0, 1]), [1, 0, 0, 1, 0, 0, 1])
    same("0x11D", 0x11D, 285)
    same("256 + 16 + 13", 256 + 16 + 13, 285)
    same("f'{100:08b}'", f"{100:08b}", "01100100")
    same("(0, 71) > (0, 9)", version("Yosys 0.71") > version("Yosys 0.9"), True)
    same("'0.71' < '0.9' as strings", "0.71" < "0.9", True)
    same("'3.10' < '3.9' as strings", "3.10" < "3.9", True)

    print("Version table cells")
    for v, spec, verdict in [
        ((2, 2, 0), "2.1", True), ((2, 0, 0), "2.1.0", False), ((2, 1, 5), "2.1.0", True),
        ((2, 9, 3), "2.1.0", False), ((3, 0, 0), "2.1.0", False), ((2, 0, 0), "2.1", False),
        ((3, 0, 0), "2.1", False),
    ]:
        same(f"{'.'.join(map(str, v))} accepted by ~={spec}", compatible(v, spec), verdict)

    print("Bits named in the preview")
    same("bits set in 100", {i for i in range(8) if 100 >> i & 1}, {6, 5, 2})
    same("bits set in 29", {i for i in range(8) if 29 >> i & 1}, {4, 3, 2, 0})
    same("(1172).bit_length()", (1172).bit_length(), 11)


# The bench --------------------------------------------------------------------

# name, command, oldest version the listings need (None where no official
# source states one), the version the book's runs used, the chapter that first
# needs it. A tool first needed in a later Part never fails the bench.
LATER = "a later Part"
TOOLS = [
    ("git", ["git", "--version"], None, "2.56.1", "Project: A Channel Simulator"),
    ("Verilator", ["verilator", "--version"], (5, 36, 0), "5.054", "Lab: A First Hardware Checker"),
    ("C++ compiler", ["c++", "--version"], None, "", "Lab: A First Hardware Checker"),
    ("make", ["make", "--version"], None, "", "Lab: A First Hardware Checker"),
    ("Yosys", ["yosys", "-V"], None, "0.71", "Lab: A First Hardware Checker"),
    ("Surfer", ["surfer", "--version"], None, "0.7.0", "Lab: A First Hardware Checker"),
    ("Icarus Verilog", ["iverilog", "-V"], (11, 0, 0), "13.0", LATER),
    ("SymbiYosys", ["sby", "--version"], None, "0.71", LATER),
    ("Z3", ["z3", "--version"], None, "5.4.1", LATER),
]

PACKAGES = [("cocotb", (2, 1, 0)), ("numpy", None), ("galois", None)]


def bench() -> int:
    from importlib import metadata

    print("\nThe lab bench on this machine")
    status = 0
    py = sys.version_info[:3]
    note = "ok" if py >= (3, 10) else "older than 3.10, too old"
    if (3, 10) <= py < (3, 12):
        note = "ok for the scripts, though NumPy 2.5 needs 3.12, so pip picks an older NumPy"
    if py >= (3, 15):
        note = "too new for cocotb 2.1, which needs 3.14 or older"
    if py < (3, 10) or py >= (3, 15):
        status = 1
    print(f"  {'Python':<16} {'.'.join(map(str, py)):<12} {note}")
    for name, cmd, floor, measured, first in TOOLS:
        path = shutil.which(cmd[0])
        if path is None:
            print(f"  {name:<16} {'missing':<12} first needed in {first}")
            if first != LATER:
                status = 1
            continue
        try:
            done = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            found = version(done.stdout + done.stderr)
        except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
            print(f"  {name:<16} {'unreadable':<12} {exc}")
            if first != LATER:
                status = 1
            continue
        shown = ".".join(map(str, found))
        if floor is None:
            note = f"no stated floor, the book used {measured}" if measured else "found"
        elif found >= floor:
            note = f"ok, at least {'.'.join(map(str, floor))}"
        else:
            note = f"older than {'.'.join(map(str, floor))}, too old"
            if first != LATER:
                status = 1
        print(f"  {name:<16} {shown:<12} {note}")
    for name, floor in PACKAGES:
        try:
            found = version(metadata.version(name))
        except metadata.PackageNotFoundError:
            print(f"  {name:<16} {'missing':<12} install it into the virtual environment")
            status = 1
            continue
        ok = floor is None or found >= floor
        status |= 0 if ok else 1
        print(f"  {name:<16} {'.'.join(map(str, found)):<12} {'ok' if ok else 'too old'}")
    try:
        import numpy as np
        m, g = np.array([1, 0, 1, 1]), np.array(G)
        same("NumPy: (m @ G) % 2", ((m @ g) % 2).tolist(), [1, 0, 1, 1, 0, 1, 0])
    except ImportError:
        pass
    try:
        import galois
        GF = galois.GF(2 ** 8)
        same("galois: GF(100) + GF(29)", int(GF(100) + GF(29)), 121)
        same("galois: GF(100) * GF(29)", int(GF(100) * GF(29)), 224)
        same("galois: the field's polynomial", str(GF.irreducible_poly), "x^8 + x^4 + x^3 + x^2 + 1")
    except ImportError:
        pass
    print("Bench passes." if status == 0 else "Bench needs attention.")
    return status


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bench", action="store_true", help="also inspect the tools installed on this machine")
    args = ap.parse_args()
    check_integers_and_bits()
    check_shift_theorem()
    check_floor_division()
    check_finite_expansion_theorem()
    check_point_one()
    check_spacing_theorem()
    check_rounding()
    check_half_unit_theorem()
    check_versions()
    check_preview()
    check_examples_and_exercises()
    print("\nAll chapter f2 numbers reproduced.")
    if args.bench:
        sys.exit(bench())


if __name__ == "__main__":
    main()
