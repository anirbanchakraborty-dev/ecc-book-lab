#!/usr/bin/env python3
"""Reproduce every number in chapter 23 of Error Correction Coding, "Weight Distributions and Undetected Errors".

Run it from the repository root:

    python3 python/ch23_weight_distributions_and_undetected_errors.py

It needs only the Python 3.10 standard library. Each check prints the value
computed here beside the value the chapter states, and the program stops with
exit status 1 at the first pair that disagrees. Whole numbers, bit patterns,
weight distributions and polynomial coefficients are compared exactly. A
probability is compared to half a unit in the last digit the chapter prints.

The small codes are built from their parity-check matrices by trying every
word of the length, so each weight distribution is counted rather than taken
from a formula. The MacWilliams transform is computed by multiplying out the
polynomials (x + y) and (x - y) with integer coefficients. The longer codes of
the last section, with five parity checks and up to 31 bits, are reached from
their 32 dual words instead, through the same transform run backwards.

A probability curve is sampled at p = j/1000 for j = 0 to 500, and every power
is formed by repeated multiplication, so the Octave twin
octave/ch23_weight_distributions_and_undetected_errors.m prints the same lines.
"""

import sys
from itertools import combinations_with_replacement
from math import comb

GRID = 1000


# Helpers ----------------------------------------------------------------------

def half_unit(stated: str) -> float:
    text = stated.replace(",", "").rstrip("%")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", "").rstrip("%"))
    if stated.endswith("%"):
        computed = 100 * float(computed)
    print(f"  {label:<62} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<62} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<62} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH in '{label}': the claim does not hold")
        sys.exit(1)


def section(title: str) -> None:
    print(f"\n{title}")


def word(v) -> str:
    return "".join(str(x) for x in v)


def listing(v) -> str:
    return ",".join(str(x) for x in v)


def parse(text: str) -> list:
    return [int(ch) for ch in text]


# Codes ------------------------------------------------------------------------

def all_words(n: int) -> list:
    """Every word of length n, in counting order, first bit most significant."""
    return [[(x >> (n - 1 - t)) & 1 for t in range(n)] for x in range(2 ** n)]


def syndrome(H: list, e: list) -> list:
    return [sum(h * b for h, b in zip(row, e)) % 2 for row in H]


def code_from_checks(H: list) -> list:
    n = len(H[0])
    return [w for w in all_words(n) if not any(syndrome(H, w))]


def span(G: list) -> list:
    """Every sum of rows of G, in counting order of the selecting message."""
    k, n = len(G), len(G[0])
    out = []
    for m in all_words(k):
        out.append([sum(m[i] * G[i][t] for i in range(k)) % 2 for t in range(n)])
    return out


def distribution(words: list, n: int) -> list:
    A = [0] * (n + 1)
    for w in words:
        A[sum(w)] += 1
    return A


def dot(u: list, v: list) -> int:
    return sum(a * b for a, b in zip(u, v)) % 2


def polymul(a: list, b: list) -> list:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def term(n: int, i: int) -> list:
    """The coefficients of y^0 to y^n in (x + y)^(n - i) (x - y)^i."""
    out = [1]
    for _ in range(n - i):
        out = polymul(out, [1, 1])
    for _ in range(i):
        out = polymul(out, [1, -1])
    return out


def transform(A: list) -> list:
    """The coefficients of W(x + y, x - y), multiplied out."""
    n = len(A) - 1
    out = [0] * (n + 1)
    for i in range(n + 1):
        if A[i]:
            t = term(n, i)
            for j in range(n + 1):
                out[j] += A[i] * t[j]
    return out


def macwilliams(A: list) -> list:
    size = sum(A)
    raw = transform(A)
    if any(c % size for c in raw):
        print("MISMATCH: W(x + y, x - y) is not divisible by the size of the code")
        sys.exit(1)
    return [c // size for c in raw]


def krawtchouk(n: int, j: int, i: int) -> int:
    return sum((-1) ** s * comb(i, s) * comb(n - i, j - s) for s in range(j + 1) if j - s <= n - i)


def columns_dual(cols: list, r: int) -> list:
    """The weight distribution of the dual code of the code whose check columns are cols."""
    n = len(cols)
    B = [0] * (n + 1)
    for a in range(2 ** r):
        B[sum(bin(a & c).count("1") % 2 for c in cols)] += 1
    return B


def from_dual(B: list, r: int) -> list:
    n = len(B) - 1
    A = []
    for j in range(n + 1):
        total = sum(B[i] * krawtchouk(n, j, i) for i in range(n + 1))
        if total % 2 ** r:
            print("MISMATCH: a Krawtchouk sum is not divisible by the size of the dual code")
            sys.exit(1)
        A.append(total // 2 ** r)
    return A


def gf2_rank(values: list) -> int:
    basis = []
    for v in values:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


# Probabilities ----------------------------------------------------------------

def powers(base: float, n: int) -> list:
    out = [1.0]
    for _ in range(n):
        out.append(out[-1] * base)
    return out


def pud(A: list, p: float) -> float:
    n = len(A) - 1
    pw = powers(p, n)
    qw = powers(1 - p, n)
    total = 0.0
    for i in range(1, n + 1):
        total += A[i] * pw[i] * qw[n - i]
    return total


def pud_dual_parts(B: list, r: int, p: float):
    n = len(B) - 1
    tw = powers(1 - 2 * p, n)
    total = 0.0
    for j in range(n + 1):
        total += B[j] * tw[j]
    return total / 2 ** r, powers(1 - p, n)[n]


def curve(A: list) -> list:
    return [pud(A, j / GRID) for j in range(GRID // 2 + 1)]


def peak(A: list):
    v = curve(A)
    j = max(range(len(v)), key=v.__getitem__)
    return v[j], j


def rising(v: list) -> bool:
    return all(v[j + 1] >= v[j] for j in range(len(v) - 1))


def thousandths(j: int) -> str:
    return f"0.{j:03d}"


# The codes of the chapter -----------------------------------------------------

H3 = [[(c >> (2 - t)) & 1 for c in range(1, 8)] for t in range(3)]
HAMMING = code_from_checks(H3)
REP3 = code_from_checks([parse("110"), parse("011")])
EVEN3 = code_from_checks([parse("111")])
EVEN4 = code_from_checks([parse("1111")])
REP4 = code_from_checks([parse("1100"), parse("0110"), parse("0011")])
SIMPLEX = span(H3)
HEXT = [row + [0] for row in H3] + [[1] * 8]
EXTENDED = code_from_checks(HEXT)
TRIPLE_H = [parse("00000011"), parse("00000101"), parse("11111001")]
TRIPLE = code_from_checks(TRIPLE_H)


def dist(code: list) -> list:
    return distribution(code, len(code[0]))


def shortened(n: int) -> list:
    return [row[:n] for row in H3]


def reach(e: list) -> int:
    """Positions from the first flip of a nonzero pattern to the last, inclusive."""
    ones = [i for i, b in enumerate(e) if b]
    return ones[-1] - ones[0] + 1


# C_5 of Distance, Detection and Correction.
C5 = [parse(w) for w in ("00000", "01101", "10110", "11011")]


# The checks -------------------------------------------------------------------

def check_silent_check() -> None:
    section("When a check stays silent: the (7,4) Hamming code as a detector")
    same("row 1 of H_3", word(H3[0]), "0001111")
    same("row 2 of H_3", word(H3[1]), "0110011")
    same("row 3 of H_3", word(H3[2]), "1010101")
    examples = [("0010000", "011"), ("0100100", "111"), ("0000111", "100"),
                ("1110000", "000"), ("1001100", "000"), ("0001111", "000"),
                ("1111111", "000")]
    for e, s in examples:
        same(f"syndrome of the error pattern {e}", word(syndrome(H3, parse(e))), s)
    patterns = [e for e in all_words(7) if any(e)]
    missed = [e for e in patterns if not any(syndrome(H3, e))]
    same("nonzero error patterns of length 7", len(patterns), 127)
    same("patterns with syndrome 000", len(missed), 15)
    same("patterns caught", len(patterns) - len(missed), 112)
    require("the missed patterns are the nonzero codewords", sorted(map(word, missed)) == sorted(word(c) for c in HAMMING if any(c)))
    A = dist(HAMMING)
    stated_total = {1: 7, 2: 21, 3: 35, 4: 35, 5: 21, 6: 7, 7: 1}
    stated_missed = {1: 0, 2: 0, 3: 7, 4: 7, 5: 0, 6: 0, 7: 1}
    for w in range(1, 8):
        same(f"patterns of weight {w}", comb(7, w), stated_total[w])
        same(f"patterns of weight {w} that are missed", A[w], stated_missed[w])
    same("weight-3 patterns caught", comb(7, 3) - A[3], 28)
    require("every pattern of weight 1 or 2 is caught", A[1] == 0 and A[2] == 0)


def check_distributions() -> None:
    section("Counting codewords by weight")
    same("weight distribution of the repetition code of length 3", listing(dist(REP3)), "1,0,0,1")
    same("weight distribution of the even-weight code of length 3", listing(dist(EVEN3)), "1,0,3,0")
    same("weight distribution of the even-weight code of length 4", listing(dist(EVEN4)), "1,0,6,0,1")
    same("weight distribution of the repetition code of length 4", listing(dist(REP4)), "1,0,0,0,1")
    same("weight distribution of the (7,4) Hamming code", listing(dist(HAMMING)), "1,0,0,7,7,0,0,1")
    same("weight distribution of the (7,3) simplex code", listing(dist(SIMPLEX)), "1,0,0,0,7,0,0,0")
    same("weight distribution of the extended (8,4) Hamming code", listing(dist(EXTENDED)), "1,0,0,0,14,0,0,0,1")
    same("W(1,1) for the (7,4) Hamming code", sum(dist(HAMMING)), 16)
    same("minimum distance of the (7,4) code, first nonzero weight", min(sum(c) for c in HAMMING if any(c)), 3)
    for name, code in (("(7,4) Hamming", HAMMING), ("extended (8,4)", EXTENDED), ("even-weight length 4", EVEN4)):
        A = dist(code)
        n = len(code[0])
        ok = True
        for c in code:
            seen = [0] * (n + 1)
            for d in code:
                seen[sum((a + b) % 2 for a, b in zip(c, d))] += 1
            ok = ok and seen == A
        require(f"distances from every word of the {name} code match A", ok)
    c = parse("0110011")
    require("0110011 is a codeword of the (7,4) code", c in HAMMING)
    seen = [0] * 8
    for d in HAMMING:
        seen[sum((a + b) % 2 for a, b in zip(c, d))] += 1
    same("distances from 0110011, counted by value", listing(seen), "1,0,0,7,7,0,0,1")
    neighbours = sorted(word(d) for d in HAMMING if sum((a + b) % 2 for a, b in zip(c, d)) == 3)
    same("the first codeword at distance 3 from 0110011", neighbours[0], "0010110")
    same("their sum, a codeword of weight 3", word([(a + b) % 2 for a, b in zip(c, parse("0010110"))]), "0100101")
    for name, code in (("repetition length 3", REP3), ("even-weight length 4", EVEN4),
                       ("(7,4) Hamming", HAMMING), ("extended (8,4)", EXTENDED)):
        A = dist(code)
        require(f"A_i = A_(n-i) for the {name} code", A == A[::-1])
    require("the even-weight code of length 3 lacks 111 and is not symmetric", dist(EVEN3) != dist(EVEN3)[::-1])


def check_probability() -> None:
    section("The probability of an undetected error")
    A = dist(HAMMING)
    p = 0.01
    q = 1 - p
    check("7 p^3 (1 - p)^4 at p = 0.01", 7 * powers(p, 3)[3] * powers(q, 4)[4], "6.72417e-6")
    check("7 p^4 (1 - p)^3 at p = 0.01", 7 * powers(p, 4)[4] * powers(q, 3)[3], "6.79209e-8")
    check("p^7 at p = 0.01", powers(p, 7)[7], "1e-14")
    check("P_ud of the (7,4) code at p = 0.01", pud(A, 0.01), "6.79209e-6")
    p = 0.1
    q = 1 - p
    check("7 p^3 (1 - p)^4 at p = 0.1", 7 * powers(p, 3)[3] * powers(q, 4)[4], "0.0045927")
    check("7 p^4 (1 - p)^3 at p = 0.1", 7 * powers(p, 4)[4] * powers(q, 3)[3], "0.0005103")
    check("p^7 at p = 0.1", powers(p, 7)[7], "1e-7")
    check("P_ud of the (7,4) code at p = 0.1", pud(A, 0.1), "0.0051031")
    check("P_ud of the (7,4) code at p = 1/2", pud(A, 0.5), "0.1171875")
    same("P_ud of the (7,4) code at p = 1/2, as a fraction", f"{sum(A) - 1}/{2 ** 7}", "15/128")
    anyerr = 1 - powers(0.99, 7)[7]
    check("probability of any error in 7 bits at p = 0.01", anyerr, "0.0679347")
    check("share of erroneous words that are missed at p = 0.01", pud(A, 0.01) / anyerr, "9.998e-5")
    check("share of erroneous words that are missed at p = 1/2", pud(A, 0.5) / (1 - 0.5 ** 7), "0.11811")
    check("the rule of thumb 2^-3", 2.0 ** -3, "0.125")
    stated = {"repetition length 3": (REP3, "1e-6", "0.001", "0.125", "1/8"),
              "even-weight length 4": (EVEN4, "5.8807e-4", "0.0487", "0.4375", "7/16"),
              "(7,4) Hamming": (HAMMING, "6.79209e-6", "0.0051031", "0.1171875", "15/128"),
              "extended (8,4)": (EXTENDED, "1.34483e-7", "9.1855e-4", "0.05859375", "15/256")}
    for name, (code, s1, s2, s3, frac) in stated.items():
        D = dist(code)
        n = len(D) - 1
        k = (len(code)).bit_length() - 1
        check(f"P_ud at p = 0.01, {name}", pud(D, 0.01), s1)
        check(f"P_ud at p = 0.1, {name}", pud(D, 0.1), s2)
        check(f"P_ud at p = 1/2, {name}", pud(D, 0.5), s3)
        same(f"(2^k - 1)/2^n for the {name} code", f"{2 ** k - 1}/{2 ** n}", frac)
        require(f"P_ud(1/2) equals that fraction for the {name} code", pud(D, 0.5) == (2 ** k - 1) / 2 ** n)
        require(f"P_ud(1/2) = 2^-(n-k) - 2^-n for the {name} code", pud(D, 0.5) == 2.0 ** -(n - k) - 2.0 ** -n)

    section("Short bursts on the (7,4) code")
    bursts = [e for e in all_words(7) if any(e) and reach(e) <= 4]
    same("nonzero patterns whose flips lie within four neighboring positions", len(bursts), 39)
    stated = {1: (7, 0), 2: (15, 0), 3: (13, 2), 4: (4, 2)}
    found = {}
    for w, (count, missed) in stated.items():
        group = [e for e in bursts if sum(e) == w]
        found[w] = (len(group), sum(1 for e in group if not any(syndrome(H3, e))))
        same(f"bursts of weight {w}", found[w][0], count)
        same(f"codewords among them, weight {w}", found[w][1], missed)
    require("no burst has weight 5 or more", all(sum(e) <= 4 for e in bursts))
    require("every burst of weight 4 is a run of four neighboring flips",
            all(reach(e) == 4 for e in bursts if sum(e) == 4))
    for w, s_list in ((3, "0010110 1110000"), (4, "0001111 0111100")):
        hits = sorted(word(e) for e in bursts if sum(e) == w and not any(syndrome(H3, e)))
        same(f"the missed bursts of weight {w}", " ".join(hits), s_list)
    check("share of weight-3 bursts missed, 2/13", found[3][1] / found[3][0], "0.154")
    check("share of all weight-3 patterns missed, 7/35", A[3] / comb(7, 3), "0.2")
    check("share of weight-4 bursts missed, 2/4", found[4][1] / found[4][0], "0.5")
    check("share of all weight-4 patterns missed, 7/35", A[4] / comb(7, 4), "0.2")
    require("weight-3 bursts are missed less often than weight-3 patterns overall",
            found[3][1] / found[3][0] < A[3] / comb(7, 3))
    require("runs of four are missed more than twice as often as weight-4 patterns overall",
            found[4][1] / found[4][0] > 2 * A[4] / comb(7, 4))
    missed_all = sum(m for _, m in found.values())
    same("short bursts that are codewords", missed_all, 4)
    check("share of the 39 short bursts that are missed, 4/39", missed_all / len(bursts), "0.103")
    check("share of all nonzero patterns that are missed, 15/127", (sum(A) - 1) / (2 ** 7 - 1), "0.118")
    require("short bursts are missed less often than nonzero patterns overall",
            missed_all / len(bursts) < (sum(A) - 1) / (2 ** 7 - 1))


def check_character_sum() -> None:
    section("The MacWilliams identity: the character sum and the double sum")
    for name, code, dual in (("repetition length 3", REP3, EVEN3), ("(7,4) Hamming", HAMMING, SIMPLEX)):
        n = len(code[0])
        ok = True
        for u in all_words(n):
            total = sum(1 - 2 * dot(u, c) for c in code)
            ok = ok and total == (len(code) if u in dual else 0)
        require(f"sum of (-1)^(u.c) over the {name} code, every u", ok)
    sums = [sum(1 - 2 * dot(u, c) for c in REP3) for u in all_words(3)]
    same("row sums of the sign grid for the repetition code", listing(sums), "2,0,0,2,0,2,2,0")
    same("words u whose row sum is 2", " ".join(word(u) for u, s in zip(all_words(3), sums) if s == 2), "000 011 101 110")
    for name, code in (("repetition length 3", REP3), ("(7,4) Hamming", HAMMING)):
        n = len(code[0])
        double = [0] * (n + 1)
        factored = True
        for c in code:
            inner = [0] * (n + 1)
            for u in all_words(n):
                inner[sum(u)] += 1 - 2 * dot(u, c)
            factored = factored and inner == term(n, sum(c))
            double = [a + b for a, b in zip(double, inner)]
        require(f"inner sum over u factors per coordinate, {name}", factored)
        dual_dist = macwilliams(dist(code))
        require(f"double sum equals |C| times the dual enumerator, {name}", double == [len(code) * b for b in dual_dist])
        require(f"double sum equals W(x + y, x - y), {name}", double == transform(dist(code)))
    n = 7
    same("K_4(0) for n = 7", krawtchouk(n, 4, 0), 35)
    same("K_4(3) for n = 7", krawtchouk(n, 4, 3), 3)
    same("K_4(4) for n = 7", krawtchouk(n, 4, 4), 3)
    same("K_4(7) for n = 7", krawtchouk(n, 4, 7), 35)
    same("B_4 = (35 + 7*3 + 7*3 + 35)/16", (35 + 7 * 3 + 7 * 3 + 35) // 16, 7)
    A = dist(HAMMING)
    B = [sum(A[i] * krawtchouk(n, j, i) for i in range(n + 1)) // 16 for j in range(n + 1)]
    same("dual distribution of the (7,4) code by Krawtchouk sums", listing(B), "1,0,0,0,7,0,0,0")
    require("K_j(i) is the coefficient of x^(n-j) y^j in (x+y)^(n-i) (x-y)^i, n = 7",
            all(term(7, i)[j] == krawtchouk(7, j, i) for i in range(8) for j in range(8)))


def check_small_identities() -> None:
    section("Checking the identity on small codes")
    t = [a + b for a, b in zip(term(3, 0), term(3, 3))]
    same("(x+y)^3 + (x-y)^3, coefficients of y^0..y^3", listing(t), "2,0,6,0")
    same("divided by |C| = 2", listing(macwilliams(dist(REP3))), "1,0,3,0")
    same("dual distribution read off directly", listing(dist(EVEN3)), "1,0,3,0")
    t = [a + 3 * b for a, b in zip(term(3, 0), term(3, 2))]
    same("(x+y)^3 + 3(x+y)(x-y)^2", listing(t), "4,0,0,4")
    same("divided by |C| = 4", listing(macwilliams(dist(EVEN3))), "1,0,0,1")
    outer = [a + b for a, b in zip(term(7, 0), term(7, 7))]
    same("(x+y)^7 + (x-y)^7", listing(outer), "2,0,42,0,70,0,14,0")
    inner = [7 * a + 7 * b for a, b in zip(term(7, 3), term(7, 4))]
    same("7(x+y)^4(x-y)^3 + 7(x+y)^3(x-y)^4", listing(inner), "14,0,-42,0,42,0,-14,0")
    same("14x(x^2 - y^2)^3, multiplied out", listing([14 * c for c in polymul([1, 0, -1], polymul([1, 0, -1], [1, 0, -1]))] + [0]), "14,0,-42,0,42,0,-14,0")
    same("the whole of W(x+y, x-y) for the (7,4) code", listing([a + b for a, b in zip(outer, inner)]), "16,0,0,0,112,0,0,0")
    same("divided by 16", listing(macwilliams(dist(HAMMING))), "1,0,0,0,7,0,0,0")
    same("the simplex code counted directly", listing(dist(SIMPLEX)), "1,0,0,0,7,0,0,0")
    require("every nonzero simplex word has weight 4", all(sum(w) == 4 for w in SIMPLEX if any(w)))
    same("back from the simplex code, divided by 8", listing(macwilliams(dist(SIMPLEX))), "1,0,0,7,7,0,0,1")
    same("the extended (8,4) code transformed, divided by 16", listing(macwilliams(dist(EXTENDED))), "1,0,0,0,14,0,0,0,1")
    require("the extended (8,4) code equals its own dual", sorted(map(word, EXTENDED)) == sorted(map(word, span(HEXT))))
    same("(x+y)^4 + 6(x+y)^2(x-y)^2 + (x-y)^4", listing(transform(dist(EVEN4))), "8,0,0,0,8")
    same("even-weight length 4 transformed, divided by 8", listing(macwilliams(dist(EVEN4))), "1,0,0,0,1")
    same("the repetition code of length 4 counted directly", listing(dist(REP4)), "1,0,0,0,1")


def check_dual_form() -> None:
    section("Undetected errors from the dual code")
    B = dist(SIMPLEX)
    first, last = pud_dual_parts(B, 3, 0.01)
    check("(1 - 2p)^4 at p = 0.01", powers(0.98, 4)[4], "0.92236816")
    check("(1/8)(1 + 7(0.98)^4)", first, "0.93207214")
    check("(0.99)^7", last, "0.93206535")
    check("their difference, P_ud of the (7,4) code at p = 0.01", first - last, "6.79209e-6")
    first, last = pud_dual_parts(B, 3, 0.1)
    check("(1/8)(1 + 7(0.8)^4)", first, "0.4834")
    check("(0.9)^7", last, "0.4782969")
    check("their difference at p = 0.1", first - last, "0.0051031")
    for j in range(GRID // 2 + 1):
        p = j / GRID
        a, b = pud_dual_parts(B, 3, p)
        if abs((a - b) - pud(dist(HAMMING), p)) > 1e-15:
            require("dual form agrees with the direct sum on the grid", False)
    require("dual form agrees with the direct sum on the grid", True)
    ok = all(7 * powers(1 - j / GRID, 6)[6] - 7 * powers(1 - 2 * (j / GRID), 3)[3] >= 0 for j in range(GRID // 2 + 1))
    require("7(1-p)^6 - 7(1-2p)^3 >= 0 on the grid", ok)
    require("P_ud of the (7,4) code rises along the grid", rising(curve(dist(HAMMING))))
    E = dist(EXTENDED)
    check("P_ud of the extended (8,4) code at 0.01, direct sum", pud(E, 0.01), "1.34483e-7")
    first, last = pud_dual_parts(E, 4, 0.01)
    check("(1/16)(1 + 14(0.98)^4 + (0.98)^8)", first, "0.922744829")
    check("(0.99)^8", last, "0.922744694")
    check("their difference", first - last, "1.34483e-7")
    check("relative gap between the two routes, powers by repeated multiplication", abs((first - last) - pud(E, 0.01)) / pud(E, 0.01), "8e-10")


def check_proper_and_good() -> None:
    section("Proper codes and good codes")
    for name, code, r in (("(7,4) Hamming", HAMMING, 3), ("extended (8,4)", EXTENDED, 4)):
        D = dist(code)
        top, j = peak(D)
        same(f"grid point of the largest P_ud, {name}", thousandths(j), "0.500")
        require(f"largest P_ud on the grid is below 2^-{r}, {name}", top < 2.0 ** -r)
    for name, dual, n, ok in (("(7,4) Hamming", SIMPLEX, 7, True), ("extended (8,4)", EXTENDED, 8, True),
                              ("(6,3) shortened Hamming", span(shortened(6)), 6, True),
                              ("triplicated parity (8,5)", span(TRIPLE_H), 8, False)):
        low = min(sum(w) for w in dual if any(w))
        same(f"smallest nonzero dual weight, {name}", low, {7: 4, 8: 4 if ok else 2, 6: 3}[n])
        same(f"is it at least n/2 = {n / 2:g}? {name}", "yes" if 2 * low >= n else "no", "yes" if ok else "no")
    same("row 1 of the triplicated-parity H", word(TRIPLE_H[0]), "00000011")
    same("row 2 of the triplicated-parity H", word(TRIPLE_H[1]), "00000101")
    same("row 3 of the triplicated-parity H", word(TRIPLE_H[2]), "11111001")
    same("number of codewords", len(TRIPLE), 32)
    require("each codeword is five bits and their parity three times",
            all(w[5] == w[6] == w[7] == sum(w[:5]) % 2 for w in TRIPLE))
    T = dist(TRIPLE)
    same("weight distribution of the triplicated-parity code", listing(T), "1,0,10,0,10,0,10,0,1")
    same("weight distribution of its dual", listing(dist(span(TRIPLE_H))), "1,0,3,0,0,0,3,0,1")
    same("the same dual by MacWilliams", listing(macwilliams(T)), "1,0,3,0,0,0,3,0,1")
    check("P_ud of the triplicated-parity code at p = 0.3", pud(T, 0.3), "0.12897")
    top, j = peak(T)
    check("its largest value on the grid", top, "0.12899")
    same("where the largest value sits", thousandths(j), "0.304")
    check("its value at p = 1/2", pud(T, 0.5), "0.12109375")
    same("as a fraction", f"{len(TRIPLE) - 1}/{2 ** 8}", "31/256")
    check("share of weight-2 patterns that are codewords, 10/28", T[2] / comb(8, 2), "0.357")
    same("number of weight-2 patterns of length 8", comb(8, 2), 28)
    require("P_ud exceeds 1/8 at p = 0.3", pud(T, 0.3) > 0.125)
    stated_count = {4: 1, 5: 3, 6: 7, 7: 14, 8: 24}
    stated_bad = {4: 0, 5: 0, 6: 0, 7: 0, 8: 2}
    found = []
    for n in range(4, 9):
        seen = set()
        for cols in combinations_with_replacement(range(1, 8), n):
            if gf2_rank(list(cols)) < 3:
                continue
            B = columns_dual(list(cols), 3)
            if B[1]:
                continue
            seen.add(tuple(from_dual(B, 3)))
        bad = []
        for A in sorted(seen):
            top, j = peak(list(A))
            if top > 0.125:
                bad.append((A, top, j))
        same(f"length {n}: distinct weight distributions", len(seen), stated_count[n])
        same(f"length {n}: distributions rising above 1/8 on the grid", len(bad), stated_bad[n])
        found = bad
    same("first not-good distribution at length 8", listing(found[0][0]), "1,0,10,0,10,0,10,0,1")
    same("second not-good distribution at length 8", listing(found[1][0]), "1,0,10,1,5,10,0,5,0")
    check("its largest value on the grid", found[1][1], "0.12925")
    same("where it sits", thousandths(found[1][2]), "0.302")


def check_lengths() -> None:
    section("Why the length matters: the (7,4) code shortened to n = 4, 5, 6, 7")
    stated = {4: ("1,0,0,1,0", "9.9e-7", "1/16", "0.0625"),
              5: ("1,0,0,2,1,0", "1.9701e-6", "3/32", "0.09375"),
              6: ("1,0,0,4,3,0,0", "3.9106e-6", "7/64", "0.109375"),
              7: ("1,0,0,7,7,0,0,1", "6.79209e-6", "15/128", "0.1171875")}
    for n, (s_dist, s_small, s_frac, s_half) in stated.items():
        code = code_from_checks(shortened(n))
        D = dist(code)
        same(f"weight distribution at n = {n}", listing(D), s_dist)
        check(f"P_ud at p = 0.01, n = {n}", pud(D, 0.01), s_small)
        same(f"P_ud at p = 1/2 as a fraction, n = {n}", f"{len(code) - 1}/{2 ** n}", s_frac)
        check(f"P_ud at p = 1/2, n = {n}", pud(D, 0.5), s_half)
        require(f"P_ud rises along the grid, n = {n}", rising(curve(D)))
    same("dual distribution at n = 5", listing(dist(span(shortened(5)))), "1,0,2,4,1,0")
    same("dual distribution at n = 6", listing(dist(span(shortened(6)))), "1,0,0,4,3,0,0")
    same("dual distribution at n = 4", listing(dist(span(shortened(4)))), "1,1,3,3,0")
    require("at n = 4 the fourth bit is 0 in every codeword", all(w[3] == 0 for w in code_from_checks(shortened(4))))
    five = code_from_checks(shortened(5))
    same("weight distribution of C_5 of Distance, Detection and Correction", listing(dist(C5)), "1,0,0,2,1,0")
    moved = sorted(word([w[i - 1] for i in (3, 1, 4, 2, 5)]) for w in C5)
    same("C_5 with its positions read in the order 3, 1, 4, 2, 5", " ".join(moved), "00000 01111 10011 11100")
    same("the n = 5 code", " ".join(sorted(word(w) for w in five)), "00000 01111 10011 11100")
    same("nonzero codewords that C_5 shares with the n = 5 code", sum(1 for w in C5 if any(w) and w in five), 0)
    check("P_ud of C_5 at p = 0.01", pud(dist(C5), 0.01), "1.9701e-6")
    grid = [j / GRID for j in range(GRID // 2 + 1)]
    require("p^2 (3 - 4p) >= 0 on the grid", all(p * p * (3 - 4 * p) >= 0 for p in grid))
    require("p^2 (6(1-p)^2 - p^2) >= 0 on the grid", all(p * p * (6 * (1 - p) ** 2 - p * p) >= 0 for p in grid))
    require("(1-p)^4 - (1-2p)^2 >= 0 on the grid", all((1 - p) ** 4 - (1 - 2 * p) ** 2 >= 0 for p in grid))
    check("6(1-p)^2 = p^2 where p = sqrt(6)/(1 + sqrt(6))", 6 ** 0.5 / (1 + 6 ** 0.5), "0.71")

    section("Why the length matters: single parity at even and odd length")
    stated = {8: (28, "0.0026368", "0.0026361", "0.99974"),
              9: (36, "0.0033566", "0.0033554", "0.99964"),
              10: (45, "0.0041543", "0.0041524", "0.99952")}
    values = {}
    for n, (pairs, s_total, s_pairs, s_share) in stated.items():
        code = code_from_checks([[1] * n])
        D = dist(code)
        total = pud(D, 0.01)
        w2 = D[2] * powers(0.01, 2)[2] * powers(0.99, n - 2)[n - 2]
        values[n] = total
        same(f"weight-2 codewords at n = {n}", D[2], pairs)
        check(f"P_ud at p = 0.01, n = {n}", total, s_total)
        check(f"the weight-2 term, n = {n}", w2, s_pairs)
        check(f"its share of P_ud, n = {n}", w2 / total, s_share)
        same(f"all-ones word is a codeword at n = {n}", "yes" if D[n] else "no", "yes" if n % 2 == 0 else "no")
        require(f"no odd-weight pattern is a codeword, n = {n}", all(D[i] == 0 for i in range(1, n + 1, 2)))
    check("P_ud at n = 9 over P_ud at n = 8", values[9] / values[8], "1.273")
    check("36/28 times 0.99", 36 / 28 * 0.99, "1.273")

    section("Why the length matters: five parity checks, n = 16 to 31")
    cols = list(range(1, 17))
    require("n = 16 has a coordinate that is 0 in every codeword", columns_dual(cols, 5)[1] > 0)
    require("n = 15 has only four independent checks", gf2_rank(list(range(1, 16))) == 4)
    stated = {17: (36, 112, "0.290", "no"), 18: (38, 126, "0.309", "no"),
              19: (41, 148, "0.500", "yes"), 20: (45, 176, "0.500", "yes"),
              21: (50, 213, "0.500", "yes"), 22: (56, 259, "0.500", "yes"),
              23: (63, 315, "0.500", "yes"), 24: (71, 371, "0.500", "yes"),
              25: (80, 438, "0.500", "yes"), 26: (90, 516, "0.500", "yes"),
              27: (101, 606, "0.500", "yes"), 28: (113, 706, "0.500", "yes"),
              29: (126, 819, "0.500", "yes"), 30: (140, 945, "0.500", "yes"),
              31: (155, 1085, "0.500", "yes")}
    ratio = {17: "1.0961", 18: "1.0188"}
    shares = []
    for n, (a3, a4, s_at, s_good) in stated.items():
        B = columns_dual(list(range(1, n + 1)), 5)
        A = from_dual(B, 5)
        shares.append(A[3] / comb(n, 3))
        top, j = peak(A)
        same(f"A_3 at n = {n}", A[3], a3)
        same(f"A_4 at n = {n}", A[4], a4)
        same(f"where the largest P_ud on the grid sits, n = {n}", thousandths(j), s_at)
        same(f"stays at or below 2^-5 on the grid, n = {n}", "yes" if top <= 2.0 ** -5 else "no", s_good)
        if n in ratio:
            check(f"2^5 times the largest P_ud on the grid, n = {n}", 32 * top, ratio[n])
        else:
            require(f"the largest value is 2^-5 - 2^-n, at p = 1/2, n = {n}", top == pud(A, 0.5) == 2.0 ** -5 - 2.0 ** -n)
    require("A_3 / C(n, 3) is above 2^-5 at every n from 17 to 31", all(s > 2.0 ** -5 for s in shares))
    A = from_dual(columns_dual(list(range(1, 18)), 5), 5)
    check("largest P_ud on the grid at n = 17", peak(A)[0], "0.034254")
    check("2^-5", 2.0 ** -5, "0.03125")
    same("weight-3 patterns of length 17", comb(17, 3), 680)
    check("share of them that are codewords, 36/680", A[3] / comb(17, 3), "0.0529")
    same("weight-3 codewords inside the first 15 positions", from_dual(columns_dual(list(range(1, 16)), 4), 4)[3], 35)
    check("their share of the weight-3 patterns there, 35/455", 35 / comb(15, 3), "0.0769")
    A = from_dual(columns_dual(list(range(1, 32)), 5), 5)
    same("A_3 at n = 31 equals the Hamming count n(n-1)/6", A[3], 31 * 30 // 6)
    same("weight-3 patterns of length 31", comb(31, 3), 4495)
    check("share of them that are codewords, 155/4495", A[3] / comb(31, 3), "0.0345")


def check_examples() -> None:
    section("Worked examples")
    D = dist(code_from_checks(shortened(6)))
    check("(6,3) code: 4 p^3 (1-p)^3 at p = 0.1", 4 * powers(0.1, 3)[3] * powers(0.9, 3)[3], "0.002916")
    check("(6,3) code: 3 p^4 (1-p)^2 at p = 0.1", 3 * powers(0.1, 4)[4] * powers(0.9, 2)[2], "0.000243")
    check("(6,3) code: P_ud at p = 0.1", pud(D, 0.1), "0.003159")
    check("extended (8,4): 14 p^4 (1-p)^4 at p = 0.01", 14 * powers(0.01, 4)[4] * powers(0.99, 4)[4], "1.34483e-7")
    check("extended (8,4): p^8 at p = 0.01", powers(0.01, 8)[8], "1e-16")
    first, last = pud_dual_parts(dist(EXTENDED), 4, 0.01)
    a, b = f"{first:.9f}", f"{last:.9f}"
    agree = 0
    while a[2 + agree] == b[2 + agree]:
        agree += 1
    same(f"leading digits shared by {a} and {b}", agree, 6)


def check_exercises() -> None:
    section("Exercises")
    code = code_from_checks(shortened(5))
    same("(5,2) code: its four codewords", " ".join(word(w) for w in code), "00000 01111 10011 11100")
    same("(5,2) code: dual distribution by MacWilliams", listing(macwilliams(dist(code))), "1,0,2,4,1,0")
    check("(5,2) code: P_ud at p = 1/2", pud(dist(code), 0.5), "0.09375")
    D = dist(code_from_checks(shortened(6)))
    same("(6,3) code: nonzero codewords", sum(D) - 1, 7)
    check("(6,3) code: P_ud at p = 1/2, 4/64 + 3/64", pud(D, 0.5), "0.109375")
    first, last = pud_dual_parts(dist(EVEN3), 2, 0.1)
    check("repetition length 3 at p = 0.1, direct p^3", pud(dist(REP3), 0.1), "0.001")
    check("repetition length 3 at p = 0.1, (1/4)(1 + 3(0.8)^2)", first, "0.73")
    check("repetition length 3 at p = 0.1, (0.9)^3", last, "0.729")
    check("repetition length 3 at p = 0.1, by the dual", first - last, "0.001")
    grid = [j / GRID for j in range(GRID // 2 + 1)]
    require("12p(1-p)(1-2p) + 4p^3 >= 0 on the grid", all(12 * p * (1 - p) * (1 - 2 * p) + 4 * p ** 3 >= 0 for p in grid))
    require("P_ud of the even-weight code of length 4 rises along the grid", rising(curve(dist(EVEN4))))
    check("even-weight length 4: P_ud(1/2) = 7/16 against 2^-3", pud(dist(EVEN4), 0.5) / 0.125, "3.5")
    same("single parity, n = 10: weight-2 codewords", comb(10, 2), 45)
    check("single parity, n = 10: P_ud at p = 0.01", pud(dist(code_from_checks([[1] * 10])), 0.01), "0.0041543")
    require("an inverted bus of 10 bits is missed, of 9 bits caught",
            dist(code_from_checks([[1] * 10]))[10] == 1 and dist(code_from_checks([[1] * 9]))[9] == 0)


def main() -> None:
    print("Chapter 23: Weight Distributions and Undetected Errors")
    check_silent_check()
    check_distributions()
    check_probability()
    check_character_sum()
    check_small_identities()
    check_dual_form()
    check_proper_and_good()
    check_lengths()
    check_examples()
    check_exercises()
    print("\nAll chapter 23 numbers reproduced.")


if __name__ == "__main__":
    main()
