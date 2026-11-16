"""Reproduce every number in chapter 27 of Error Correction Coding, "Bounds on Codes".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch27_bounds_on_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
The chapter's numbers are whole numbers, ratios printed to a fixed number of
places, or lists of whole numbers, so every comparison is exact once the
value is printed in the chapter's own format.

A binary word of length n is held as an integer whose bit j is position j of
the word. A linear code is listed in full by XOR-ing every subset of the rows
of its generator matrix, and its minimum distance is read as its smallest
nonzero weight. The cyclic codes of length 15 have generator matrices built
from the shifts g(x), x g(x), ..., x^(k-1) g(x) of their generator
polynomials. The module also holds the functions that ch27_charts.py draws,
so the figures and this check share one source.
"""

import functools
import sys
from fractions import Fraction
from math import comb, floor, log2


def same(label: str, computed, stated) -> None:
    print(f"  {label:<62} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<62} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH: '{label}' does not hold")
        sys.exit(1)


def section(title: str) -> None:
    print(f"\n{title}")


def join(values) -> str:
    return " ".join(str(v) for v in values)


# ---------------------------------------------------------------- words and codes

def weight(x: int) -> int:
    return bin(x).count("1")


def span(rows):
    """Every codeword of the binary linear code whose generator rows are given."""
    words = [0]
    for row in rows:
        words = words + [w ^ row for w in words]
    return words


def min_weight(words) -> int:
    return min(weight(w) for w in words if w)


def min_pair_distance(words) -> int:
    return min(weight(a ^ b) for i, a in enumerate(words) for b in words[i + 1:])


def poly_rows(g: int, n: int):
    """Rows g(x), x g(x), ... of a cyclic code of length n, g given as a bit mask."""
    deg = g.bit_length() - 1
    return [g << i for i in range(n - deg)]


def parse(bits: str) -> int:
    """A word written left to right as position 0, 1, 2, ..."""
    return sum(1 << j for j, b in enumerate(bits) if b == "1")


def show(x: int, n: int) -> str:
    return "".join("1" if (x >> j) & 1 else "0" for j in range(n))


# ---------------------------------------------------------------- the bounds

def volume(n: int, t: int, q: int = 2) -> int:
    return sum(comb(n, i) * (q - 1) ** i for i in range(t + 1))


def hamming_bound(n: int, d: int) -> int:
    """Largest M the sphere-packing bound allows for binary (n, M, d) codes."""
    return 2 ** n // volume(n, (d - 1) // 2)


def largest_k(m: int) -> int:
    """Largest k with 2^k <= m."""
    return m.bit_length() - 1


def hamming_k(n: int, d: int) -> int:
    return largest_k(hamming_bound(n, d))


def singleton_k(n: int, d: int) -> int:
    return n - d + 1


def plotkin_simple(n: int, d: int) -> int:
    """Floor of 2d / (2d - n), valid when 2d > n."""
    return (2 * d) // (2 * d - n)


def plotkin_sharp(n: int, d: int) -> int:
    """2 floor(d / (2d - n)), valid when 2d > n."""
    return 2 * (d // (2 * d - n))


def plotkin_best(n: int, d: int) -> int:
    """The sharpened Plotkin bound, through the extended length when d is odd."""
    best = plotkin_sharp(n, d)
    if d % 2 == 1 and 2 * (d + 1) > n + 1:
        best = min(best, plotkin_sharp(n + 1, d + 1))
    return best


def griesmer_sum(k: int, d: int, q: int = 2) -> int:
    return sum(-(-d // q ** i) for i in range(k))


def griesmer_k(n: int, d: int) -> int:
    k = 1
    while griesmer_sum(k + 1, d) <= n:
        k += 1
    return k


def gilbert_value(n: int, d: int) -> float:
    return 2 ** n / volume(n, d - 1)


def varshamov_k(n: int, d: int) -> int:
    """Largest k with V(n-1, d-2) < 2^(n-k), at least 1 when d <= n."""
    v = volume(n - 1, d - 2)
    k = n - 1
    while k > 1 and not v < 2 ** (n - k):
        k -= 1
    return k


def forbidden_sums(columns, d: int) -> set:
    """Every sum of at most d - 2 of the given columns, the empty sum 0 included."""
    sums = {0}
    for _ in range(d - 2):
        sums |= {s ^ c for s in sums for c in columns}
    return sums


def varshamov_greedy(r: int, d: int):
    """Choose r-bit columns one at a time, each the smallest one that no sum of
    at most d - 2 earlier columns equals. Returns the forbidden set seen at each
    step and the columns chosen, stopping at the first step with nothing left."""
    chosen, steps = [], []
    while True:
        forbidden = forbidden_sums(chosen, d)
        steps.append(sorted(forbidden))
        free = [c for c in range(2 ** r) if c not in forbidden]
        if not free:
            return steps, chosen
        chosen.append(free[0])


def columns_code(columns, r: int):
    """The code whose parity-check matrix has the given r-bit columns, every word listed."""
    n = len(columns)
    return [x for x in range(2 ** n)
            if not functools.reduce(lambda a, j: a ^ columns[j], (j for j in range(n) if x >> j & 1), 0)]


def elias_bound(n: int, d: int):
    """Smallest value of (nd/2) / (r^2 - nr + nd/2) * 2^n / V(n, r) over 1 <= r <= n/2."""
    best = None
    for r in range(1, n // 2 + 1):
        den = Fraction(r * r) - n * r + Fraction(n * d, 2)
        if den > 0:
            value = Fraction(n * d, 2) / den * Fraction(2 ** n, volume(n, r))
            if best is None or value < best[0]:
                best = (value, r)
    return best


def elias_floor(n: int, d: int) -> int:
    return floor(elias_bound(n, d)[0])


def elias_k(n: int, d: int) -> int:
    return largest_k(elias_floor(n, d))


# ---------------------------------------------------------------- the codes of length 15

N15 = 15
G_HAMMING = 0b10011                      # x^4 + x + 1
G_EXP_HAMMING = 0b110101                 # (x + 1)(x^4 + x + 1)
G_BCH2 = parse("100010111")              # 1 + x^4 + x^6 + x^7 + x^8
G_BCH2_EVEN = None                       # filled below as (x + 1) g(x)
G_BCH3 = parse("11101100101")            # 1 + x + x^2 + x^4 + x^5 + x^8 + x^10


def times_x_plus_1(g: int) -> int:
    return g ^ (g << 1)


G_BCH2_EVEN = times_x_plus_1(G_BCH2)


def simplex_rows(m: int):
    n = 2 ** m - 1
    return [sum(1 << j for j in range(n) if ((j + 1) >> i) & 1) for i in range(m)]


def code_15_2_10():
    return [sum(1 << j for j in range(10)), sum(1 << j for j in range(5, 15))]


def codes_15():
    """(label, k, d) of every code of length 15 the chapter names, each enumerated."""
    table = [
        ("even-weight code", [(1 << j) | (1 << 14) for j in range(14)]),
        ("Hamming code, x^4+x+1", poly_rows(G_HAMMING, N15)),
        ("even-weight subcode, (x+1)(x^4+x+1)", poly_rows(G_EXP_HAMMING, N15)),
        ("degree-8 generator", poly_rows(G_BCH2, N15)),
        ("(x+1) times the degree-8 generator", poly_rows(G_BCH2_EVEN, N15)),
        ("degree-10 generator", poly_rows(G_BCH3, N15)),
        ("simplex code", simplex_rows(4)),
        ("two rows of weight 10", code_15_2_10()),
        ("repetition code", [(1 << 15) - 1]),
    ]
    out = []
    for label, rows in table:
        words = span(rows)
        k = largest_k(len(words))
        out.append((label, k, min_weight(words), len(words)))
    return out


# The best code size the chapter exhibits at each d of length 15 (largest k).
EXHIBITED_K = {2: 14, 3: 11, 4: 10, 5: 7, 6: 6, 7: 5, 8: 4, 9: 2, 10: 2,
               11: 1, 12: 1, 13: 1, 14: 1, 15: 1}

D_RANGE = list(range(2, 16))

# A code of 20 words of length 8 at distance 3, found once by a branch-and-bound
# clique search over the 256 words with the all-zero word fixed.
CODE_8_20_3 = [0, 7, 25, 30, 42, 53, 75, 84, 97, 108, 114, 127,
               140, 147, 166, 169, 176, 194, 197, 216]

CODE_5_2_3 = [parse(w) for w in ("00000", "10101", "01111", "11010")]
ROWS_6_3_3 = [parse(w) for w in ("100011", "010101", "001110")]
ROWS_7_4_3 = [parse(w) for w in ("1000110", "0100101", "0010011", "0001111")]
ROWS_8_4_4 = [parse(w) for w in ("10001101", "01001011", "00100111", "00011110")]


# ---------------------------------------------------------------- hardware

def sec_bits(k: int) -> int:
    r = 1
    while 2 ** r < k + r + 1:
        r += 1
    return r


def dec_bits(k: int) -> int:
    r = 1
    while 2 ** r < 1 + (k + r) + comb(k + r, 2):
        r += 1
    return r


DATA_WIDTHS = [4, 8, 16, 32, 64, 128, 256]


# ---------------------------------------------------------------- the checks

def check_easy_values() -> None:
    section("Easy values of A_2(n,d)")
    even = span([(1 << j) | (1 << 14) for j in range(14)])
    same("even-weight code of length 15, number of words", len(even), 16384)
    same("even-weight code of length 15, minimum distance", min_weight(even), 2)
    same("2^(n-1) at n = 15", 2 ** 14, 16384)
    rep = span([(1 << 15) - 1])
    same("repetition code of length 15, number of words", len(rep), 2)
    same("repetition code of length 15, minimum distance", min_weight(rep), 15)


def check_small_exact() -> None:
    section("Small exact values")
    same("(5,3): simple Plotkin 2d/(2d-n)", plotkin_simple(5, 3), 6)
    same("(5,3): sharpened Plotkin on A_2(6,4), 2 floor(4/2)", plotkin_sharp(6, 4), 4)
    same("(5,3): words in the exhibited code", len(CODE_5_2_3), 4)
    same("(5,3): its minimum distance", min_pair_distance(CODE_5_2_3), 3)
    require("(5,3): the exhibited code is closed under XOR",
            all(a ^ b in CODE_5_2_3 for a in CODE_5_2_3 for b in CODE_5_2_3))
    same("(6,3): sharpened Plotkin on A_2(7,4), 2 floor(4/1)", plotkin_sharp(7, 4), 8)
    c633 = span(ROWS_6_3_3)
    same("(6,3): words in the shortened Hamming code", len(c633), 8)
    same("(6,3): its minimum distance", min_weight(c633), 3)
    same("(7,3): Hamming bound 2^7/(1+7)", hamming_bound(7, 3), 16)
    c743 = span(ROWS_7_4_3)
    same("(7,3): words in the Hamming code", len(c743), 16)
    same("(7,3): its minimum distance", min_weight(c743), 3)
    same("exact values A_2(5,3), A_2(6,3), A_2(7,3)",
         join([len(CODE_5_2_3), len(c633), len(c743)]), "4 8 16")


def check_nonlinear() -> None:
    section("A nonlinear code beats every linear one at (8,3)")
    same("words in the listed (8,3) code", len(set(CODE_8_20_3)), 20)
    pairs = sum(1 for i in range(20) for _ in range(i + 1, 20))
    same("pairs checked", pairs, 190)
    same("smallest distance over all pairs", min_pair_distance(CODE_8_20_3), 3)
    require("every listed word has 8 bits", all(0 <= w < 256 for w in CODE_8_20_3))
    ext = [(w << 1) | (weight(w) & 1) for w in CODE_8_20_3]
    same("extended by parity to length 9: words, distance",
         join([len(set(ext)), min_pair_distance(ext)]), "20 4")
    same("Hamming bound for (8,3), floor(256/9)", hamming_bound(8, 3), 28)
    same("[8,5,3]: 2^5 (1 + 8)", 2 ** 5 * (1 + 8), 288)
    require("[8,5,3]: 288 > 2^8 = 256, so no such code", 2 ** 5 * 9 > 2 ** 8)
    c844 = span(ROWS_8_4_4)
    same("extended Hamming code: words", len(c844), 16)
    same("extended Hamming code: minimum distance", min_weight(c844), 4)
    same("largest linear code at (8,3), words", 2 ** hamming_k(8, 3), 16)


def check_hamming() -> None:
    section("Hamming bound at n = 15")
    same("V_2(15,t) for t = 0..7", join(volume(15, t) for t in range(8)),
         "1 16 121 576 1941 4944 9949 16384")
    same("largest k, d = 3..12",
         join(hamming_k(15, d) for d in range(3, 13)), "11 11 8 8 5 5 4 4 2 2")
    same("largest k, d = 13..15", join(hamming_k(15, d) for d in range(13, 16)), "1 1 1")
    same("Hamming bound on M at d = 5", hamming_bound(15, 5), 270)
    same("Hamming bound on M at d = 9", hamming_bound(15, 9), 16)
    same("Hamming bound on M at d = 11", hamming_bound(15, 11), 6)
    same("sharpened Plotkin on M at d = 11", plotkin_best(15, 11), 2)


def check_singleton() -> None:
    section("Singleton bound")
    same("largest k at n = 15, d = 2..15",
         join(singleton_k(15, d) for d in D_RANGE), "14 13 12 11 10 9 8 7 6 5 4 3 2 1")
    for label, k, d in (("repetition code (15,1,15)", 1, 15),
                        ("even-weight code (15,14,2)", 14, 2)):
        require(f"{label} meets k = n - d + 1", k == 15 - d + 1)
    require("Hamming (7,4,3) misses it: 4 < 7 - 3 + 1 = 5", 4 < 5)
    same("[5,3,3]: 2^3 V_2(5,1) against 2^5", join([2 ** 3 * volume(5, 1), 2 ** 5]), "48 32")


def check_plotkin() -> None:
    section("Plotkin bound")
    rows = CODE_5_2_3
    cols = [sum((w >> j) & 1 for w in rows) for j in range(5)]
    same("(5,3) code: ones per column", join(cols), "2 2 2 2 2")
    contrib = [2 * w * (4 - w) for w in cols]
    same("(5,3) code: 2w(M-w) per column", join(contrib), "8 8 8 8 8")
    ordered = sum(weight(a ^ b) for a in rows for b in rows)
    same("(5,3) code: sum of distances over ordered pairs", ordered, 40)
    same("(5,3) code: column total", sum(contrib), 40)
    same("lower bound M(M-1)d", 4 * 3 * 3, 36)
    same("upper bound n M^2 / 2", 5 * 16 // 2, 40)
    same("simple form at n = 15, d = 8, 9, 10",
         join(plotkin_simple(15, d) for d in (8, 9, 10)), "16 6 4")
    same("sharpened, d = 9 through A_2(16,10)", plotkin_sharp(16, 10), 4)
    same("sharpened, d = 8 and d = 10",
         join(plotkin_sharp(15, d) for d in (8, 10)), "16 4")
    same("sharpened bound at n = 15, d = 8..15",
         join(plotkin_best(15, d) for d in range(8, 16)), "16 4 4 2 2 2 2 2")
    simplex = span(simplex_rows(4))
    same("simplex code: words", len(simplex), 16)
    same("simplex code: minimum distance", min_weight(simplex), 8)
    require("simplex code: every nonzero weight is 8",
            all(weight(w) == 8 for w in simplex if w))
    for d in (2, 3, 4, 5, 6, 7, 8):
        require(f"A_2(2d,d) <= 4d at d = {d}: shortened bound times 2",
                2 * plotkin_sharp(2 * d - 1, d) == 4 * d)


def check_plotkin_extension() -> None:
    section("Exercise: Plotkin for n >= 2d at n = 15")
    ext = []
    vals = []
    for d in range(2, 8):
        dd, nn = (d, 15) if d % 2 == 0 else (d + 1, 16)
        vals.append(dd * 2 ** (nn - 2 * dd + 2))
        ext.append(largest_k(vals[-1]))
    same("bound on words, d = 2..7", join(vals), "16384 4096 2048 384 192 32")
    same("largest k from d 2^(n-2d+2), d = 2..7", join(ext), "14 12 11 8 7 5")
    same("Griesmer largest k, d = 2..7", join(griesmer_k(15, d) for d in range(2, 8)),
         "14 12 11 8 7 5")


def check_griesmer() -> None:
    section("Griesmer bound")
    stated = [((7, 4, 3), 7), ((8, 4, 4), 8), ((15, 4, 8), 15), ((16, 5, 8), 16),
              ((32, 6, 16), 32), ((15, 5, 7), 15), ((15, 11, 3), 14), ((15, 7, 5), 14),
              ((23, 12, 7), 22), ((24, 12, 8), 23)]
    for (n, k, d), s in stated:
        same(f"sum for ({n},{k},{d})", griesmer_sum(k, d), s)
    meets = [f"({n},{k},{d})" for (n, k, d), s in stated if griesmer_sum(k, d) == n]
    same("codes that meet it with equality", " ".join(meets),
         "(7,4,3) (8,4,4) (15,4,8) (16,5,8) (32,6,16) (15,5,7)")
    same("(15,5,7) sum written out", "+".join(str(-(-7 // 2 ** i)) for i in range(5)),
         "7+4+2+1+1")
    same("largest k at n = 15, d = 2..10",
         join(griesmer_k(15, d) for d in range(2, 11)), "14 12 11 8 7 5 4 2 2")
    below = [d for d in D_RANGE if griesmer_k(15, d) < hamming_k(15, d)]
    equal = [d for d in D_RANGE if griesmer_k(15, d) == hamming_k(15, d)]
    above = [d for d in D_RANGE if griesmer_k(15, d) > hamming_k(15, d)]
    same("n = 15: Griesmer below Hamming at d", join(below), "2 6 8 9 10 11 12")
    same("n = 15: Griesmer equal to Hamming at d", join(equal), "4 5 7 13 14 15")
    same("n = 15: Griesmer above Hamming at d", join(above), "3")
    same("(15,2,10): sum 10 + 5", griesmer_sum(2, 10), 15)
    same("Golay (23,12,7): volume V_2(23,3)", volume(23, 3), 2048)
    require("Golay (23,12,7) is perfect: 2^12 * 2048 = 2^23", 2 ** 12 * volume(23, 3) == 2 ** 23)


def check_gv() -> None:
    section("Gilbert-Varshamov bound")
    same("V_2(7,2)", volume(7, 2), 29)
    same("Gilbert at (7,3), 128/29", f"{gilbert_value(7, 3):.2f}", "4.41")
    same("so at least this many words", -(-128 // 29), 5)
    same("Varshamov at d = 3, n - k = 3: largest n with n < 8", max(n for n in range(1, 20) if n < 2 ** 3), 7)
    same("V_2(14, d-2) for d = 2..15", join(volume(14, d - 2) for d in D_RANGE),
         "1 15 106 470 1471 3473 6476 9908 12911 14913 15914 16278 16369 16383")
    same("Varshamov largest k at n = 15, d = 2..8",
         join(varshamov_k(15, d) for d in range(2, 9)), "14 11 8 6 4 3 2")
    same("Varshamov largest k at n = 15, d = 9..15",
         join(varshamov_k(15, d) for d in range(9, 16)), "1 1 1 1 1 1 1")
    same("Gilbert log2 at n = 15, d = 2..8",
         " ".join(f"{log2(gilbert_value(15, d)):.3f}" for d in range(2, 9)),
         "11.000 8.081 5.830 4.077 2.729 1.720 1.000")
    steps, cols = varshamov_greedy(3, 3)
    same("greedy columns at d = 3, r = 3", join(cols), "1 2 3 4 5 6 7")
    same("forbidden count at steps 1..8", join(len(f) for f in steps), "1 2 3 4 5 6 7 8")
    words = columns_code(cols, 3)
    same("greedy code at d = 3, r = 3: (n, k, d)",
         f"({len(cols)},{largest_k(len(words))},{min_weight(words)})", "(7,4,3)")
    same("Varshamov at d = 4, r = 4: largest n with V(n-1,2) < 16",
         max(n for n in range(2, 20) if volume(n - 1, 2) < 2 ** 4), 5)
    steps, cols = varshamov_greedy(4, 4)
    same("greedy columns at d = 4, r = 4", join(cols), "1 2 4 7 8 11 13 14")
    words = columns_code(cols, 4)
    same("greedy code at d = 4, r = 4: (n, k, d)",
         f"({len(cols)},{largest_k(len(words))},{min_weight(words)})", "(8,4,4)")


def check_elias() -> None:
    section("Elias bound at n = 15")
    same("Elias check at (10,3), floor", elias_floor(10, 3), 232)
    same("(10,3): best r and V_2(10,1)", join([elias_bound(10, 3)[1], volume(10, 1)]), "1 11")
    same("(15,4): best r and r^2 - 15r + 30 there",
         join([elias_bound(15, 4)[1], 2 ** 2 - 15 * 2 + 30]), "2 4")
    same("floor of the Elias bound, d = 2..8",
         join(elias_floor(15, d) for d in range(2, 9)), "30720 5421 2031 883 284 104 30")
    same("floor of the Elias bound, d = 9..15",
         join(elias_floor(15, d) for d in range(9, 16)), "11 7 6 5 4 4 3")
    same("best r, d = 2..15", join(elias_bound(15, d)[1] for d in D_RANGE),
         "1 1 2 2 3 4 7 7 7 7 7 7 7 7")
    same("largest k, d = 2..15", join(elias_k(15, d) for d in D_RANGE),
         "14 12 10 9 8 6 4 3 2 2 2 2 2 1")
    require("Elias never below an exhibited code at n = 15",
            all(elias_floor(15, d) >= 2 ** EXHIBITED_K[d] for d in D_RANGE))
    require("Elias at d = 4 rules out [15,11,4]: 2031 < 2048", elias_floor(15, 4) < 2 ** 11)


def check_codes_15() -> None:
    section("Real codes at n = 15, every codeword enumerated")
    stated = {
        "even-weight code": (14, 2, 16384),
        "Hamming code, x^4+x+1": (11, 3, 2048),
        "even-weight subcode, (x+1)(x^4+x+1)": (10, 4, 1024),
        "degree-8 generator": (7, 5, 128),
        "(x+1) times the degree-8 generator": (6, 6, 64),
        "degree-10 generator": (5, 7, 32),
        "simplex code": (4, 8, 16),
        "two rows of weight 10": (2, 10, 4),
        "repetition code": (1, 15, 2),
    }
    for label, k, d, size in codes_15():
        same(f"{label}: (k, d)", f"({k},{d})", f"({stated[label][0]},{stated[label][1]})")
        require(f"{label}: 2^k distinct words", size == 2 ** k == stated[label][2])
    for d, k in EXHIBITED_K.items():
        require(f"exhibited k = {k} at d = {d} within every upper bound",
                k <= min(singleton_k(15, d), hamming_k(15, d), griesmer_k(15, d), elias_k(15, d))
                and (2 * d <= 15 or 2 ** k <= plotkin_best(15, d)))
        require(f"exhibited k = {k} at d = {d} at least Varshamov's k",
                k >= varshamov_k(15, d))
    same("[15,11,4] punctured to [14,11,3]: 2^(14-11)", 2 ** 3, 8)
    require("8 < 1 + 14 = 15, so no [14,11,3] and no [15,11,4]", 2 ** 3 < 15)
    same("Hamming at d = 3: 2^11 * 16", 2 ** 11 * 16, 2 ** 15)
    same("simplex: Griesmer sum", griesmer_sum(4, 8), 15)
    same("simplex: Plotkin value at (15,8)", plotkin_sharp(15, 8), 16)
    same("(15,5,7): 32768/576", f"{32768 / 576:.1f}", "56.9")
    same("(15,5,7): largest k from Hamming", hamming_k(15, 7), 5)
    same("(15,5,7): Griesmer sum", griesmer_sum(5, 7), 15)
    same("Hamming largest k at d = 2", hamming_k(15, 2), 15)
    same("Griesmer largest k, d = 11..15", join(griesmer_k(15, d) for d in range(11, 16)), "1 1 1 1 1")
    same("Plotkin largest k, d = 8..15",
         join(largest_k(plotkin_best(15, d)) for d in range(8, 16)), "4 2 2 1 1 1 1 1")
    upper = []
    for d in D_RANGE:
        u = min(singleton_k(15, d), hamming_k(15, d), griesmer_k(15, d), elias_k(15, d))
        if 2 * d > 15:
            u = min(u, largest_k(plotkin_best(15, d)))
        upper.append(u)
    same("smallest upper bound on k, d = 2..15", join(upper), "14 11 10 8 7 5 4 2 2 1 1 1 1 1")
    same("exhibited k, d = 2..15", join(EXHIBITED_K[d] for d in D_RANGE), "14 11 10 7 6 5 4 2 2 1 1 1 1 1")
    same("distances where no bound applied directly meets a code",
         join(d for d in D_RANGE if upper[d - 2] > EXHIBITED_K[d]), "5 6")
    same("A_2(15,4) = A_2(14,3) <= floor(2^14 / V_2(14,1))", 2 ** 14 // volume(14, 1), 1092)
    require("1092 < 2048, and below the Elias value 2031",
            2 ** 14 // volume(14, 1) < 2 ** 11 and 2 ** 14 // volume(14, 1) < elias_floor(15, 4))


def check_residual() -> None:
    section("Residual codes close d = 5 and d = 6 at n = 15")
    same("[15,8,5] residual: length, dimension, distance", join([15 - 5, 8 - 1, -(-5 // 2)]),
         "10 7 3")
    same("[10,7,3]: 2^7 * V_2(10,1) against 2^10", join([2 ** 7 * volume(10, 1), 2 ** 10]),
         "1408 1024")
    same("[15,7,6] residual: length, dimension, distance", join([15 - 6, 7 - 1, -(-6 // 2)]),
         "9 6 3")
    same("[9,6,3]: 2^6 * V_2(9,1) against 2^9", join([2 ** 6 * volume(9, 1), 2 ** 9]),
         "640 512")
    require("2^7 (1 + m) > 2^m for m = 1..10, so no [15,8,5]",
            all(2 ** 7 * (1 + m) > 2 ** m for m in range(1, 11)))
    require("2^6 (1 + m) > 2^m for m = 1..9, so no [15,7,6]",
            all(2 ** 6 * (1 + m) > 2 ** m for m in range(1, 10)))
    require("2^m / (1 + m) increases for m = 1..10",
            all(2 ** m * (2 + m) < 2 ** (m + 1) * (1 + m) for m in range(1, 10)))
    same("Griesmer sums for k = 8, d = 5 and k = 7, d = 6",
         join([griesmer_sum(8, 5), griesmer_sum(7, 6)]), "15 15")


def check_hardware() -> None:
    section("Check bits for 64 data bits")
    same("SEC check bits at k = 64", sec_bits(64), 7)
    same("SECDED check bits at k = 64", sec_bits(64) + 1, 8)
    same("DEC check bits at k = 64, sphere-packing floor", dec_bits(64), 12)
    same("DEC at r = 12: 1 + 76 + C(76,2)", 1 + 76 + comb(76, 2), 2927)
    same("DEC-TED floor at k = 64: DEC floor plus 1", dec_bits(64) + 1, 13)
    same("DEC check bits at k = 32", dec_bits(32), 10)
    same("DEC at k = 32, r = 10: 1 + 42 + C(42,2) against 2^10",
         join([1 + 42 + comb(42, 2), 2 ** 10]), "904 1024")
    same("DEC at k = 32, r = 9: 1 + 41 + C(41,2) against 2^9",
         join([1 + 41 + comb(41, 2), 2 ** 9]), "862 512")
    same("SEC at k = 4..256", join(sec_bits(k) for k in DATA_WIDTHS), "3 4 5 6 7 8 9")
    same("SECDED at k = 4..256", join(sec_bits(k) + 1 for k in DATA_WIDTHS), "4 5 6 7 8 9 10")
    same("DEC floor at k = 4..256", join(dec_bits(k) for k in DATA_WIDTHS), "6 7 9 10 12 14 16")
    same("SECDED with 7 checks, punctured to [70,64,3]: 2^6", 2 ** 6, 64)
    require("64 < 1 + 70, so 7 check bits cannot give SECDED", 2 ** 6 < 71)
    same("(8,3) nonlinear code: log2 20 to three places", f"{log2(20):.3f}", "4.322")


def check_examples() -> None:
    section("Worked examples and exercises")
    same("(5,3): Hamming bound floor(32/6)", hamming_bound(5, 3), 5)
    dists = sorted(weight(a ^ b) for i, a in enumerate(CODE_5_2_3) for b in CODE_5_2_3[i + 1:])
    same("(5,3) code: the six pairwise distances", join(dists), "3 3 3 3 4 4")
    same("Griesmer sum for k = 6, d = 7", griesmer_sum(6, 7), 16)
    same("(15,11): 32768/4944 to two places", f"{32768 / 4944:.2f}", "6.63")
    same("128 data bits: SEC, SECDED, DEC floor",
         join([sec_bits(128), sec_bits(128) + 1, dec_bits(128)]), "8 9 14")
    same("SEC at k = 128, r = 7: 2^7 against 128 + 7 + 1", join([2 ** 7, 128 + 7 + 1]), "128 136")
    same("DEC at k = 128, r = 14: 1 + 142 + C(142,2) against 2^14",
         join([1 + 142 + comb(142, 2), 2 ** 14]), "10154 16384")
    same("DEC at k = 128, r = 13: 1 + 141 + C(141,2) against 2^13",
         join([1 + 141 + comb(141, 2), 2 ** 13]), "10012 8192")
    same("Hamming bound at (10,3), floor(1024/11)", hamming_bound(10, 3), 93)
    same("Singleton largest k at (5,3)", singleton_k(5, 3), 3)
    same("[15,3,10]: Griesmer sum 10 + 5 + 3", griesmer_sum(3, 10), 18)
    same("[15,2,10]: Plotkin value at (15,10)", plotkin_best(15, 10), 4)


def main() -> None:
    check_easy_values()
    check_small_exact()
    check_nonlinear()
    check_hamming()
    check_singleton()
    check_plotkin()
    check_plotkin_extension()
    check_griesmer()
    check_gv()
    check_elias()
    check_codes_15()
    check_residual()
    check_hardware()
    check_examples()
    print("\nAll chapter 27 numbers reproduced.")


if __name__ == "__main__":
    main()
