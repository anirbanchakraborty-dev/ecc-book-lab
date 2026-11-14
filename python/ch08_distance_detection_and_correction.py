"""Reproduce every number in chapter 8 of Error Correction Coding, "Distance, Detection and Correction".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch08_distance_detection_and_correction.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. Small codes are
checked exhaustively: every codeword against every error pattern the chapter
names, every combination of errors and erasures inside the bound, all 2^16
rectangular parity arrays, and every reflected Gray code up to 8 bits. The
check-digit schemes are tested at every position against every digit pair, on
numbers built to be valid, and the Luhn and ISBN results are confirmed again on
every valid number of five digits. A run that reaches the final line has
reproduced every number in the chapter.

The module also holds the encoders, decoders and checkers that ch08_charts.py
plots, so the charts and this check share one source.
"""

import random
import sys
from fractions import Fraction
from itertools import combinations, product
from math import comb, log2


# Helpers ------------------------------------------------------------------------

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
        print(f"MISMATCH: '{label}' does not hold")
        sys.exit(1)


def section(title: str) -> None:
    print(f"\n{title}")


# Words, distance and codes ----------------------------------------------------

def weight(x: int) -> int:
    return bin(x).count("1")


def distance(x: int, y: int) -> int:
    return weight(x ^ y)


def word(s: str) -> int:
    return int(s, 2)


def show(x: int, n: int) -> str:
    return format(x, f"0{n}b")


def code_of(*words: str) -> list:
    return [word(s) for s in words]


def min_distance(code: list) -> int:
    return min(distance(a, b) for a, b in combinations(code, 2))


def distance_distribution(code: list) -> dict:
    counts = {}
    for a, b in combinations(code, 2):
        k = distance(a, b)
        counts[k] = counts.get(k, 0) + 1
    return dict(sorted(counts.items()))


def closed_under_xor(code: list) -> bool:
    members = set(code)
    return all((a ^ b) in members for a in code for b in code)


def ball_size(n: int, t: int) -> int:
    return sum(comb(n, i) for i in range(t + 1))


def ball(center: int, n: int, t: int) -> set:
    return {y for y in range(1 << n) if distance(center, y) <= t}


def cube_edges(n: int) -> list:
    return [(x, x ^ (1 << i)) for x in range(1 << n) for i in range(n) if x < x ^ (1 << i)]


def nearest(code: list, r: int) -> list:
    best = min(distance(c, r) for c in code)
    return [c for c in code if distance(c, r) == best]


def decode_nearest(code: list, r: int):
    """Nearest-codeword decoding. Returns None on a tie."""
    winners = nearest(code, r)
    return winners[0] if len(winners) == 1 else None


def decode_bounded(code: list, r: int, t: int):
    """Correct when exactly one codeword lies within t, otherwise flag (None)."""
    inside = [c for c in code if distance(c, r) <= t]
    return inside[0] if len(inside) == 1 else None


def decode_erasures(code: list, n: int, received: str):
    """Nearest codeword on the positions that are not erased ('?'). None on a tie."""
    kept = [i for i, ch in enumerate(received) if ch != "?"]
    def gap(c):
        bits = show(c, n)
        return sum(bits[i] != received[i] for i in kept)
    best = min(gap(c) for c in code)
    winners = [c for c in code if gap(c) == best]
    return winners[0] if len(winners) == 1 else None


def apply_errors_erasures(c: int, n: int, flips, erased) -> str:
    bits = list(show(c, n))
    for i in flips:
        bits[i] = "1" if bits[i] == "0" else "0"
    for i in erased:
        bits[i] = "?"
    return "".join(bits)


def random_code(seed: int, words: int, n: int) -> list:
    rng = random.Random(seed)
    return [rng.getrandbits(n) for _ in range(words)]


def codeword_errors(code: list, n: int, p: float) -> list:
    """Exact error probability of each codeword under ML decoding on BSC(p), ties split evenly."""
    pw = [p**k * (1 - p) ** (n - k) for k in range(n + 1)]
    correct = [0.0] * len(code)
    for y in range(1 << n):
        dists = [distance(c, y) for c in code]
        best = min(dists)
        winners = [i for i, dd in enumerate(dists) if dd == best]
        for i in winners:
            correct[i] += pw[best] / len(winners)
    return [1 - v for v in correct]


def ml_block_error(code: list, n: int, p: float) -> float:
    errs = codeword_errors(code, n, p)
    return sum(errs) / len(errs)


def undetected_probability(code: list, n: int, p: float) -> float:
    """For a code closed under XOR: the chance the error pattern is a nonzero codeword."""
    return sum(p ** weight(c) * (1 - p) ** (n - weight(c)) for c in code if c)


# Single parity -------------------------------------------------------------------

def parity_alarm(n: int, p):
    return (1 - (1 - 2 * p) ** n) / 2


def any_error(n: int, p):
    return 1 - (1 - p) ** n


def parity_coverage(n: int, p):
    return parity_alarm(n, p) / any_error(n, p)


def parity_alarm_by_enumeration(n: int, p) -> Fraction:
    return sum(p ** weight(e) * (1 - p) ** (n - weight(e)) for e in range(1 << n) if weight(e) % 2 == 1)


# Check digits ---------------------------------------------------------------------

def digits(s: str) -> list:
    return [10 if ch in "Xx" else int(ch) for ch in s if ch not in "- "]


def isbn10_sum(ds: list) -> int:
    return sum((10 - i) * d for i, d in enumerate(ds))


def isbn10_check(payload: list) -> int:
    return (-isbn10_sum(payload)) % 11


def isbn10_valid(ds: list) -> bool:
    return len(ds) == 10 and isbn10_sum(ds) % 11 == 0


def isbn13_weights(length: int = 13) -> list:
    return [1 if i % 2 == 0 else 3 for i in range(length)]


def isbn13_sum(ds: list) -> int:
    return sum(w * d for w, d in zip(isbn13_weights(len(ds)), ds))


def isbn13_check(payload: list) -> int:
    return (-isbn13_sum(payload)) % 10


def isbn13_valid(ds: list) -> bool:
    return len(ds) == 13 and isbn13_sum(ds) % 10 == 0


def luhn_double(x: int) -> int:
    return 2 * x if 2 * x <= 9 else 2 * x - 9


def luhn_terms(ds: list) -> list:
    """The value each digit adds, counting positions from the right, where position 1 is the check digit."""
    n = len(ds)
    return [luhn_double(d) if (n - i) % 2 == 0 else d for i, d in enumerate(ds)]


def luhn_sum(ds: list) -> int:
    return sum(luhn_terms(ds))


def luhn_check(payload: list) -> int:
    return next(c for c in range(10) if luhn_sum(payload + [c]) % 10 == 0)


def luhn_valid(ds: list) -> bool:
    return luhn_sum(ds) % 10 == 0


def make_valid(valid, base: list, fixed: dict):
    """A number equal to base except at the fixed positions and two free positions, that passes valid()."""
    free = [i for i in range(len(base)) if i not in fixed]
    for i, j in combinations(free, 2):
        for a, b in product(range(10), repeat=2):
            ds = list(base)
            for k, v in fixed.items():
                ds[k] = v
            ds[i], ds[j] = a, b
            if valid(ds):
                return ds
    raise ValueError("no valid number")


def swap_escapes(valid, base: list, i: int, j: int, a: int, b: int) -> bool:
    ds = make_valid(valid, base, {i: a, j: b})
    ds[i], ds[j] = ds[j], ds[i]
    return valid(ds)


def substitution_escapes(valid, base: list, i: int, a: int, b: int) -> bool:
    ds = make_valid(valid, base, {i: a})
    ds[i] = b
    return valid(ds)


def scheme_rates(valid, base: list) -> dict:
    """Fractions of single substitutions, adjacent swaps and jump swaps (two apart) that are caught."""
    n = len(base)
    pairs = [(a, b) for a in range(10) for b in range(10) if a != b]
    def rate(cases):
        caught = sum(not escaped for escaped in cases)
        return Fraction(caught, len(cases))
    subs = [substitution_escapes(valid, base, i, a, b) for i in range(n) for a, b in pairs]
    adj = [swap_escapes(valid, base, i, i + 1, a, b) for i in range(n - 1) for a, b in pairs]
    jump = [swap_escapes(valid, base, i, i + 2, a, b) for i in range(n - 2) for a, b in pairs]
    return {"substitution": rate(subs), "adjacent": rate(adj), "jump": rate(jump)}


def escaping_adjacent_pairs(valid, base: list) -> set:
    return {(a, b) for i in range(len(base) - 1) for a in range(10) for b in range(10)
            if a != b and swap_escapes(valid, base, i, i + 1, a, b)}


# Rectangular parity ----------------------------------------------------------------

def rect_encode(rows: list) -> list:
    """Append a parity bit to each row, then a parity row, so every row and column of the result is even."""
    with_rows = [r + [sum(r) % 2] for r in rows]
    return with_rows + [[sum(col) % 2 for col in zip(*with_rows)]]


def rect_failures(array: list) -> tuple:
    bad_rows = [i for i, r in enumerate(array) if sum(r) % 2]
    bad_cols = [j for j, col in enumerate(zip(*array)) if sum(col) % 2]
    return bad_rows, bad_cols


def rect_decode(array: list):
    """Correct one error, flag anything else that breaks a parity. Returns (array or None, verdict)."""
    bad_rows, bad_cols = rect_failures(array)
    if not bad_rows and not bad_cols:
        return array, "clean"
    if len(bad_rows) == 1 and len(bad_cols) == 1:
        fixed = [list(r) for r in array]
        fixed[bad_rows[0]][bad_cols[0]] ^= 1
        return fixed, "corrected"
    return None, "flagged"


def flatten(array: list) -> int:
    value = 0
    for bit in (b for r in array for b in r):
        value = (value << 1) | bit
    return value


# Gray codes ------------------------------------------------------------------------

def gray(b: int) -> int:
    return b ^ (b >> 1)


def gray_inverse(g: int) -> int:
    b = 0
    while g:
        b ^= g
        g >>= 1
    return b


def reflected(n: int) -> list:
    if n == 1:
        return ["0", "1"]
    prev = reflected(n - 1)
    return ["0" + s for s in prev] + ["1" + s for s in reversed(prev)]


# Checks ----------------------------------------------------------------------------

C5 = code_of("00000", "01101", "10110", "11011")
C6 = code_of("000000", "011011", "101101", "110110")
E6 = code_of("000000", "001111", "110011", "111100")
R3 = code_of("000", "111")
R4 = code_of("0000", "1111")
R5 = code_of("00000", "11111")


def check_distance_and_weight() -> None:
    section("Distance and Weight")
    x, y = word("10110"), word("11011")
    same("10110 XOR 11011", show(x ^ y, 5), "01101")
    same("weight of 01101", weight(x ^ y), 3)
    same("d(10110, 11011)", distance(x, y), 3)
    same("weight of 10110", weight(x), 3)
    same("weight of 11011", weight(y), 4)
    same("10110 AND 11011", show(x & y, 5), "10010")
    same("w(x) + w(y) - 2 w(x AND y)", weight(x) + weight(y) - 2 * weight(x & y), 3)
    require("w(x XOR y) = w(x) + w(y) - 2 w(x AND y) for all 5-bit pairs",
            all(weight(a ^ b) == weight(a) + weight(b) - 2 * weight(a & b) for a in range(32) for b in range(32)))
    require("triangle inequality holds for all 6-bit triples",
            all(distance(a, c) <= distance(a, b) + distance(b, c)
                for a in range(64) for b in range(64) for c in range(64)))
    a, b, c = word("00000"), word("01101"), word("11011")
    same("d(00000, 11011)", distance(a, c), 4)
    same("d(00000, 01101) + d(01101, 11011)", distance(a, b) + distance(b, c), 6)
    table = {(show(p, 5), show(q, 5)): distance(p, q) for p, q in combinations(C5, 2)}
    same("C5 distance 00000-01101", table[("00000", "01101")], 3)
    same("C5 distance 00000-10110", table[("00000", "10110")], 3)
    same("C5 distance 00000-11011", table[("00000", "11011")], 4)
    same("C5 distance 01101-10110", table[("01101", "10110")], 4)
    same("C5 distance 01101-11011", table[("01101", "11011")], 3)
    same("C5 distance 10110-11011", table[("10110", "11011")], 3)
    same("pairs in C5", len(table), 6)
    same("minimum distance of C5", min_distance(C5), 3)
    require("C5 is closed under XOR", closed_under_xor(C5))
    same("01101 XOR 10110", show(word("01101") ^ word("10110"), 5), "11011")
    same("smallest nonzero weight in C5", min(weight(c) for c in C5 if c), 3)


def check_minimum_distance() -> None:
    section("Minimum Distance")
    same("d of all 2-bit words", min_distance(list(range(4))), 1)
    same("d of 3-bit even-weight words", min_distance([v for v in range(8) if weight(v) % 2 == 0]), 2)
    same("d of R3", min_distance(R3), 3)
    code = random_code(1, 16, 12)
    same("random code, distinct codewords", len(set(code)), 16)
    same("random code, first codeword", show(code[0], 12), "001000100110")
    same("random code, last codeword", show(code[-1], 12), "001101011011")
    same("random code, pairs", comb(16, 2), 120)
    dist = distance_distribution(code)
    same("random code, distance distribution", dist,
         {2: 2, 3: 6, 4: 15, 5: 23, 6: 17, 7: 29, 8: 22, 9: 5, 10: 1})
    same("random code, minimum distance", min_distance(code), 2)
    close = sorted(tuple(sorted((show(a, 12), show(b, 12)))) for a, b in combinations(code, 2) if distance(a, b) == 2)
    same("random code, pairs at distance 2", close,
         [("110000101100", "110000111000"), ("110010011110", "110011010110")])
    mean = Fraction(sum(k * v for k, v in dist.items()), 120)
    check("random code, mean pairwise distance", mean, "6.11")
    check("random pair, expected distance n/2", Fraction(12, 2), "6")
    check("P(two random 12-bit words within distance 2)", Fraction(1 + 12 + 66, 4096), "0.0193")
    same("words within distance 2 of a 12-bit word, 1 + 12 + 66", ball_size(12, 2), 79)
    check("expected pairs within distance 2 among 120", Fraction(120 * 79, 4096), "2.31")
    expected = {k: 120 * comb(12, k) / 4096 for k in range(13)}
    check("expected pairs at distance 6", expected[6], "27.1")
    check("expected pairs at distance 7", expected[7], "23.2")
    check("ML block error, BSC(0.1)", ml_block_error(code, 12, 0.1), "0.10539")
    errs = codeword_errors(code, 12, 0.01)
    block = sum(errs) / 16
    check("ML block error, BSC(0.01)", block, "0.0032306")
    in_pairs = {word(s) for pair in close for s in pair}
    share = sum(e for c, e in zip(code, errs) if c in in_pairs) / 16
    check("share of the four close-pair codewords at p = 0.01", share / block, "0.82")
    check("first-order estimate p/4 at p = 0.01", 0.01 / 4, "0.0025")
    check("uncoded 4 bits, BSC(0.1)", 1 - 0.9**4, "0.3439")


def check_cube() -> None:
    section("The Geometry of the Binary Cube")
    same("vertices of the 3-cube", 1 << 3, 8)
    same("edges of the 3-cube", len(cube_edges(3)), 12)
    same("edges of the 4-cube", len(cube_edges(4)), 32)
    require("edges of the n-cube = n 2^(n-1), n = 1..8",
            all(len(cube_edges(n)) == n * 2 ** (n - 1) for n in range(1, 9)))
    require("graph distance on the 4-cube equals Hamming distance", cube_graph_distance_matches(4))
    b0, b1 = ball(word("000"), 3, 1), ball(word("111"), 3, 1)
    same("B(000, 1)", sorted(show(v, 3) for v in b0), ["000", "001", "010", "100"])
    same("B(111, 1)", sorted(show(v, 3) for v in b1), ["011", "101", "110", "111"])
    require("the two radius-1 balls of R3 are disjoint and cover the 3-cube", not (b0 & b1) and len(b0 | b1) == 8)
    same("V(3, 1)", ball_size(3, 1), 4)
    same("V(4, 1)", ball_size(4, 1), 5)
    c0, c1 = ball(word("0000"), 4, 1), ball(word("1111"), 4, 1)
    same("words in the two radius-1 balls of R4", len(c0 | c1), 10)
    outside = sorted(show(v, 4) for v in range(16) if v not in c0 | c1)
    same("words of the 4-cube in neither ball", outside, ["0011", "0101", "0110", "1001", "1010", "1100"])
    require("each of them is at distance 2 from both codewords",
            all(distance(word(s), 0) == 2 and distance(word(s), 15) == 2 for s in outside))
    same("4-cube words by weight", [sum(weight(v) == k for v in range(16)) for k in range(5)], [1, 4, 6, 4, 1])
    same("V(5, 1)", ball_size(5, 1), 6)
    same("C5: M V(5, 1)", 4 * ball_size(5, 1), 24)
    stray = [v for v in range(32) if min(distance(v, c) for c in C5) >= 2]
    same("5-bit words in no radius-1 ball of C5", len(stray), 8)
    require("each lies at distance 2 from exactly two codewords of C5",
            all(sorted(distance(v, c) for c in C5) == [2, 2, 3, 3] for v in stray))
    require("V(n, t) counts the ball for n <= 6, t <= n",
            all(len(ball(0, n, t)) == ball_size(n, t) for n in range(1, 7) for t in range(n + 1)))
    same("V(6, 1)", ball_size(6, 1), 7)
    same("16 V(6, 1)", 16 * ball_size(6, 1), 112)
    same("2^6", 2**6, 64)
    same("V(7, 1)", ball_size(7, 1), 8)
    same("16 V(7, 1)", 16 * ball_size(7, 1), 128)
    hamming = hamming_7_4()
    same("Hamming (7,4) codewords", len(set(hamming)), 16)
    same("Hamming (7,4) minimum distance", min_distance(hamming), 3)
    require("Hamming (7,4) radius-1 balls tile the 7-cube",
            len(set().union(*(ball(c, 7, 1) for c in hamming))) == 128)


def cube_graph_distance_matches(n: int) -> bool:
    adjacency = {v: [v ^ (1 << i) for i in range(n)] for v in range(1 << n)}
    for start in range(1 << n):
        seen, frontier, depth = {start: 0}, [start], 0
        while frontier:
            depth += 1
            frontier = [u for v in frontier for u in adjacency[v] if u not in seen]
            for u in frontier:
                seen.setdefault(u, depth)
        if any(seen[v] != distance(start, v) for v in range(1 << n)):
            return False
    return True


def hamming_7_4() -> list:
    rows = [word(r) for r in ("1000110", "0100101", "0010011", "0001111")]
    out = []
    for m in range(16):
        c = 0
        for i in range(4):
            if m >> (3 - i) & 1:
                c ^= rows[i]
        out.append(c)
    return out


def check_detection() -> None:
    section("Detecting Errors")
    sent, pattern = word("01101"), word("00110")
    r = sent ^ pattern
    same("01101 with errors 00110", show(r, 5), "01011")
    same("distances from 01011 to C5", [distance(r, c) for c in C5], [3, 2, 4, 1])
    require("01011 is not a codeword", r not in C5)
    same("01101 with errors 10110", show(sent ^ word("10110"), 5), "11011")
    require("C5 flags every pattern of weight 1 or 2 from every codeword",
            all((c ^ e) not in C5 for c in C5 for e in range(1, 32) if weight(e) <= 2))
    require("detection bound: d - 1 detected, some weight-d pattern missed, on C5, C6, R4, R5",
            all(detects_exactly(code, n) for code, n in ((C5, 5), (C6, 6), (R4, 4), (R5, 5))))
    same("weight distribution of C5", [sum(weight(c) == k for c in C5) for k in range(6)], [1, 0, 0, 2, 1, 0])
    check("C5 undetected-error probability, p = 0.01", undetected_probability(C5, 5, 0.01), "1.9701e-6")
    check("2 p^3 (1-p)^2 + p^4 (1-p)", 2 * 0.01**3 * 0.99**2 + 0.01**4 * 0.99, "1.9701e-6")
    check("C5 chance of any error, p = 0.01", any_error(5, 0.01), "0.049010")
    same("error patterns of weight 3 on 5 bits", comb(5, 3), 10)
    same("weight-3 patterns that escape, from each codeword of C5",
         sorted({sum((c ^ e) in C5 for e in range(32) if weight(e) == 3) for c in C5}), [2])


def detects_exactly(code: list, n: int) -> bool:
    d = min_distance(code)
    caught = all((c ^ e) not in code for c in code for e in range(1, 1 << n) if weight(e) < d)
    missed = any((c ^ e) in code for c in code for e in range(1, 1 << n) if weight(e) == d)
    return caught and missed


def check_correction() -> None:
    section("Correcting Errors")
    r = word("10111")
    same("distances from 10111 to C5", [distance(r, c) for c in C5], [4, 3, 1, 2])
    same("10111 decodes to", show(decode_nearest(C5, r), 5), "10110")
    r = word("10000")
    same("distances from 10000 to C5", [distance(r, c) for c in C5], [1, 4, 2, 3])
    same("10000 decodes to", show(decode_nearest(C5, r), 5), "00000")
    same("t for d = 1..7", [(d - 1) // 2 for d in range(1, 8)], [0, 0, 1, 1, 2, 2, 3])
    same("d - 1 for d = 1..7", [d - 1 for d in range(1, 8)], [0, 1, 2, 3, 4, 5, 6])
    require("C5 corrects every single error from every codeword",
            all(decode_nearest(C5, c ^ (1 << i)) == c for c in C5 for i in range(5)))
    require("R5 corrects every pattern of weight 2 or less",
            all(decode_nearest(R5, c ^ e) == c for c in R5 for e in range(32) if weight(e) <= 2))
    require("correction bound is tight: some weight-(t+1) pattern fails, C5 and R5",
            any(decode_nearest(C5, c ^ e) != c for c in C5 for e in range(32) if weight(e) == 2)
            and any(decode_nearest(R5, c ^ e) != c for c in R5 for e in range(32) if weight(e) == 3))
    same("R4: 0001 with t = 1", show(decode_bounded(R4, word("0001"), 1), 4), "0000")
    same("R4: 0011 with t = 1", decode_bounded(R4, word("0011"), 1), None)
    require("R4 and C6 correct 1 and flag 2 from every codeword",
            all(correct_and_detect(code, n, 1, 2) for code, n in ((R4, 4), (C6, 6))))
    require("R5 corrects 1 and flags 3 (t + s = 4 = d - 1)", correct_and_detect(R5, 5, 1, 3))
    require("R5 with t = 2 cannot also flag 3", not correct_and_detect(R5, 5, 2, 3))


def correct_and_detect(code: list, n: int, t: int, s: int) -> bool:
    for c in code:
        for e in range(1, 1 << n):
            out = decode_bounded(code, c ^ e, t)
            if weight(e) <= t and out != c:
                return False
            if t < weight(e) <= s and out is not None:
                return False
    return True


def check_erasures() -> None:
    section("Errors and Erasures Together")
    same("minimum distance of E6", min_distance(E6), 4)
    received = apply_errors_erasures(word("110011"), 6, flips=[5], erased=[2])
    same("110011, position 3 erased, position 6 flipped", received, "11?010")
    kept = [0, 1, 3, 4, 5]
    gaps = [sum(show(c, 6)[i] != received[i] for i in kept) for c in E6]
    same("distances on the five kept positions", gaps, [3, 4, 1, 2])
    same("11?010 decodes to", show(decode_erasures(E6, 6, received), 6), "110011")
    pairs_d4 = [(t, e) for t in range(3) for e in range(5) if 2 * t + e < 4]
    same("(t, e) allowed by d = 4", pairs_d4, [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1)])
    pairs_d5 = [(t, e) for t in range(3) for e in range(6) if 2 * t + e < 5]
    same("(t, e) allowed by d = 5", pairs_d5, [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (2, 0)])
    require("E6 decodes every combination with 2t + e < 4", erasure_bound_holds(E6, 6))
    require("R5 decodes every combination with 2t + e < 5", erasure_bound_holds(R5, 5))
    require("C5 decodes every combination with 2t + e < 3", erasure_bound_holds(C5, 5))
    require("some combination with 2t + e = d fails, for E6, R5 and C5",
            all(erasure_bound_tight(code, n) for code, n in ((E6, 6), (R5, 5), (C5, 5))))
    require("one erasure in an even-parity word is always filled (d = 2, t = 0, e = 1)",
            erasure_bound_holds([v for v in range(16) if weight(v) % 2 == 0], 4))


def erasure_bound_holds(code: list, n: int) -> bool:
    d = min_distance(code)
    for c in code:
        for e_count in range(n + 1):
            for t in range(n + 1):
                if 2 * t + e_count >= d:
                    continue
                for erased in combinations(range(n), e_count):
                    rest = [i for i in range(n) if i not in erased]
                    for flips in combinations(rest, t):
                        if decode_erasures(code, n, apply_errors_erasures(c, n, flips, erased)) != c:
                            return False
    return True


def erasure_bound_tight(code: list, n: int) -> bool:
    d = min_distance(code)
    for t in range(d // 2 + 1):
        e_count = d - 2 * t
        found = any(decode_erasures(code, n, apply_errors_erasures(c, n, flips, erased)) != c
                    for c in code for erased in combinations(range(n), e_count)
                    for flips in combinations([i for i in range(n) if i not in erased], t))
        if not found:
            return False
    return True


def check_parity() -> None:
    section("Single Parity")
    p = Fraction(1, 100)
    require("closed form equals enumeration for n = 1..10, p = 1/100",
            all(parity_alarm(n, p) == parity_alarm_by_enumeration(n, p) for n in range(1, 11)))
    require("closed form equals enumeration for n = 8, p = 1/3",
            parity_alarm(8, Fraction(1, 3)) == parity_alarm_by_enumeration(8, Fraction(1, 3)))
    check("n = 8, p = 0.01: any error", any_error(8, p), "0.077255")
    check("n = 8, p = 0.01: 0.98^8", Fraction(98, 100) ** 8, "0.850763")
    check("n = 8, p = 0.01: alarm", parity_alarm(8, p), "0.074618")
    check("n = 8, p = 0.01: undetected", any_error(8, p) - parity_alarm(8, p), "0.0026368")
    check("n = 8, p = 0.01: coverage", parity_coverage(8, p), "0.96587")
    check("n = 8, p = 0.01: coverage in percent", parity_coverage(8, p), "96.6%")
    check("small-p estimate 1 - (n-1)p/2, n = 8", 1 - 7 * 0.01 / 2, "0.965")
    check("n = 8, p = 1/2: coverage", parity_coverage(8, Fraction(1, 2)), "0.50196")
    same("n = 8, p = 1/2: coverage as a fraction", parity_coverage(8, Fraction(1, 2)), Fraction(128, 255))
    check("n = 32, p = 0.01: coverage", parity_coverage(32, p), "0.86561")
    check("n = 32, p = 0.01: alarm", parity_alarm(32, p), "0.23806")
    check("n = 32, p = 0.01: any error", any_error(32, p), "0.27502")
    check("n = 32, p = 0.01: corrupted words that slip through", 1 - parity_coverage(32, p), "13%")
    check("n = 8, p = 0.01: corrupted words that slip through", 1 - parity_coverage(8, p), "3.4%")
    check("n = 64, p = 0.01: coverage", parity_coverage(64, p), "0.76469")
    check("n = 64, p = 0.001: coverage", parity_coverage(64, Fraction(1, 1000)), "0.96944")
    require("injected flips: odd counts always caught, even never, n = 8",
            all(weight(0b10110110 ^ e) % 2 != weight(0b10110110) % 2 for e in range(256) if weight(e) % 2)
            and all(weight(0b10110110 ^ e) % 2 == weight(0b10110110) % 2 for e in range(256) if weight(e) % 2 == 0))


def check_isbn() -> None:
    section("Check Digits: ISBN-10")
    payload = digits("0-13-042672")
    products = [(10 - i) * d for i, d in enumerate(payload)]
    same("ISBN-10 products", products, [0, 9, 24, 0, 24, 10, 24, 21, 4])
    same("ISBN-10 payload sum", sum(products), 116)
    same("116 mod 11", 116 % 11, 6)
    same("ISBN-10 check digit", isbn10_check(payload), 5)
    require("0-13-042672-5 is valid", isbn10_valid(digits("0-13-042672-5")))
    swapped = digits("0-13-042627-5")
    same("swapped ISBN-10 weighted sum", isbn10_sum(swapped), 116)
    same("swapped ISBN-10 sum mod 11", isbn10_sum(swapped) % 11, 6)
    require("the swapped ISBN-10 is caught", not isbn10_valid(swapped))
    same("payload 0-12-345672: sum", isbn10_sum(digits("0-12-345672")), 144)
    same("payload 0-12-345672: check", isbn10_check(digits("0-12-345672")), 10)
    require("0-12-345672-X is valid", isbn10_valid(digits("0-12-345672-X")))
    rates = scheme_rates(isbn10_valid, digits("0-13-042672-5"))
    same("ISBN-10 substitutions caught", rates["substitution"], Fraction(1))
    same("ISBN-10 adjacent swaps caught", rates["adjacent"], Fraction(1))
    same("ISBN-10 jump swaps caught", rates["jump"], Fraction(1))
    require("ISBN-10 catches swaps at every distance",
            all(not swap_escapes(isbn10_valid, digits("0-13-042672-5"), i, j, a, b)
                for i, j in combinations(range(10), 2) for a in range(10) for b in range(10) if a != b))

    section("Check Digits: ISBN-13")
    payload = digits("978-0-13-042672")
    terms = [w * d for w, d in zip(isbn13_weights(12), payload)]
    same("ISBN-13 terms", terms, [9, 21, 8, 0, 1, 9, 0, 12, 2, 18, 7, 6])
    same("ISBN-13 payload sum", sum(terms), 93)
    same("ISBN-13 check digit", isbn13_check(payload), 7)
    require("978-0-13-042672-7 is valid", isbn13_valid(digits("978-0-13-042672-7")))
    swapped = digits("978-0-13-042627-7")
    same("swapped ISBN-13 sum", isbn13_sum(swapped), 110)
    require("the swapped ISBN-13 slips through", isbn13_valid(swapped))
    same("3 x 7 mod 10", 3 * 7 % 10, 1)
    rates = scheme_rates(isbn13_valid, digits("978-0-13-042672-7"))
    same("ISBN-13 substitutions caught", rates["substitution"], Fraction(1))
    same("ISBN-13 adjacent swaps caught", rates["adjacent"], Fraction(80, 90))
    check("ISBN-13 adjacent swaps caught, percent", rates["adjacent"], "88.9%")
    same("ISBN-13 jump swaps caught", rates["jump"], Fraction(0))
    escaping = escaping_adjacent_pairs(isbn13_valid, digits("978-0-13-042672-7"))
    same("ISBN-13 adjacent swaps that escape", sorted(escaping),
         sorted((a, b) for a in range(10) for b in range(10) if abs(a - b) == 5))
    same("count of escaping ordered pairs", len(escaping), 10)


def check_luhn() -> None:
    section("Check Digits: the Luhn formula")
    same("doubling map D(0..9)", [luhn_double(x) for x in range(10)], [0, 2, 4, 6, 8, 1, 3, 5, 7, 9])
    require("D is a permutation of 0..9", sorted(luhn_double(x) for x in range(10)) == list(range(10)))
    same("f(x) = D(x) - x", [luhn_double(x) - x for x in range(10)], [0, 1, 2, 3, 4, -4, -3, -2, -1, 0])
    payload = digits("314159265")
    same("Luhn check digit of 314159265", luhn_check(payload), 3)
    number = digits("3141592653")
    same("Luhn terms of 3141592653, left to right", luhn_terms(number), [6, 1, 8, 1, 1, 9, 4, 6, 1, 3])
    same("Luhn sum of 3141592653", luhn_sum(number), 40)
    sub = digits("3141582653")
    same("Luhn sum after 9 -> 8", luhn_sum(sub), 39)
    require("the substitution is caught", not luhn_valid(sub))
    swap = digits("3141952653")
    same("Luhn terms after swapping 5 and 9", luhn_terms(swap), [6, 1, 8, 1, 9, 5, 4, 6, 1, 3])
    same("Luhn sum after swapping 5 and 9", luhn_sum(swap), 44)
    require("the swap is caught", not luhn_valid(swap))
    for base in (digits("3141592653"), digits("4000123412341234")):
        base = make_valid(luhn_valid, base, {})
        rates = scheme_rates(luhn_valid, base)
        same(f"Luhn substitutions caught, length {len(base)}", rates["substitution"], Fraction(1))
        same(f"Luhn adjacent swaps caught, length {len(base)}", rates["adjacent"], Fraction(88, 90))
        same(f"Luhn jump swaps caught, length {len(base)}", rates["jump"], Fraction(0))
        same(f"Luhn adjacent swaps that escape, length {len(base)}",
             sorted(escaping_adjacent_pairs(luhn_valid, base)), [(0, 9), (9, 0)])
    check("Luhn adjacent swaps caught, percent", Fraction(88, 90), "97.8%")
    exercise = digits("50961")
    same("Luhn check digit of 50961", luhn_check(exercise), 2)
    same("Luhn sum of 509612", luhn_sum(digits("509612")), 20)
    same("Luhn sum of 590612", luhn_sum(digits("590612")), 20)
    require("509612 and 590612 both pass", luhn_valid(digits("509612")) and luhn_valid(digits("590612")))

    section("Check digits confirmed on every valid number of five digits")
    luhn5 = [list(map(int, f"{v:05d}")) for v in range(100000) if luhn_valid(list(map(int, f"{v:05d}")))]
    same("valid 5-digit Luhn numbers", len(luhn5), 10000)
    missed_subs = sum(luhn_valid(ds[:i] + [b] + ds[i + 1:]) for ds in luhn5 for i in range(5) for b in range(10) if b != ds[i])
    same("single substitutions that pass", missed_subs, 0)
    escaped = {(ds[i], ds[i + 1]) for ds in luhn5 for i in range(4)
               if ds[i] != ds[i + 1] and luhn_valid(ds[:i] + [ds[i + 1], ds[i]] + ds[i + 2:])}
    same("adjacent swaps that pass", sorted(escaped), [(0, 9), (9, 0)])
    w13 = [1, 3, 1, 3, 1]
    toy13 = [list(map(int, f"{v:05d}")) for v in range(100000)
             if sum(w * d for w, d in zip(w13, map(int, f"{v:05d}"))) % 10 == 0]
    escaped13 = {(ds[i], ds[i + 1]) for ds in toy13 for i in range(4)
                 if ds[i] != ds[i + 1]
                 and sum(w * d for w, d in zip(w13, ds[:i] + [ds[i + 1], ds[i]] + ds[i + 2:])) % 10 == 0}
    same("weights 1,3,1,3,1: adjacent swaps that pass", len(escaped13), 10)


def check_rectangular() -> None:
    section("Rectangular Parity")
    data = [[1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1], [0, 0, 1, 0]]
    array = rect_encode(data)
    same("row parities", [r[4] for r in array[:4]], [1, 0, 1, 1])
    same("column parities", array[4][:4], [0, 0, 1, 0])
    same("corner bit", array[4][4], 1)
    same("parity of the row parities", sum(r[4] for r in array[:4]) % 2, 1)
    same("parity of the column parities", sum(array[4][:4]) % 2, 1)
    hit = [list(r) for r in array]
    hit[2][1] ^= 1
    same("failures after an error at row 3, column 2", rect_failures(hit), ([2], [1]))
    fixed, verdict = rect_decode(hit)
    require("the decoder flips row 3, column 2 back", verdict == "corrected" and fixed == array)
    same("bits in the array", 25, 25)
    check("rate 16/25", Fraction(16, 25), "0.64")
    weights = {}
    for value in range(1 << 16):
        rows = [[value >> (15 - 4 * i - j) & 1 for j in range(4)] for i in range(4)]
        k = weight(flatten(rect_encode(rows)))
        weights[k] = weights.get(k, 0) + 1
    same("smallest nonzero weight over all 2^16 arrays", min(k for k in weights if k), 4)
    same("arrays of weight 4", weights[4], 100)
    same("C(5,2)^2 rectangles", comb(5, 2) ** 2, 100)
    zero = rect_encode([[0] * 4 for _ in range(4)])
    def with_errors(cells):
        a = [list(r) for r in zero]
        for i, j in cells:
            a[i][j] ^= 1
        return a
    cells = [(i, j) for i in range(5) for j in range(5)]
    singles = [rect_decode(with_errors([c]))[1] for c in cells]
    doubles = [rect_decode(with_errors(list(pair)))[1] for pair in combinations(cells, 2)]
    triples = [rect_failures(with_errors(list(tri))) for tri in combinations(cells, 3)]
    same("single errors corrected", singles.count("corrected"), 25)
    same("double errors flagged", doubles.count("flagged"), 300)
    same("triple errors that break some parity", sum(bool(r or c) for r, c in triples), comb(25, 3))
    same("two errors in row 1", rect_failures(with_errors([(0, 0), (0, 2)])), ([], [0, 2]))
    same("errors at (1,1) and (3,4)", rect_failures(with_errors([(0, 0), (2, 3)])), ([0, 2], [0, 3]))
    mimic = rect_failures(with_errors([(0, 0), (0, 2), (2, 0)]))
    same("three errors that mimic one", mimic, ([2], [2]))
    check("R4 rate", Fraction(1, 4), "0.25")
    check("8 x 8 array rate 64/81", Fraction(64, 81), "0.790")


def check_gray() -> None:
    section("Gray Codes")
    g3 = reflected(3)
    same("3-bit reflected Gray code", g3, ["000", "001", "011", "010", "110", "111", "101", "100"])
    g4 = reflected(4)
    same("4-bit reflected Gray code", g4,
         ["0000", "0001", "0011", "0010", "0110", "0111", "0101", "0100",
          "1100", "1101", "1111", "1110", "1010", "1011", "1001", "1000"])
    require("reflected code equals b XOR (b >> 1), n = 1..8",
            all(reflected(n) == [show(gray(b), n) for b in range(1 << n)] for n in range(1, 9)))
    require("neighbors differ in one bit, wrap included, n = 1..8",
            all(distance(gray(b), gray((b + 1) % (1 << n))) == 1 for n in range(1, 9) for b in range(1 << n)))
    require("the bit that flips from b to b+1 is the count of trailing ones of b",
            all(gray(b) ^ gray(b + 1) == 1 << trailing_ones(b) for b in range(255)))
    require("gray_inverse undoes gray for all 12-bit values", all(gray_inverse(gray(b)) == b for b in range(4096)))
    same("1011 XOR 0101", show(word("1011") ^ word("0101"), 4), "1110")
    same("Gray code of 1011", show(gray(word("1011")), 4), "1110")
    same("position 11 of the 4-bit list", g4[11], "1110")
    same("binary of Gray 1110", show(gray_inverse(word("1110")), 4), "1011")
    same("Gray code of 1101", show(gray(word("1101")), 4), "1011")
    same("binary of Gray 0111", show(gray_inverse(word("0111")), 4), "0101")
    same("position 5 of the 4-bit list", g4[5], "0111")
    same("bits that change from 011 to 100 in binary", distance(word("011"), word("100")), 3)
    same("bits that change from 3 to 4 in Gray", distance(gray(3), gray(4)), 1)
    same("Gray codes of 3 and 4", (show(gray(3), 3), show(gray(4), 3)), ("010", "110"))
    readings = {int("".join(bits), 2) for bits in product(*[("0", "1")] * 3)}
    same("binary readings possible mid-transition 3 -> 4", len(readings), 8)
    same("Gray readings possible mid-transition 3 -> 4", sorted(gray_inverse(v) for v in (word("010"), word("110"))), [3, 4])
    same("binary positions of 111 and 000", (word("111"), word("000")), (7, 0))
    same("Gray reading 111 decodes to", gray_inverse(word("111")), 5)
    same("Gray mid-transition readings, 6 -> 7", sorted({gray_inverse(gray(6) & gray(7)), gray_inverse(gray(6) | gray(7))}), [6, 7])


def trailing_ones(b: int) -> int:
    count = 0
    while b & 1:
        count += 1
        b >>= 1
    return count


def check_examples() -> None:
    section("Worked example: from distance three to four")
    extended = [(c << 1) | (weight(c) % 2) for c in C5]
    same("C5 with an overall parity bit", sorted(show(c, 6) for c in extended), sorted(show(c, 6) for c in C6))
    same("C6 distances", sorted(distance(a, b) for a, b in combinations(C6, 2)), [4, 4, 4, 4, 4, 4])
    same("minimum distance of C6", min_distance(C6), 4)
    require("appending parity to any code of odd d gives d + 1, all 4-word codes from 5-bit words with d = 3",
            parity_extension_holds())
    same("C6 corrects 000001 to", show(decode_bounded(C6, word("000001"), 1), 6), "000000")
    same("C6 with two errors 000011", decode_bounded(C6, word("000011"), 1), None)

    section("Worked example: R5 under errors and erasures")
    same("1?01? decodes to", show(decode_erasures(R5, 5, "1?01?"), 5), "11111")
    same("0?01? decodes to", show(decode_erasures(R5, 5, "0?01?"), 5), "00000")
    same("????1 decodes to", show(decode_erasures(R5, 5, "????1"), 5), "11111")
    same("2t + e for one error, two erasures", 2 * 1 + 2, 4)
    same("2t + e for two errors, two erasures", 2 * 2 + 2, 6)
    same("10?0? decodes to", show(decode_erasures(R5, 5, "10?0?"), 5), "00000")

    section("Worked example: rectangular parity with two errors")
    data = [[1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1], [0, 0, 1, 0]]
    array = rect_encode(data)
    hit = [list(r) for r in array]
    hit[0][1] ^= 1
    hit[0][3] ^= 1
    same("errors at (1,2) and (1,4): failures", rect_failures(hit), ([], [1, 3]))
    same("verdict", rect_decode(hit)[1], "flagged")
    hit = [list(r) for r in array]
    hit[0][1] ^= 1
    hit[2][3] ^= 1
    same("errors at (1,2) and (3,4): failures", rect_failures(hit), ([0, 2], [1, 3]))
    same("verdict", rect_decode(hit)[1], "flagged")
    other = [list(r) for r in array]
    other[0][3] ^= 1
    other[2][1] ^= 1
    require("errors at (1,4) and (3,2) give the same failures", rect_failures(other) == rect_failures(hit))


def parity_extension_holds() -> bool:
    for code in combinations(range(32), 4):
        code = list(code)
        d = min_distance(code)
        if d % 2 == 0:
            continue
        extended = [(c << 1) | (weight(c) % 2) for c in code]
        if min_distance(extended) != d + 1:
            return False
    return True


def check_exercises() -> None:
    section("Exercises")
    code = code_of("00000", "00111", "11001", "11110")
    same("exercise code distances", sorted(distance(a, b) for a, b in combinations(code, 2)), [3, 3, 3, 3, 4, 4])
    same("exercise code minimum distance", min_distance(code), 3)
    same("V(7, 2)", ball_size(7, 2), 29)
    same("floor(128 / 29)", 128 // 29, 4)
    check("128 / 29", Fraction(128, 29), "4.41")
    sums = {distance(a, b) + distance(b, c) + distance(a, c) for a, b, c in combinations(range(128), 3)}
    require("three 7-bit words: the three distances sum to an even number at most 14",
            all(s % 2 == 0 for s in sums) and max(sums) == 14)
    require("no three 7-bit words lie pairwise at distance 5 or more",
            not any(min(distance(a, b), distance(b, c), distance(a, c)) >= 5
                    for a, b, c in combinations(range(128), 3)))
    same("(t, e) allowed by d = 7",
         [(t, max(e for e in range(7) if 2 * t + e < 7)) for t in range(4)], [(0, 6), (1, 4), (2, 2), (3, 0)])
    check("8 x 8 rectangular rate", Fraction(64, 81), "0.790")
    same("V(12, 1)", ball_size(12, 1), 13)
    same("floor(4096 / 13)", 4096 // 13, 315)
    require("three words of any length have d(x,y) + d(y,z) + d(x,z) <= 2n, n = 7",
            max(distance(a, b) + distance(b, c) + distance(a, c) for a, b, c in combinations(range(128), 3)) == 14)


def main() -> None:
    check_distance_and_weight()
    check_minimum_distance()
    check_cube()
    check_detection()
    check_correction()
    check_erasures()
    check_parity()
    check_isbn()
    check_luhn()
    check_rectangular()
    check_gray()
    check_examples()
    check_exercises()
    print("\nAll chapter 8 numbers reproduced.")


if __name__ == "__main__":
    main()
