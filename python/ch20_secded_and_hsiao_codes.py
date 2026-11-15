"""Reproduce every number in chapter 20 of Error Correction Coding, "SECDED and Hsiao Codes".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch20_secded_and_hsiao_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in the chapter is a whole number, a bit pattern or a percentage
printed to one decimal place, so every comparison is exact. The distance of
every code is checked by testing every set of one, two and three columns of its
check matrix, and the triple-error counts come from all 59,640 three-bit errors
of the 72-bit codes, not from a sample.

A column of a check matrix with r rows is held as an integer whose bit i is the
entry in row i, and a bit pattern is printed with row 0 first. The module also
holds the code constructions that ch20_charts.py draws, so the figures and this
check share one source.
"""

import sys
from itertools import combinations
from math import comb


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


# ---------------------------------------------------------------- helpers

def weight(v: int) -> int:
    return bin(v).count("1")


def col(rows) -> int:
    """The column with ones in the given rows."""
    v = 0
    for i in rows:
        v |= 1 << i
    return v


def bits(v: int, r: int) -> str:
    """A column or syndrome as a string, row 0 first."""
    return "".join("1" if (v >> i) & 1 else "0" for i in range(r))


def joined(values) -> str:
    return ",".join(str(x) for x in values)


def row_loads(cols, r: int) -> list[int]:
    return [sum((c >> i) & 1 for c in cols) for i in range(r)]


def levels(w: int, fan_in: int = 2) -> int:
    """Levels of fan_in-input gates in a balanced tree over w signals."""
    d, reach = 0, 1
    while reach < w:
        reach *= fan_in
        d += 1
    return d


def ceil_div(a: int, b: int) -> int:
    return -(-a // b)


def pct(part: int, whole: int) -> str:
    return f"{100 * part / whole:.1f}"


# ---------------------------------------------------------------- codes
# Each code is a list of columns, data columns first and check columns last.

def hsiao_8_4() -> list[int]:
    data = [col(c) for c in combinations(range(4), 3)]
    return data + [1 << i for i in range(4)]


def hsiao_13_8() -> list[int]:
    dropped = {col((0, 1, 2)), col((0, 3, 4))}
    data = [col(c) for c in combinations(range(5), 3) if col(c) not in dropped]
    return data + [1 << i for i in range(5)]


def hsiao_22_16() -> list[int]:
    dropped = {col((0, 1, 2)), col((0, 3, 4)), col((1, 3, 5)), col((2, 4, 5))}
    data = [col(c) for c in combinations(range(6), 3) if col(c) not in dropped]
    return data + [1 << i for i in range(6)]


def hsiao_39_32() -> list[int]:
    dropped = {col((0, 1, 2)), col((0, 3, 4)), col((1, 5, 6))}
    data = [col(c) for c in combinations(range(7), 3) if col(c) not in dropped]
    return data + [1 << i for i in range(7)]


def hsiao_72_64() -> list[int]:
    """All 56 weight-3 columns, then the 8 cyclic runs of five rows."""
    data = [col(c) for c in combinations(range(8), 3)]
    data += [col((i + j) % 8 for j in range(5)) for i in range(8)]
    return data + [1 << i for i in range(8)]


def hsiao_72_64_dictionary() -> list[int]:
    """All 56 weight-3 columns, then the first 8 weight-5 columns in dictionary order."""
    data = [col(c) for c in combinations(range(8), 3)]
    data += [col(c) for c in list(combinations(range(8), 5))[:8]]
    return data + [1 << i for i in range(8)]


def hsiao_137_128() -> list[int]:
    """All 84 weight-3 columns of length 9, then 44 weight-5 columns chosen greedily.

    Each step takes the unused weight-5 word whose rows hold the fewest ones so
    far, the first in dictionary order on a tie, which keeps the rows even.
    """
    data = [col(c) for c in combinations(range(9), 3)]
    loads = [1 + comb(8, 2)] * 9
    used = []
    for _ in range(128 - len(data)):
        best = min((c for c in combinations(range(9), 5) if c not in used),
                   key=lambda c: sum(loads[i] for i in c))
        used.append(best)
        for i in best:
            loads[i] += 1
    return data + [col(c) for c in used] + [1 << i for i in range(9)]


def extended_hamming(n: int, rh: int) -> list[int]:
    """Positions 0 to n-1 of the extended Hamming code with rh Hamming rows.

    Rows 0 to rh-1 hold the binary form of the position number, and row rh
    is the overall parity row of all ones. Position 0 is the parity bit.
    """
    return [j | (1 << rh) for j in range(n)]


# ---------------------------------------------------------------- properties

def distance_at_least_4(cols) -> bool:
    if any(c == 0 for c in cols):
        return False
    if len(set(cols)) != len(cols):
        return False
    present = set(cols)
    for a, b in combinations(cols, 2):
        x = a ^ b
        if x in present and x != a and x != b:
            return False
    return True


def has_weight_4_codeword(cols) -> bool:
    index = {}
    for a, b in combinations(range(len(cols)), 2):
        x = cols[a] ^ cols[b]
        for (c, d) in index.get(x, []):
            if len({a, b, c, d}) == 4:
                return True
        index.setdefault(x, []).append((a, b))
    return False


def syndrome(cols, positions) -> int:
    s = 0
    for p in positions:
        s ^= cols[p]
    return s


def triple_outcomes(cols) -> tuple[int, int]:
    present = set(cols)
    mis = det = 0
    for a, b, c in combinations(cols, 3):
        if a ^ b ^ c in present:
            mis += 1
        else:
            det += 1
    return mis, det


def triple_by_weight(cols, r: int) -> list[tuple[int, int]]:
    """For each syndrome weight 0..r, (miscorrected, detected) triple counts."""
    present = set(cols)
    table = [[0, 0] for _ in range(r + 1)]
    for a, b, c in combinations(cols, 3):
        x = a ^ b ^ c
        table[weight(x)][0 if x in present else 1] += 1
    return [tuple(t) for t in table]


def decoder_gates(cols, r: int, k: int, parity_row: int | None) -> dict:
    """Two-input gate counts of the full decoder.

    parity_row is None for a Hsiao code, whose single-error signal is the
    parity of the whole syndrome, and the index of the all-ones row for an
    extended Hamming code, whose single-error signal is that row itself.
    Besides correcting and raising the double-error flag, the decoder matches
    the syndrome against the n - k check columns too, ORs all n match signals
    into M, and raises the no-match flag for an odd syndrome with M = 0.
    """
    loads = row_loads(cols, r)
    n = len(cols)
    out = {}
    out["syndrome"] = sum(w - 1 for w in loads)
    out["syndrome_depth"] = max(levels(w) for w in loads)
    if parity_row is None:
        zero_inputs = r
        out["parity"] = r - 1
        parity_depth = out["syndrome_depth"] + levels(r)
        zero_depth = out["syndrome_depth"] + levels(r)
    else:
        zero_inputs = r - 1
        out["parity"] = 0
        parity_depth = levels(loads[parity_row])
        zero_depth = max(levels(loads[i]) for i in range(r) if i != parity_row) + levels(r - 1)
    out["zero"] = zero_inputs - 1
    out["match"] = k * (r - 1)
    out["correct"] = k
    out["ded"] = 1
    out["check_match"] = (n - k) * (r - 1)
    out["any_match"] = n - 1
    out["no_match"] = 1
    out["total"] = (out["syndrome"] + out["zero"] + out["parity"]
                    + out["match"] + out["correct"] + out["ded"]
                    + out["check_match"] + out["any_match"] + out["no_match"])
    out["inverters"] = r + 2
    match_depth = out["syndrome_depth"] + levels(r)
    out["path"] = match_depth + 1
    out["ded_depth"] = max(parity_depth, zero_depth) + 1
    out["m_depth"] = match_depth + levels(n)
    out["no_match_depth"] = max(out["m_depth"], parity_depth) + 1
    return out


def decoder_actions(cols, r: int, parity_row: int | None, weight_of_error: int) -> dict:
    """Run the full decoder's logic on every error of the given weight.

    Returns how many errors end in each outcome: no flag and no flip,
    a flip of the one bit in error, a flip of some other bit, the
    double-error flag, and the no-match flag.
    """
    present = {c: j for j, c in enumerate(cols)}
    tally = {"silent": 0, "corrected": 0, "miscorrected": 0, "double": 0,
             "no_match": 0, "conflict": 0}
    for err in combinations(range(len(cols)), weight_of_error):
        s = syndrome(cols, err)
        if parity_row is None:
            nonzero, odd = s != 0, weight(s) % 2 == 1
        else:
            nonzero = (s & ~(1 << parity_row)) != 0
            odd = (s >> parity_row) & 1 == 1
        ded = nonzero and not odd
        match = present.get(s)
        no_match = odd and match is None
        if ded + no_match + (match is not None) > 1:
            tally["conflict"] += 1
        elif ded:
            tally["double"] += 1
        elif no_match:
            tally["no_match"] += 1
        elif match is None:
            tally["silent"] += 1
        elif err == (match,):
            tally["corrected"] += 1
        else:
            tally["miscorrected"] += 1
    return tally


# ---------------------------------------------------------------- checks

def smallest_r(test) -> int:
    r = 1
    while not test(r):
        r += 1
    return r


def check_bit_budget() -> None:
    section("How many check bits a 64-bit word needs")
    k = 64
    same("SEC test 2^r >= k + r + 1 at r = 6", f"{2**6} >= {k + 6 + 1}", "64 >= 71")
    same("SEC test at r = 7", f"{2**7} >= {k + 7 + 1}", "128 >= 72")
    same("least r for single correction, k = 64", smallest_r(lambda r: 2**r >= k + r + 1), 7)
    same("SECDED test 2^(r-1) >= k + r at r = 7", f"{2**6} >= {k + 7}", "64 >= 71")
    same("SECDED test at r = 8", f"{2**7} >= {k + 8}", "128 >= 72")
    same("least r for SECDED, k = 64", smallest_r(lambda r: 2**(r - 1) >= k + r), 8)
    same("codeword length 64 + 8", 64 + 8, 72)
    same("bytes in a 72-bit word", 72 // 8, 9)
    same("check-bit overhead 8/64, percent", pct(8, 64), "12.5")
    for kk, rr, nn in [(4, 4, 8), (8, 5, 13), (32, 7, 39), (128, 9, 137)]:
        r = smallest_r(lambda r: 2**(r - 1) >= kk + r)
        same(f"least r for SECDED, k = {kk}", r, rr)
        same(f"length for k = {kk}", kk + r, nn)
    same("least r for single correction, k = 32", smallest_r(lambda r: 2**r >= 32 + r + 1), 6)
    same("SECDED test for k = 32 at r = 6", f"{2**5} >= {32 + 6}", "32 >= 38")
    same("check-bit overhead 7/32, percent", pct(7, 32), "21.9")
    same("longest SECDED length with r = 8", 2**7, 128)
    same("most data bits with r = 8", 2**7 - 8, 120)
    same("longest SECDED length with r = 7", 2**6, 64)
    same("most data bits with r = 7", 2**6 - 7, 57)
    same("SECDED test for k = 128 at r = 8", f"{2**7} >= {128 + 8}", "128 >= 136")
    same("check-bit overhead 9/128, percent", pct(9, 128), "7.0")
    same("least r for single correction, k = 16", smallest_r(lambda r: 2**r >= 16 + r + 1), 5)
    same("least r for SECDED, k = 16", smallest_r(lambda r: 2**(r - 1) >= 16 + r), 6)
    same("longest SECDED length with r = 5", 2**4, 16)
    same("most data bits with r = 5", 2**4 - 5, 11)
    for kk, r1, r2, nn, over in [(4, 3, 4, 8, "100.0"), (8, 4, 5, 13, "62.5"), (16, 5, 6, 22, "37.5"),
                                 (32, 6, 7, 39, "21.9"), (64, 7, 8, 72, "12.5"), (128, 8, 9, 137, "7.0")]:
        a = smallest_r(lambda r: 2**r >= kk + r + 1)
        b = smallest_r(lambda r: 2**(r - 1) >= kk + r)
        same(f"table row k = {kk}: SEC r, SECDED r, n, overhead", f"{a},{b},{kk + b},{pct(b, kk)}",
             f"{r1},{r2},{nn},{over}")
    require("SECDED needs exactly one check bit more than SEC, k = 1 to 4096",
            all(smallest_r(lambda r: 2**(r - 1) >= kk + r) == smallest_r(lambda r: 2**r >= kk + r + 1) + 1
                for kk in range(1, 4097)))


def check_extended_hamming() -> None:
    section("SECDED from the extended Hamming code")
    full = extended_hamming(128, 7)
    require("(128,120) extended Hamming check matrix has distance 4",
            distance_at_least_4(full) and has_weight_4_codeword(full))
    h = extended_hamming(72, 7)
    same("positions kept by the shortening", len(h), 72)
    same("positions deleted from the (128,120) code", 128 - 72, 56)
    checks = [j for j in range(1, 72) if weight(j) == 1]
    same("Hamming check positions among 1..71", joined(checks), "1,2,4,8,16,32,64")
    same("data positions among 1..71", 71 - len(checks), 64)
    require("shortened (72,64) check matrix has distance 4",
            distance_at_least_4(h) and has_weight_4_codeword(h))
    loads = row_loads(h, 8)
    same("ones in each row, rows 0 to 7", joined(loads), "36,36,36,32,32,32,8,72")
    same("ones in the extended Hamming (72,64) check matrix", sum(loads), 284)
    # the decision table on (s', p)
    s = syndrome(h, [37])
    same("syndrome of an error at position 37, rows 0 to 6", bits(s & 127, 7), "1010010")
    same("overall parity of that error", s >> 7, 1)
    s = syndrome(h, [5, 9])
    same("syndrome of errors at positions 5 and 9, rows 0 to 6", bits(s & 127, 7), "0011000")
    same("overall parity of that pair", s >> 7, 0)
    s = syndrome(h, [0])
    same("syndrome of an error in the parity bit, rows 0 to 6", bits(s & 127, 7), "0000000")
    # every single and double error falls in its row of the table
    ok = True
    for i in range(72):
        s = h[i]
        ok &= (s >> 7) == 1
    for i, j in combinations(range(72), 2):
        s = h[i] ^ h[j]
        ok &= (s >> 7) == 0 and (s & 127) != 0
    require("every single error has p = 1, every double p = 0 and s' nonzero", ok)
    h39 = extended_hamming(39, 6)
    require("extended Hamming (39,32) has distance 4",
            distance_at_least_4(h39) and has_weight_4_codeword(h39))
    loads = row_loads(h39, 7)
    same("ones in each row of extended Hamming (39,32)", joined(loads), "19,19,19,16,16,7,39")
    same("ones in the extended Hamming (39,32) check matrix", sum(loads), 135)
    # replace the all-ones row by the sum of all eight rows: the same code, odd columns
    odd_form = [(j & 127) | (((1 + weight(j & 127)) % 2) << 7) for j in h]
    require("odd-column form: every column has odd weight", all(weight(c) % 2 for c in odd_form))
    require("odd-column form still has distance 4",
            distance_at_least_4(odd_form) and has_weight_4_codeword(odd_form))
    same("positions 0..71 whose index has even weight", sum(1 for j in range(72) if weight(j) % 2 == 0), 36)
    same("ones in the odd-column form, 212 + 36", sum(row_loads(odd_form, 8)), 248)


def check_odd_weight() -> None:
    section("Odd-weight columns")
    for r, n in [(4, 8), (5, 16), (6, 32), (7, 64), (8, 128)]:
        same(f"odd-weight words of length {r}", sum(1 for v in range(2**r) if weight(v) % 2), n)
    by_weight = [sum(1 for v in range(256) if weight(v) == w) for w in (1, 3, 5, 7)]
    same("8-bit odd words of weight 1, 3, 5, 7", joined(by_weight), "8,56,56,8")
    ok = True
    for a in range(1, 64):
        for b in range(1, 64):
            w = weight(a ^ b)
            ok &= w == weight(a) + weight(b) - 2 * weight(a & b)
    require("wt(a+b) = wt(a) + wt(b) - 2 wt(a AND b), all 6-bit pairs", ok)


def check_hsiao_8_4() -> None:
    section("The Hsiao (8,4) code")
    h = hsiao_8_4()
    same("columns d0..d3 c0..c3", " ".join(bits(c, 4) for c in h),
         "1110 1101 1011 0111 1000 0100 0010 0001")
    require("every column has odd weight", all(weight(c) % 2 for c in h))
    require("distance is exactly 4", distance_at_least_4(h) and has_weight_4_codeword(h))
    same("ones in each row", joined(row_loads(h, 4)), "4,4,4,4")
    same("ones in the check matrix", sum(row_loads(h, 4)), 16)
    same("columns used out of the 8 odd words", len(set(h)), 8)
    data = [1, 0, 1, 1]
    checks = [sum(data[j] for j in range(4) if (h[j] >> i) & 1) % 2 for i in range(4)]
    same("check bits for data 1011", "".join(map(str, checks)), "0010")
    same("syndrome of an error in d2", bits(syndrome(h, [2]), 4), "1011")
    same("syndrome of errors in d0 and c1", bits(syndrome(h, [0, 5]), 4), "1010")
    counts = {}
    for i, j in combinations(range(8), 2):
        s = h[i] ^ h[j]
        counts[s] = counts.get(s, 0) + 1
    same("double errors", sum(counts.values()), 28)
    same("distinct double-error syndromes", len(counts), 7)
    require("every double-error syndrome is even and nonzero",
            all(s != 0 and weight(s) % 2 == 0 for s in counts))
    same("pairs behind each double-error syndrome", joined(sorted(set(counts.values()))), "4")
    pairs = sorted((i, j) for i, j in combinations(range(8), 2) if h[i] ^ h[j] == col((0, 1)))
    names = ["d0", "d1", "d2", "d3", "c0", "c1", "c2", "c3"]
    classes = {}
    for i, j in combinations(range(8), 2):
        classes.setdefault(bits(h[i] ^ h[j], 4), []).append(names[i] + names[j])
    for syn, stated in [("1111", "d0c3 d1c2 d2c1 d3c0"), ("1100", "d0c2 d1c3 d2d3 c0c1"),
                        ("1010", "d0c1 d1d3 d2c3 c0c2"), ("1001", "d0d3 d1c1 d2c2 c0c3"),
                        ("0110", "d0c0 d1d2 d3c3 c1c2"), ("0101", "d0d2 d1c0 d3c2 c1c3"),
                        ("0011", "d0d1 d2c0 d3c1 c2c3")]:
        same(f"pairs with syndrome {syn}", " ".join(classes[syn]), stated)
    same("syndromes: 1 zero + 8 single + 7 double", 1 + 8 + len(classes), 16)


def check_hsiao_13_8() -> None:
    section("The Hsiao (13,8) code")
    h = hsiao_13_8()
    same("data columns d0..d7", " ".join(bits(c, 5) for c in h[:8]),
         "11010 11001 10110 10101 01110 01101 01011 00111")
    require("every column has odd weight", all(weight(c) % 2 for c in h))
    require("distance is exactly 4", distance_at_least_4(h) and has_weight_4_codeword(h))
    loads = row_loads(h, 5)
    same("ones in each row", joined(loads), "5,6,6,6,6")
    same("ones in the check matrix", sum(loads), 29)
    same("data bits in each check", joined(w - 1 for w in loads), "4,5,5,5,5")
    odd_unused = [v for v in range(32) if weight(v) % 2 and v not in h]
    same("odd words that are not columns", " ".join(bits(v, 5) for v in odd_unused),
         "11100 10011 11111")
    data = [1, 0, 1, 1, 0, 0, 1, 0]
    checks = [sum(data[j] for j in range(8) if (h[j] >> i) & 1) % 2 for i in range(5)]
    same("check bits for data 10110010", "".join(map(str, checks)), "10010")
    same("weight of that codeword", sum(data) + sum(checks), 6)
    same("syndrome of an error in d5", bits(syndrome(h, [5]), 5), "01101")
    same("syndrome of errors in d0 and c2", bits(syndrome(h, [0, 10]), 5), "11110")
    s = syndrome(h, [0, 1, 2])
    same("syndrome of errors in d0, d1, d2", bits(s, 5), "10101")
    same("column it matches", h.index(s), 3)
    same("syndrome of errors in d0, d7, c3", bits(syndrome(h, [0, 7, 11]), 5), "11111")
    same("syndrome of errors in d0, c2, c3", bits(syndrome(h, [0, 10, 11]), 5), "11100")
    same("syndrome of d0, c0, c1, c3", bits(syndrome(h, [0, 8, 9, 11]), 5), "00000")
    doubles = {a ^ b for a, b in combinations(h, 2)}
    same("double errors in 13 bits, C(13,2)", comb(13, 2), 78)
    same("distinct double-error syndromes", len(doubles), 15)
    require("every double-error syndrome is even and nonzero",
            all(s != 0 and weight(s) % 2 == 0 for s in doubles))
    same("syndromes: 1 zero + 13 single + 3 unused odd + 15 even", 1 + 13 + 3 + len(doubles), 32)


def check_hsiao_72_64() -> None:
    section("The Hsiao (72,64) and (39,32) codes")
    h = hsiao_72_64()
    same("columns", len(h), 72)
    require("all 72 columns distinct and odd", len(set(h)) == 72 and all(weight(c) % 2 for c in h))
    require("distance is exactly 4", distance_at_least_4(h) and has_weight_4_codeword(h))
    same("data columns of weight 3 and 5", f"{sum(weight(c) == 3 for c in h)},{sum(weight(c) == 5 for c in h)}", "56,8")
    same("weight-5 columns", " ".join(bits(c, 8) for c in h[56:64]),
         "11111000 01111100 00111110 00011111 10001111 11000111 11100011 11110001")
    loads = row_loads(h, 8)
    same("ones in each row", joined(loads), "27,27,27,27,27,27,27,27")
    same("ones in the check matrix, 8 + 3*56 + 5*8", sum(loads), 216)
    same("8 + 3*56 + 5*8", 8 + 3 * 56 + 5 * 8, 216)
    same("ones in the data part of each row", joined(w - 1 for w in loads), "26,26,26,26,26,26,26,26")
    same("weight-3 columns that contain one given row, C(7,2)", comb(7, 2), 21)
    same("cyclic runs of five rows that contain row 0", sum((h[56 + i] >> 0) & 1 for i in range(8)), 5)
    same("least possible largest row, ceil(216/8)", ceil_div(216, 8), 27)
    same("syndrome of the four-bit error d0 c0 c1 c2", bits(syndrome(h, [0, 64, 65, 66]), 8), "00000000")
    hd = hsiao_72_64_dictionary()
    require("dictionary-order choice has distance 4", distance_at_least_4(hd) and has_weight_4_codeword(hd))
    loads_d = row_loads(hd, 8)
    same("ones in each row, dictionary order", joined(loads_d), "30,30,30,26,26,25,25,24")
    same("ones in the dictionary-order matrix", sum(loads_d), 216)
    same("two-input levels for 27 and 30 inputs", f"{levels(27)},{levels(30)}", "5,5")
    same("three-input levels for 27 and 30 inputs", f"{levels(27, 3)},{levels(30, 3)}", "3,4")
    same("two-input levels for 72 inputs", levels(72), 7)
    h39 = hsiao_39_32()
    require("Hsiao (39,32) has distance exactly 4", distance_at_least_4(h39) and has_weight_4_codeword(h39))
    same("weight-3 columns of length 7, C(7,3)", comb(7, 3), 35)
    loads39 = row_loads(h39, 7)
    same("ones in each row of Hsiao (39,32)", joined(loads39), "14,14,15,15,15,15,15")
    same("ones in Hsiao (39,32)", sum(loads39), 103)
    same("largest row: ext. Hamming, dictionary, balanced (72,64)",
         f"{max(row_loads(extended_hamming(72, 7), 8))},{max(loads_d)},{max(loads)}", "72,30,27")
    same("largest row: ext. Hamming and Hsiao (39,32)",
         f"{max(row_loads(extended_hamming(39, 6), 7))},{max(loads39)}", "39,15")
    same("lower bound 3k + r for k = 32, r = 7", 3 * 32 + 7, 103)
    same("least possible largest row, ceil(103/7)", ceil_div(103, 7), 15)
    same("weight-3 columns of length 7 that contain a row, C(6,2)", comb(6, 2), 15)
    h22 = hsiao_22_16()
    require("Hsiao (22,16) has distance exactly 4", distance_at_least_4(h22) and has_weight_4_codeword(h22))
    same("ones in each row of Hsiao (22,16)", joined(row_loads(h22, 6)), "9,9,9,9,9,9")
    same("ones in Hsiao (22,16), 6 + 3*16", sum(row_loads(h22, 6)), 54)
    same("weight-3 columns of length 6, C(6,3)", comb(6, 3), 20)
    # a 128-bit word
    same("weight-3 columns of length 9, C(9,3)", comb(9, 3), 84)
    same("weight-5 columns needed for k = 128", 128 - comb(9, 3), 44)
    same("ones in Hsiao (137,128), 9 + 3*84 + 5*44", 9 + 3 * 84 + 5 * 44, 481)
    same("least possible largest row, ceil(481/9)", ceil_div(481, 9), 54)
    same("two-input levels for 54 inputs", levels(54), 6)
    h137 = hsiao_137_128()
    require("Hsiao (137,128) columns are 137 distinct odd words",
            len(set(h137)) == 137 and all(weight(c) % 2 for c in h137))
    same("ones in each row of a balanced Hsiao (137,128)", joined(row_loads(h137, 9)),
         "54,54,54,53,54,53,53,53,53")
    # systematic SECDED codes have no data column of weight 2
    ok = True
    for i, j in combinations(range(8), 2):
        ok &= not distance_at_least_4([1 << t for t in range(8)] + [col((i, j))])
    require("a weight-2 column beside I_8 always breaks distance 4", ok)


def check_decoder() -> None:
    section("The decoder in gates")
    same("gates in a tree over w inputs, w = 27", 27 - 1, 26)
    same("levels of a tree over 27 inputs", levels(27), 5)
    hs = decoder_gates(hsiao_72_64(), 8, 64, None)
    eh = decoder_gates(extended_hamming(72, 7), 8, 64, 7)
    same("Hsiao (72,64) syndrome XOR gates", hs["syndrome"], 208)
    same("Hsiao (72,64) syndrome levels", hs["syndrome_depth"], 5)
    same("Hsiao (72,64) zero detector OR gates", hs["zero"], 7)
    same("Hsiao (72,64) parity detector XOR gates", hs["parity"], 7)
    same("Hsiao (72,64) match AND gates, 64 x 7", hs["match"], 448)
    same("Hsiao (72,64) correction XOR gates", hs["correct"], 64)
    same("Hsiao (72,64) double-error AND gate", hs["ded"], 1)
    same("Hsiao (72,64) two-input gates in all", hs["total"], 863)
    same("Hsiao (72,64) inverters", hs["inverters"], 10)
    same("Hsiao (72,64) levels to corrected data", hs["path"], 9)
    same("Hsiao (72,64) levels to the double-error flag", hs["ded_depth"], 9)
    same("Hsiao (72,64) encoder XOR gates, 8 x 25", 8 * (26 - 1), 200)
    same("extended Hamming (72,64) syndrome XOR gates", eh["syndrome"], 276)
    same("extended Hamming (72,64) syndrome levels", eh["syndrome_depth"], 7)
    same("extended Hamming (72,64) zero detector OR gates", eh["zero"], 6)
    same("extended Hamming (72,64) parity detector gates", eh["parity"], 0)
    same("extended Hamming (72,64) two-input gates in all", eh["total"], 923)
    same("extended Hamming (72,64) levels to corrected data", eh["path"], 11)
    same("extended Hamming (72,64) levels to the double-error flag", eh["ded_depth"], 10)
    same("syndrome gates saved by Hsiao (72,64)", eh["syndrome"] - hs["syndrome"], 68)
    same("gates saved by Hsiao (72,64) in all", eh["total"] - hs["total"], 60)
    same("gates Hsiao (72,64) spends that extended Hamming does not", 68 - (eh["total"] - hs["total"]), 8)
    hd = decoder_gates(hsiao_72_64_dictionary(), 8, 64, None)
    same("dictionary-order syndrome XOR gates and levels", f"{hd['syndrome']},{hd['syndrome_depth']}", "208,5")
    h39 = decoder_gates(hsiao_39_32(), 7, 32, None)
    e39 = decoder_gates(extended_hamming(39, 6), 7, 32, 6)
    same("Hsiao (39,32) syndrome XOR gates", h39["syndrome"], 96)
    same("Hsiao (39,32) syndrome levels", h39["syndrome_depth"], 4)
    same("Hsiao (39,32) zero detector OR gates", h39["zero"], 6)
    same("Hsiao (39,32) parity detector XOR gates", h39["parity"], 6)
    same("Hsiao (39,32) match AND gates, 32 x 6", h39["match"], 192)
    same("Hsiao (39,32) two-input gates in all", h39["total"], 414)
    same("Hsiao (39,32) levels to corrected data", h39["path"], 8)
    same("Hsiao (39,32) levels to the double-error flag", h39["ded_depth"], 8)
    same("extended Hamming (39,32) syndrome XOR gates", e39["syndrome"], 128)
    same("extended Hamming (39,32) syndrome levels", e39["syndrome_depth"], 6)
    same("extended Hamming (39,32) two-input gates in all", e39["total"], 439)
    same("extended Hamming (39,32) zero detector OR gates", e39["zero"], 5)
    same("extended Hamming (39,32) levels to corrected data", e39["path"], 10)
    same("extended Hamming (39,32) levels to the double-error flag", e39["ded_depth"], 9)
    same("Hsiao (72,64) check-column match AND gates, 8 x 7", hs["check_match"], 56)
    same("Hsiao (72,64) any-match OR gates, 72 - 1", hs["any_match"], 71)
    same("Hsiao (72,64) no-match AND gate", hs["no_match"], 1)
    same("gates the no-match flag adds, 56 + 71 + 1",
         hs["check_match"] + hs["any_match"] + hs["no_match"], 128)
    same("Hsiao (72,64) gates before the no-match flag", hs["total"] - 128, 735)
    same("Hsiao (72,64) levels to the any-match signal M", hs["m_depth"], 15)
    same("Hsiao (72,64) levels to the no-match flag", hs["no_match_depth"], 16)
    same("extended Hamming (72,64) levels to M and to the no-match flag",
         f"{eh['m_depth']},{eh['no_match_depth']}", "17,18")
    same("Hsiao (39,32) gates the no-match flag adds, 42 + 38 + 1",
         h39["check_match"] + h39["any_match"] + h39["no_match"], 81)
    same("Hsiao (39,32) levels to the no-match flag", h39["no_match_depth"], 14)
    same("extended Hamming (39,32) levels to the no-match flag", e39["no_match_depth"], 16)
    same("Hsiao (39,32) gates before the no-match flag", h39["total"] - 81, 333)
    same("check match, any-match OR and no-match AND of the four decoders",
         " ".join(f"{d['check_match']},{d['any_match']},{d['no_match']}" for d in (hs, eh, h39, e39)),
         "56,71,1 56,71,1 42,38,1 42,38,1")
    for name, cols, r, prow in (("Hsiao (72,64)", hsiao_72_64(), 8, None),
                                ("extended Hamming (72,64)", extended_hamming(72, 7), 8, 7)):
        one = decoder_actions(cols, r, prow, 1)
        two = decoder_actions(cols, r, prow, 2)
        same(f"{name} decoder: single errors corrected", one["corrected"], 72)
        same(f"{name} decoder: double errors flagged, C(72,2)", two["double"], 2556)
        require(f"{name} decoder: no single or double error ends any other way",
                sum(one.values()) == one["corrected"] and sum(two.values()) == two["double"])


def check_triples() -> None:
    section("Three errors")
    ok = True
    h = hsiao_72_64()
    for a, b, c in combinations(h, 3):
        ok &= weight(a ^ b ^ c) % 2 == 1 and (a ^ b ^ c) not in (a, b, c)
    require("a triple's syndrome is odd and is none of its three columns", ok)
    same("three-bit errors in 72 bits, C(72,3)", comb(72, 3), 59640)
    mis, det = triple_outcomes(h)
    same("Hsiao (72,64) triples miscorrected", mis, 33568)
    same("Hsiao (72,64) triples detected", det, 26072)
    same("Hsiao (72,64) percent miscorrected", pct(mis, mis + det), "56.3")
    same("72 columns out of 128 odd words, percent", f"{100 * 72 / 128:.2f}", "56.25")
    tw = triple_by_weight(h, 8)
    same("Hsiao (72,64) triples by syndrome weight 1, 3, 5, 7",
         joined(tw[w][0] + tw[w][1] for w in (1, 3, 5, 7)), "3952,26112,26040,3536")
    same("Hsiao (72,64) weight-1 and weight-3 triples detected", tw[1][1] + tw[3][1], 0)
    same("Hsiao (72,64) weight-5 triples miscorrected", tw[5][0], 3504)
    same("Hsiao (72,64) weight-5 triples detected", tw[5][1], 22536)
    same("Hsiao (72,64) weight-7 triples miscorrected", tw[7][0], 0)
    mis_e, det_e = triple_outcomes(extended_hamming(72, 7))
    same("extended Hamming (72,64) triples miscorrected", mis_e, 45304)
    same("extended Hamming (72,64) triples detected", det_e, 14336)
    same("extended Hamming (72,64) percent miscorrected", pct(mis_e, mis_e + det_e), "76.0")
    eh = extended_hamming(72, 7)
    present = set(eh)
    inside = sum(1 for a, b, c in combinations(eh[:64], 3) if a ^ b ^ c in present)
    same("triples inside positions 0..63, C(64,3)", comb(64, 3), 41664)
    same("triples inside positions 0..63 that miscorrect", inside, 41664)
    mis_d, det_d = triple_outcomes(hsiao_72_64_dictionary())
    same("dictionary-order Hsiao (72,64) triples miscorrected", mis_d, 34032)
    same("dictionary-order Hsiao (72,64) triples detected", det_d, 25608)
    same("dictionary-order Hsiao (72,64) percent miscorrected", pct(mis_d, mis_d + det_d), "57.1")
    same("extended Hamming miscorrected triples touching 64..71", mis_e - 41664, 3640)
    same("more miscorrected triples in dictionary order than balanced", mis_d - mis, 464)
    mis39, det39 = triple_outcomes(hsiao_39_32())
    same("three-bit errors in 39 bits, C(39,3)", comb(39, 3), 9139)
    same("Hsiao (39,32) triples miscorrected", mis39, 5452)
    same("Hsiao (39,32) triples detected", det39, 3687)
    same("Hsiao (39,32) percent miscorrected", pct(mis39, 9139), "59.7")
    mis39e, det39e = triple_outcomes(extended_hamming(39, 6))
    same("extended Hamming (39,32) triples miscorrected", mis39e, 6332)
    same("extended Hamming (39,32) triples detected", det39e, 2807)
    same("extended Hamming (39,32) percent miscorrected", pct(mis39e, 9139), "69.3")
    mis13, det13 = triple_outcomes(hsiao_13_8())
    same("three-bit errors in 13 bits, C(13,3)", comb(13, 3), 286)
    same("Hsiao (13,8) triples miscorrected", mis13, 220)
    same("Hsiao (13,8) triples detected", det13, 66)
    same("Hsiao (13,8) percent miscorrected", pct(mis13, 286), "76.9")
    mis8, det8 = triple_outcomes(hsiao_8_4())
    same("Hsiao (8,4) triples miscorrected", mis8, 56)
    same("Hsiao (8,4) triples detected", det8, 0)
    for name, cols, r, prow, stated in (
            ("Hsiao (72,64)", h, 8, None, "33568,26072"),
            ("dictionary-order Hsiao (72,64)", hsiao_72_64_dictionary(), 8, None, "34032,25608"),
            ("extended Hamming (72,64)", extended_hamming(72, 7), 8, 7, "45304,14336"),
            ("Hsiao (39,32)", hsiao_39_32(), 7, None, "5452,3687"),
            ("extended Hamming (39,32)", extended_hamming(39, 6), 7, 6, "6332,2807"),
            ("Hsiao (13,8)", hsiao_13_8(), 5, None, "220,66"),
            ("Hsiao (8,4)", hsiao_8_4(), 4, None, "56,0")):
        t = decoder_actions(cols, r, prow, 3)
        same(f"{name} full decoder on triples: miscorrected, no-match flag",
             f"{t['miscorrected']},{t['no_match']}", stated)
        require(f"{name} full decoder: no triple passes unflagged or raises the double-error flag",
                t["silent"] + t["double"] + t["conflict"] == 0)


def main() -> None:
    check_bit_budget()
    check_extended_hamming()
    check_odd_weight()
    check_hsiao_8_4()
    check_hsiao_13_8()
    check_hsiao_72_64()
    check_decoder()
    check_triples()
    print("\nAll chapter 20 numbers reproduced.")


if __name__ == "__main__":
    main()
