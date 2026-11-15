"""Reproduce every number in chapter 19 of Error Correction Coding, "Hamming Codes and Perfect Codes".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch19_hamming_codes_and_perfect_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in the chapter is a whole number, a binary or ternary word, a
matrix or a code rate printed to three places, so every comparison is exact.
The small codes are checked exhaustively: every one of the 2,048 messages of
the (15,11) code is encoded and every single error on every codeword is
corrected, the 16 balls of radius one around the (7,4) codewords are shown to
cover all 128 words of length 7 exactly once, every double error on the
extended (8,4) and (16,11) codes is shown to be detected, and the gate-level
decoder is simulated on every correctable received word of both binary codes.

The binary Hamming code of redundancy m is built here with its columns in
counting order: column j of H is the number j written in m bits, most
significant bit in the top row. A word is a list of bits, and position j of a
word is index j - 1 of the list. The module also holds the functions that
ch19_charts.py draws, so the figures and this check share one source.
"""

import sys
from itertools import combinations, product
from math import comb


# Helpers ----------------------------------------------------------------------

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


def word(v) -> str:
    return "".join(str(x) for x in v)


def show(rows) -> str:
    return " ".join(word(row) for row in rows)


def parse(text: str) -> list:
    return [int(ch) for ch in text]


def weight(v) -> int:
    return sum(1 for x in v if x != 0)


def rate(k: int, n: int) -> str:
    return f"{k / n:.3f}"


# Binary Hamming codes in counting order ------------------------------------------

def hamming_H(m: int) -> list:
    """The m-by-(2^m - 1) parity-check matrix whose column j is j in binary."""
    n = 2 ** m - 1
    return [[(j >> (m - 1 - i)) & 1 for j in range(1, n + 1)] for i in range(m)]


def columns_H(n: int, m: int) -> list:
    """H with columns 1, 2, ..., n in m bits, for a code shorter than 2^m - 1."""
    return [[(j >> (m - 1 - i)) & 1 for j in range(1, n + 1)] for i in range(m)]


def syndrome(H: list, v: list) -> list:
    return [sum(h * x for h, x in zip(row, v)) % 2 for row in H]


def index_of(bits: list) -> int:
    value = 0
    for b in bits:
        value = 2 * value + b
    return value


def check_positions(n: int) -> list:
    return [j for j in range(1, n + 1) if j & (j - 1) == 0]


def data_positions(n: int) -> list:
    return [j for j in range(1, n + 1) if j & (j - 1) != 0]


def encode(n: int, data: list) -> list:
    """Place the data at the non-power-of-two positions and fill each check bit."""
    c = [0] * (n + 1)
    for p, b in zip(data_positions(n), data):
        c[p] = b
    for p in check_positions(n):
        c[p] = sum(c[j] for j in range(1, n + 1) if j & p and j != p) % 2
    return c[1:]


def extract(n: int, c: list) -> list:
    return [c[p - 1] for p in data_positions(n)]


def decode(n: int, H: list, v: list) -> list:
    """Flip the bit the syndrome names, if any."""
    j = index_of(syndrome(H, v))
    out = v[:]
    if j:
        out[j - 1] ^= 1
    return out


def flip(v: list, *positions) -> list:
    out = v[:]
    for p in positions:
        out[p - 1] ^= 1
    return out


def all_codewords(n: int) -> list:
    k = len(data_positions(n))
    return [encode(n, list(msg)) for msg in product((0, 1), repeat=k)]


def ball(v: list, radius: int) -> list:
    words = []
    for t in range(radius + 1):
        for pos in combinations(range(1, len(v) + 1), t):
            words.append(flip(v, *pos))
    return words


def volume(n: int, t: int, q: int = 2) -> int:
    return sum(comb(n, i) * (q - 1) ** i for i in range(t + 1))


def min_redundancy(k: int) -> int:
    m = 1
    while k + m + 1 > 2 ** m:
        m += 1
    return m


H7 = hamming_H(3)
H15 = hamming_H(4)
C7 = all_codewords(7)
C15 = all_codewords(15)
FAMILY = [(m, 2 ** m - 1, 2 ** m - 1 - m) for m in range(2, 8)]
M15 = parse("10110011010")
CW15 = encode(15, M15)
RX15 = flip(CW15, 11)


# The ternary Hamming code of length 4 ------------------------------------------------

H3 = [[0, 1, 1, 1], [1, 0, 1, 2]]


def syndrome_q(H: list, v: list, q: int) -> list:
    return [sum(h * x for h, x in zip(row, v)) % q for row in H]


def ternary_codewords() -> list:
    return [v for v in product(range(3), repeat=4) if syndrome_q(H3, list(v), 3) == [0, 0]]


def ternary_decode(v: list) -> tuple:
    """Return (position, value, corrected word) for a word with at most one error."""
    s = syndrome_q(H3, v, 3)
    if s == [0, 0]:
        return 0, 0, v[:]
    for j in range(4):
        col = [H3[0][j], H3[1][j]]
        for a in (1, 2):
            if [(a * x) % 3 for x in col] == s:
                out = v[:]
                out[j] = (out[j] - a) % 3
                return j + 1, a, out
    raise ValueError("no column matches")


def projective_columns(m: int, q: int) -> list:
    """One column from each one-dimensional subspace: first nonzero entry 1."""
    cols = []
    for v in product(range(q), repeat=m):
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            cols.append(list(v))
    return cols


# Extended Hamming codes ----------------------------------------------------------

def extend(c: list) -> list:
    """Put the overall parity bit in front, as position 0."""
    return [sum(c) % 2] + c


def extended_decision(v: list, H: list) -> str:
    s = index_of(syndrome(H, v[1:]))
    p = sum(v) % 2
    if s == 0 and p == 0:
        return "no error"
    if p == 1 and s == 0:
        return "flip the parity bit"
    if p == 1:
        return f"flip position {s}"
    return "double error detected"


# Gate-level decoder --------------------------------------------------------------

def gate_decoder(H: list, v: list) -> tuple:
    """Syndrome XOR trees, an AND gate per position with inverted inputs where
    the column has a 0, and one correction XOR per position."""
    m, n = len(H), len(v)
    s = [sum(H[i][j] * v[j] for j in range(n)) % 2 for i in range(m)]
    lines = []
    for j in range(n):
        a = 1
        for i in range(m):
            a &= s[i] if H[i][j] else 1 - s[i]
        lines.append(a)
    return lines, [x ^ a for x, a in zip(v, lines)]


def encoder_counts(n: int) -> tuple:
    weights = [len([j for j in data_positions(n) if j & p]) for p in check_positions(n)]
    gates = sum(w - 1 for w in weights)
    depth = max((w - 1).bit_length() for w in weights)
    return weights, gates, depth


def decoder_counts(H: list) -> dict:
    m, n = len(H), len(H[0])
    weights = [sum(row) for row in H]
    return {
        "syndrome": sum(w - 1 for w in weights),
        "syndrome depth": max((w - 1).bit_length() for w in weights),
        "and": n,
        "and2": n * (m - 1),
        "not": m,
        "xor": n,
    }


# Checks -------------------------------------------------------------------------

def check_counting() -> None:
    section("How many check bits one error needs")
    same("syndromes of 3 check bits", 2 ** 3, 8)
    same("outcomes to tell apart for n = 7: no error plus 7 positions", 7 + 1, 8)
    table = [(k, min_redundancy(k), k + min_redundancy(k)) for k in (1, 4, 8, 11, 16, 26, 32, 57, 64, 120)]
    same("least m for k = 1 4 8 11 16 26 32 57 64 120",
         " ".join(str(m) for _, m, _ in table), "2 3 4 4 5 5 6 6 7 7")
    same("the matching lengths n = k + m",
         " ".join(str(n) for _, _, n in table), "3 7 12 15 21 31 38 63 71 127")
    require("k = 4 with m = 3 meets n + 1 = 2^m exactly", 4 + 3 + 1 == 2 ** 3)
    require("k = 8 with m = 3 fails, 12 > 8", 8 + 3 + 1 > 2 ** 3)


def check_hamming_7_4() -> None:
    section("The (7,4) code in counting order")
    same("H for m = 3", show(H7), "0001111 0110011 1010101")
    same("check positions", " ".join(map(str, check_positions(7))), "1 2 4")
    same("data positions", " ".join(map(str, data_positions(7))), "3 5 6 7")
    same("codewords", len(C7), 16)
    require("every codeword has zero syndrome", all(syndrome(H7, c) == [0, 0, 0] for c in C7))
    same("minimum weight", min(weight(c) for c in C7 if any(c)), 3)
    c = encode(7, parse("1011"))
    same("message 1011 encodes to", word(c), "0110011")
    rx = flip(c, 6)
    same("error at position 6 gives", word(rx), "0110001")
    same("its syndrome", word(syndrome(H7, rx)), "110")
    same("which is the binary for", index_of(syndrome(H7, rx)), 6)
    require("every single error on every codeword is corrected",
            all(decode(7, H7, flip(c, j)) == c for c in C7 for j in range(1, 8)))
    same("cosets of the (7,4) code, distinct syndromes of the 128 words",
         len({tuple(syndrome(H7, list(v))) for v in product((0, 1), repeat=7)}), 8)


def check_distance() -> None:
    section("Dimension and distance three")
    for m in (2, 3, 4, 5):
        H = hamming_H(m)
        n = 2 ** m - 1
        cols = [tuple(H[i][j] for i in range(m)) for j in range(n)]
        require(f"m = {m}: no zero column and no two columns equal",
                all(any(c) for c in cols) and len(set(cols)) == n)
        require(f"m = {m}: columns 1, 2 and 3 add to zero",
                all((cols[0][i] + cols[1][i] + cols[2][i]) % 2 == 0 for i in range(m)))
        units = [cols[(1 << i) - 1] for i in range(m)]
        require(f"m = {m}: the columns at 1, 2, 4, ... are the m unit vectors",
                sorted(units) == sorted(tuple(int(i == u) for i in range(m)) for u in range(m)))
    for m, code in ((2, all_codewords(3)), (3, C7), (4, C15)):
        n = 2 ** m - 1
        same(f"m = {m}: log2 of the number of codewords, n - m", len(code).bit_length() - 1, n - m)
        same(f"m = {m}: minimum distance by enumeration", min(weight(c) for c in code if any(c)), 3)
    same("weight-3 codewords of the (7,4) code", sum(1 for c in C7 if weight(c) == 3), 7)
    same("C(7,2) / 3", comb(7, 2) // 3, 7)
    w3 = flip([0] * 7, 1, 2, 3)
    same("the weight-3 word on columns 1, 2, 3", word(w3), "1110000")
    require("and it is a codeword, with syndrome zero", not any(syndrome(H7, w3)))


def check_hamming_15_11() -> None:
    section("The (15,11) code by hand")
    same("H for m = 4", show(H15), "000000011111111 000111100001111 011001100110011 101010101010101")
    same("data positions", " ".join(map(str, data_positions(15))), "3 5 6 7 9 10 11 12 13 14 15")
    same("data bits feeding c1 c2 c4 c8", " ".join(str(w) for w in encoder_counts(15)[0]), "7 7 7 7")
    same("message", word(M15), "10110011010")
    same("check bits c1 c2 c4 c8", word([CW15[p - 1] for p in (1, 2, 4, 8)]), "1101")
    same("codeword", word(CW15), "111001110011010")
    same("weight of the codeword", weight(CW15), 9)
    same("received with position 11 flipped", word(RX15), "111001110001010")
    same("syndrome", word(syndrome(H15, RX15)), "1011")
    same("which names position", index_of(syndrome(H15, RX15)), 11)
    same("decoded data", word(extract(15, decode(15, H15, RX15))), "10110011010")
    require("all 2,048 codewords have zero syndrome", all(syndrome(H15, c) == [0] * 4 for c in C15))
    require("all 30,720 single errors are corrected",
            all(decode(15, H15, flip(c, j)) == c for c in C15 for j in range(1, 16)))
    same("weight-3 codewords of the (15,11) code", sum(1 for c in C15 if weight(c) == 3), 35)


def check_family() -> None:
    section("The Hamming family")
    same("lengths", " ".join(str(n) for _, n, _ in FAMILY), "3 7 15 31 63 127")
    same("dimensions", " ".join(str(k) for _, _, k in FAMILY), "1 4 11 26 57 120")
    same("rates", " ".join(rate(k, n) for _, n, k in FAMILY), "0.333 0.571 0.733 0.839 0.905 0.945")
    require("2^k (n + 1) = 2^n for every member", all(2 ** k * (n + 1) == 2 ** n for _, n, k in FAMILY))
    same("the (3,1) code", show([c for c in all_codewords(3)]), "000 111")


def check_sphere_packing() -> None:
    section("The sphere-packing bound")
    same("words within distance 1 of a word of length 7", volume(7, 1), 8)
    same("16 balls of 8", len(C7) * volume(7, 1), 128)
    covered = [word(v) for c in C7 for v in ball(c, 1)]
    require("the 16 balls of radius 1 are disjoint", len(set(covered)) == len(covered))
    require("and cover all 128 words of length 7", len(set(covered)) == 128)
    same("words within distance 1 of a word of length 15", volume(15, 1), 16)
    same("2^11 times 16", 2 ** 11 * 16, 2 ** 15)
    same("C(23,1) C(23,2) C(23,3)", f"{comb(23, 1)} {comb(23, 2)} {comb(23, 3)}", "23 253 1771")
    same("1 + 23 + 253 + 1771", volume(23, 3), 2048)
    require("2^12 times 2^11 = 2^23", 2 ** 12 * volume(23, 3) == 2 ** 23)
    require("repetition codes of length 3 to 15: 2 V(n, (n - 1)/2) = 2^n",
            all(2 * volume(n, (n - 1) // 2) == 2 ** n for n in range(3, 16, 2)))
    same("V(5,2)", volume(5, 2), 16)
    same("1 + 90 + 4005", volume(90, 2), 4096)
    require("and 4096 = 2^12", volume(90, 2) == 2 ** 12)
    same("ternary Golay: 1 + 11*2 + 55*4", volume(11, 2, 3), 243)
    require("3^6 times 243 = 3^11", 3 ** 6 * volume(11, 2, 3) == 3 ** 11)


def check_ternary() -> None:
    section("Hamming codes over GF(3)")
    same("one-dimensional subspaces of GF(3)^2", (3 ** 2 - 1) // (3 - 1), 4)
    same("their representatives as columns", " ".join(word(c) for c in projective_columns(2, 3)), "01 10 11 12")
    same("H", show(H3), "0111 1012")
    code = ternary_codewords()
    same("codewords", len(code), 9)
    same("the nine codewords", " ".join(word(c) for c in code),
         "0000 0111 0222 1012 1120 1201 2021 2102 2210")
    require("every nonzero codeword has weight 3", all(weight(c) == 3 for c in code if any(c)))
    same("words within distance 1 of a ternary word of length 4", volume(4, 1, 3), 9)
    same("9 balls of 9", 9 * volume(4, 1, 3), 3 ** 4)
    rx = [1, 2, 1, 1]
    same("syndrome of 1211", word(syndrome_q(H3, rx, 3)), "11")
    pos, val, out = ternary_decode(rx)
    same("position, value, corrected word", f"{pos} {val} {word(out)}", "3 1 1201")
    rx = [1, 1, 0, 1]
    same("syndrome of 1101", word(syndrome_q(H3, rx, 3)), "20")
    pos, val, out = ternary_decode(rx)
    same("position, value, corrected word", f"{pos} {val} {word(out)}", "2 2 1201")
    ok = True
    for c in code:
        for j in range(4):
            for a in (1, 2):
                v = list(c)
                v[j] = (v[j] + a) % 3
                ok &= ternary_decode(v)[2] == list(c)
    require("all 72 single errors on the 9 codewords are corrected", ok)
    params = []
    for q, m in ((3, 2), (3, 3), (4, 2), (5, 2), (7, 2)):
        n = (q ** m - 1) // (q - 1)
        params.append((q, n, n - m))
        require(f"q = {q}, m = {m}: q^k (1 + n(q - 1)) = q^n", q ** (n - m) * volume(n, 1, q) == q ** n)
    same("(n, k) for q, m = 3 2, 3 3, 4 2, 5 2, 7 2",
         " ".join(f"({n},{k})" for _, n, k in params), "(4,2) (13,10) (5,3) (6,4) (8,6)")
    cols = projective_columns(3, 3)
    same("columns of the ternary code with m = 3", len(cols), 13)
    proportional = any(all((a * x) % 3 == y for x, y in zip(u, v))
                       for u, v in combinations(cols, 2) for a in (1, 2))
    require("no two of them are multiples of each other", not proportional)


def check_extended() -> None:
    section("The extended Hamming code")
    E8 = [extend(c) for c in C7]
    require("every extended (8,4) codeword has even weight", all(weight(c) % 2 == 0 for c in E8))
    same("minimum distance of the (8,4) code", min(weight(c) for c in E8 if any(c)), 4)
    same("weight-4 codewords of the (8,4) code", sum(1 for c in E8 if weight(c) == 4), 14)
    E16 = [extend(c) for c in C15]
    same("minimum distance of the (16,11) code", min(weight(c) for c in E16 if any(c)), 4)
    same("extended lengths and dimensions",
         " ".join(f"({n + 1},{k})" for _, n, k in FAMILY), "(4,1) (8,4) (16,11) (32,26) (64,57) (128,120)")
    single = all(extended_decision(flip_ext(c, j), H7) == expected_single(j)
                 for c in E8 for j in range(8))
    require("(8,4): every single error is named, parity bit included", single)
    single = all(extended_decision(flip_ext(c, j), H15) == expected_single(j)
                 for c in E16 for j in range(16))
    require("(16,11): every single error is named, parity bit included", single)
    doubles = list(combinations(range(8), 2))
    same("double-error patterns of length 8", len(doubles), 28)
    require("(8,4): every double error on every codeword is detected",
            all(extended_decision(flip_ext(c, i, j), H7) == "double error detected"
                for c in E8 for i, j in doubles))
    doubles = list(combinations(range(16), 2))
    same("double-error patterns of length 16", len(doubles), 120)
    require("(16,11): every double error on every codeword is detected",
            all(extended_decision(flip_ext(c, i, j), H15) == "double error detected"
                for c in E16 for i, j in doubles))
    x = extend(CW15)
    same("parity bit of 111001110011010", x[0], 1)
    same("extended codeword c0 c1 ... c15", word(x), "1111001110011010")
    rx = flip_ext(x, 3, 12)
    same("syndrome with positions 3 and 12 flipped", word(syndrome(H15, rx[1:])), "1111")
    same("overall parity of that word", sum(rx) % 2, 0)
    same("decision", extended_decision(rx, H15), "double error detected")
    wrong = decode(15, H15, rx[1:])
    same("positions where the plain decoder's output differs from c", weight([a ^ b for a, b in zip(wrong, CW15)]), 3)


def flip_ext(v: list, *positions) -> list:
    out = v[:]
    for p in positions:
        out[p] ^= 1
    return out


def expected_single(j: int) -> str:
    return "flip the parity bit" if j == 0 else f"flip position {j}"


def check_hardware() -> None:
    section("The encoder and decoder as gates")
    w, g, d = encoder_counts(7)
    same("(7,4) encoder: data bits per check, gates, depth", f"{' '.join(map(str, w))} {g} {d}", "3 3 3 6 2")
    w, g, d = encoder_counts(15)
    same("(15,11) encoder: gates and depth", f"{g} {d}", "24 3")
    gates = []
    depths = []
    for m in range(3, 8):
        _, g, d = encoder_counts(2 ** m - 1)
        require(f"m = {m}: gates equal m(2^(m - 1) - 2)", g == m * (2 ** (m - 1) - 2))
        gates.append(g)
        depths.append(d)
    same("encoder gates for m = 3 to 7", " ".join(map(str, gates)), "6 24 70 180 434")
    same("encoder depths for m = 3 to 7", " ".join(map(str, depths)), "2 3 4 5 6")
    for H, name, stated in ((H7, "(7,4)", "9 2 7 14 3 7"), (H15, "(15,11)", "28 3 15 45 4 15")):
        c = decoder_counts(H)
        same(f"{name} decoder: syndrome XORs, depth, ANDs, AND2s, NOTs, XORs",
             f"{c['syndrome']} {c['syndrome depth']} {c['and']} {c['and2']} {c['not']} {c['xor']}", stated)
    for H, name, total, depth in ((H7, "(7,4)", 33, 5), (H15, "(15,11)", 92, 6)):
        c = decoder_counts(H)
        same(f"{name} decoder: two-input gates in all", c["syndrome"] + c["and2"] + c["not"] + c["xor"], total)
        same(f"{name} decoder: gate levels, tree + NOT + AND + XOR", c["syndrome depth"] + 3, depth)
    for H, code, n in ((H7, C7, 7), (H15, C15, 15)):
        ok = True
        for c in code:
            lines, out = gate_decoder(H, c)
            ok &= sum(lines) == 0 and out == c
            for j in range(1, n + 1):
                lines, out = gate_decoder(H, flip(c, j))
                ok &= sum(lines) == 1 and lines[j - 1] == 1 and out == c
        require(f"gate-level decoder of length {n} on every correctable word", ok)
    same("data-only decoder ANDs and XORs, (7,4) and (15,11)", f"{len(data_positions(7))} {len(data_positions(15))}", "4 11")


def check_examples() -> None:
    section("Worked examples")
    H12 = columns_H(12, 4)
    same("H of the (12,8) code, columns 1 to 12", show(H12), "000000011111 000111100001 011001100110 101010101010")
    c = encode(12, parse("11001010"))
    same("byte 11001010 encodes to", word(c), "001110001010")
    require("its syndrome is zero", syndrome(H12, c) == [0] * 4)
    code12 = all_codewords(12)
    same("(12,8) codewords and minimum distance", f"{len(code12)} {min(weight(x) for x in code12 if any(x))}", "256 3")
    same("(7,4): data read from 0110011", word(extract(7, decode(7, H7, parse("0110001")))), "1011")
    same("V(10,2)", volume(10, 2), 56)
    same("floor(1024 / 56)", 1024 // 56, 18)
    same("largest k with 2^k <= 18", max(k for k in range(11) if 2 ** k <= 18), 4)
    same("five data bits in 8 bits: 32 times 9", 32 * volume(8, 1), 288)
    same("five data bits in 9 bits: 32 times 10", 32 * volume(9, 1), 320)
    x = extend(encode(7, parse("1011")))
    same("extended codeword of 1011", word(x), "00110011")
    rx = flip_ext(x, 2, 5)
    same("positions 2 and 5 flipped", word(rx), "00010111")
    same("its syndrome and overall parity", f"{word(syndrome(H7, rx[1:]))} {sum(rx) % 2}", "111 0")
    wrong = decode(7, H7, rx[1:])
    same("the (7,4) decoder alone outputs", word(wrong), "0010110")
    require("which is a codeword", syndrome(H7, wrong) == [0, 0, 0])
    same("carrying the data", word(extract(7, wrong)), "1110")


def check_exercises() -> None:
    section("Exercises")
    same("least m for k = 26 and k = 27", f"{min_redundancy(26)} {min_redundancy(27)}", "5 6")
    same("C(15,2) / 3", comb(15, 2) // 3, 35)
    rx = flip(CW15, 6)
    same("received word", word(rx), "111000110011010")
    same("its syndrome", word(syndrome(H15, rx)), "0110")
    same("H for m = 2", show(hamming_H(2)), "011 101")
    same("C(90,2)", comb(90, 2), 4005)
    n3 = (3 ** 3 - 1) // (3 - 1)
    same("ternary m = 3: n, k, 1 + 2n", f"{n3} {n3 - 3} {1 + 2 * n3}", "13 10 27")
    _, g, d = encoder_counts(31)
    same("(31,26) encoder gates and depth", f"{g} {d}", "70 4")
    rx = parse("01110011")
    same("syndrome and parity of 01110011", f"{word(syndrome(H7, rx[1:]))} {sum(rx) % 2}", "001 1")
    same("decision", extended_decision(rx, H7), "flip position 1")


def main() -> None:
    check_counting()
    check_hamming_7_4()
    check_distance()
    check_hamming_15_11()
    check_family()
    check_sphere_packing()
    check_ternary()
    check_extended()
    check_hardware()
    check_examples()
    check_exercises()
    print("\nAll chapter 19 numbers reproduced.")


if __name__ == "__main__":
    main()
