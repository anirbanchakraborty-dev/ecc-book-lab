"""Reproduce every number in chapter 24 of Error Correction Coding, "Modifying and Combining Codes".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch24_modifying_and_combining_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in this chapter is a whole number, a binary word, a binary matrix,
a weight distribution or a rate printed to three places, so every comparison is
exact. Every code is small enough to list in full, so each minimum distance is
found by computing the weight of every codeword rather than taken from a
theorem. The largest lists are the 65,536 codewords of the (25,16) rectangular
parity code and of the (49,16) product of the Hamming code with itself.

A word is a tuple of bits, a matrix is a list of rows, and a code is the sorted
list of its distinct codewords. The module also holds the operations that
ch24_charts.py draws, so the figures and this check share one source.
"""

import sys
from itertools import product


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


# Words, matrices and codes over GF(2) ---------------------------------------------

def parse(text: str) -> list:
    """Rows written as digit strings separated by spaces, such as "1101 1011"."""
    return [tuple(int(ch) for ch in row) for row in text.split()]


def show(M: list) -> str:
    return " ".join("".join(str(v) for v in row) for row in M)


def word(v) -> str:
    return "".join(str(b) for b in v)


def weight(v) -> int:
    return sum(v)


def xor(a, b) -> tuple:
    return tuple(x ^ y for x, y in zip(a, b))


def mat_mul(A: list, B: list) -> list:
    """The product of two binary matrices, reduced modulo 2."""
    return [tuple(sum(a * b for a, b in zip(row, col)) % 2 for col in zip(*B)) for row in A]


def transpose(A: list) -> list:
    return [tuple(col) for col in zip(*A)]


def rank(A: list) -> int:
    """Rank over GF(2) by elimination."""
    rows = [list(r) for r in A]
    r = 0
    for c in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [x ^ y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r


def span(G: list) -> list:
    """Every codeword m G, in message order, as tuples."""
    n = len(G[0])
    out = []
    for m in product([0, 1], repeat=len(G)):
        c = (0,) * n
        for bit, row in zip(m, G):
            if bit:
                c = xor(c, row)
        out.append(c)
    return out


def code(G: list) -> list:
    return sorted(set(span(G)))


def wdist(C: list) -> str:
    counts = {}
    for c in C:
        counts[weight(c)] = counts.get(weight(c), 0) + 1
    return " ".join(f"{w}:{counts[w]}" for w in sorted(counts))


def dmin(C: list) -> int:
    return min(weight(c) for c in C if any(c))


def dimension(C: list) -> int:
    return len(C).bit_length() - 1


def params(C: list) -> str:
    return f"({len(C[0])},{dimension(C)},{dmin(C)})"


def puncture(C: list, i: int) -> list:
    """Delete position i, counted from 1."""
    return sorted(set(c[:i - 1] + c[i:] for c in C))


def shorten(C: list, i: int) -> list:
    """Keep the codewords that are 0 at position i, counted from 1, then delete it."""
    return sorted(set(c[:i - 1] + c[i:] for c in C if c[i - 1] == 0))


def extend(C: list) -> list:
    return sorted(c + (weight(c) % 2,) for c in C)


def expurgate(C: list) -> list:
    return sorted(c for c in C if weight(c) % 2 == 0)


def augment(C: list) -> list:
    ones = (1,) * len(C[0])
    return sorted(set(C) | {xor(c, ones) for c in C})


def lengthen(C: list) -> list:
    return extend(augment(C))


def direct_sum(C1: list, C2: list) -> list:
    return sorted(u + v for u in C1 for v in C2)


def u_u_plus_v(C1: list, C2: list) -> list:
    return sorted(u + xor(u, v) for u in C1 for v in C2)


def product_array(G1: list, G2: list, M: list) -> list:
    """The n1-by-n2 array G1^T M G2: rows encoded by the row code C2, columns by C1."""
    return mat_mul(mat_mul(transpose(G1), M), G2)


def rows_then_columns(G1: list, G2: list, M: list) -> list:
    """Encode the k1 data rows with C2, then every column with C1."""
    rows = mat_mul(M, G2)
    return transpose(mat_mul(transpose(rows), G1))


def columns_then_rows(G1: list, G2: list, M: list) -> list:
    """Encode the k2 data columns with C1, then every row with C2."""
    cols = mat_mul(transpose(G1), M)
    return mat_mul(cols, G2)


def product_generator(G1: list, G2: list) -> list:
    """Rows are the outer products of a row of G1 and a row of G2, read row by row."""
    return [tuple(a * b for a in g1 for b in g2) for g1 in G1 for g2 in G2]


def big_wdist(G: list) -> dict:
    """The weight distribution of a larger code, with words held as integers."""
    words = [0]
    for row in G:
        g = int(word(row), 2)
        words += [w ^ g for w in words]
    counts = {}
    for w in words:
        b = bin(w).count("1")
        counts[b] = counts.get(b, 0) + 1
    return counts


def gates(checks: list) -> int:
    """XOR gates for a set of parity equations, one tree per distinct equation."""
    return sum(max(len(eq) - 1, 0) for eq in set(frozenset(e) for e in checks))


def check_positions(G: list) -> list:
    """Positions (from 1) whose column of a systematic G is not a unit vector, so they hold checks."""
    return [j + 1 for j, col in enumerate(zip(*G)) if sum(col) != 1]


def check_equations(G: list, positions: list, kept=None) -> list:
    """The message bits (from 1) that each check position sums, read down its column of G.
    Message bits outside kept are tied to 0 and drop out of every equation."""
    return [[r + 1 for r in range(len(G)) if G[r][j - 1] and (kept is None or r + 1 in kept)]
            for j in positions]


def rate(k: int, n: int) -> str:
    return f"{k / n:.3f}"


# The codes of the chapter ---------------------------------------------------------

G7 = parse("1000110 0100011 0010101 0001111")
H7 = parse("1011100 1101010 0111001")
C7 = code(G7)
E8 = extend(C7)
SIMPLEX = code(H7)
EVEN4 = code(parse("1001 0101 0011"))
REP4 = code(parse("1111"))
REP3 = code(parse("111"))
G32 = parse("101 011")
G43 = parse("1001 0101 0011")
D7 = code(parse("1110000 0001111"))
S6 = shorten(C7, 1)
S5 = shorten(S6, 1)
DIRECT8 = direct_sum(EVEN4, REP4)
UUV8 = u_u_plus_v(EVEN4, REP4)
PRODUCT9 = sorted(tuple(x for row in product_array(G32, G32, [m[:2], m[2:]]) for x in row)
                  for m in product([0, 1], repeat=4))

# Every code the rate-against-distance chart plots: label, code.
CHART_CODES = [
    ("Hamming", C7),
    ("extended", E8),
    ("punctured", puncture(C7, 7)),
    ("shortened", S6),
    ("twice shortened", S5),
    ("expurgated", expurgate(C7)),
    ("direct sum", DIRECT8),
    ("(u | u + v)", UUV8),
    ("product", PRODUCT9),
]


# Checks -------------------------------------------------------------------------

def check_hamming_code() -> None:
    section("The (7,4) Hamming code of Linear Codes")
    same("generator G", show(G7), "1000110 0100011 0010101 0001111")
    same("parity-check matrix H", show(H7), "1011100 1101010 0111001")
    require("H G^T = 0", all(v == 0 for row in mat_mul(H7, transpose(G7)) for v in row))
    same("number of codewords", len(C7), 16)
    same("weight distribution", wdist(C7), "0:1 3:7 4:7 7:1")
    same("minimum distance", dmin(C7), 3)
    require("the all-ones word 1111111 is a codeword", (1,) * 7 in C7)
    same("weights of the rows of H", " ".join(str(weight(h)) for h in H7), "4 4 4")
    same("1s shared by rows 1 and 2, 1 and 3, 2 and 3 of H",
         " ".join(str(sum(x & y for x, y in zip(H7[i], H7[j]))) for i, j in [(0, 1), (0, 2), (1, 2)]), "2 2 2")
    require("H H^T = 0", all(v == 0 for row in mat_mul(H7, transpose(H7)) for v in row))
    same("rank of H", rank(H7), 3)


def check_memory_word() -> None:
    section("A 64-bit memory word")
    same("Hamming code with m = 6 checks, (n,k)", f"({2**6 - 1},{2**6 - 1 - 6})", "(63,57)")
    same("Hamming code with m = 7 checks, (n,k)", f"({2**7 - 1},{2**7 - 1 - 7})", "(127,120)")
    same("data positions to remove from 120", 120 - 64, 56)
    same("shortened Hamming code, (n,k)", f"({127 - 56},{120 - 56})", "(71,64)")
    same("extended shortened code, (n,k)", f"({127 - 56 + 1},{120 - 56})", "(72,64)")
    same("extended Hamming code, then shortened, (n,k)", f"({128 - 56},{120 - 56})", "(72,64)")
    same("check bits of the (72,64) code", 72 - 64, 8)
    same("C(72,2)", 72 * 71 // 2, 2556)
    same("error patterns of weight at most 2 in 72 bits", 1 + 72 + 72 * 71 // 2, 2629)
    same("syndromes of 8 check bits", 2**8, 256)
    require("2629 patterns exceed 256 syndromes, so d is at most 4", 1 + 72 + 72 * 71 // 2 > 2**8)
    same("rate of the (127,120) code", rate(120, 127), "0.945")
    same("rate of the (72,64) code", rate(64, 72), "0.889")


def check_extending() -> None:
    section("Extending")
    GE = [g + (weight(g) % 2,) for g in G7]
    HE = [h + (0,) for h in H7] + [(1,) * 8]
    same("extended generator", show(GE), "10001101 01000111 00101011 00011110")
    same("extended parity-check matrix", show(HE), "10111000 11010100 01110010 11111111")
    require("extended H times extended G^T = 0",
            all(v == 0 for row in mat_mul(HE, transpose(GE)) for v in row))
    same("rank of the extended parity-check matrix", rank(HE), 4)
    same("extended code, weight distribution", wdist(E8), "0:1 4:14 8:1")
    same("extended code, minimum distance", dmin(E8), 4)
    require("every extended codeword has even weight", all(weight(c) % 2 == 0 for c in E8))
    require("the extended code equals the span of the extended generator", code(GE) == E8)
    same("weight-3 codeword 1000110 extended", word(GE[0]), "10001101")
    same("overall parity bit as message coefficients", word(row[7] for row in GE), "1110")
    require("c5 XOR c6 XOR c7 = m4 for every message",
            all(c[4] ^ c[5] ^ c[6] == m[3] for m, c in zip(product([0, 1], repeat=4), span(G7))))
    require("overall parity equals m1 XOR m2 XOR m3 for every message",
            all(weight(c) % 2 == m[0] ^ m[1] ^ m[2] for m, c in zip(product([0, 1], repeat=4), span(G7))))
    same("XOR gates in the overall parity tree", gates(check_equations(GE, [8])), 2)
    GEV = [g + (weight(g) % 2,) for g in G43]
    same("even-weight (4,3) code, extended generator", show(GEV), "10010 01010 00110")
    require("every appended bit of the even-weight code is 0", all(c[4] == 0 for c in extend(EVEN4)))
    same("even-weight code, distance before and after", f"{dmin(EVEN4)} {dmin(extend(EVEN4))}", "2 2")


def check_puncturing() -> None:
    section("Puncturing")
    require("puncturing the extended code at position 8 gives the (7,4) code", puncture(E8, 8) == C7)
    P6 = puncture(C7, 7)
    same("(7,4) punctured at position 7, (n,k,d)", params(P6), "(6,4,2)")
    same("its weight distribution", wdist(P6), "0:1 2:3 3:8 4:3 6:1")
    same("weight-2 word from the row 0100011", word(G7[1][:6]), "010001")
    three = [c for c in C7 if weight(c) == 3]
    same("weight-3 codewords covering each position",
         " ".join(str(sum(c[i] for c in three)) for i in range(7)), "3 3 3 3 3 3 3")
    same("code D, weight distribution", wdist(D7), "0:1 3:1 4:1 7:1")
    same("code D, (n,k,d)", params(D7), "(7,2,3)")
    PD = puncture(D7, 7)
    same("code D punctured at position 7, (n,k,d)", params(PD), "(6,2,3)")
    same("its weight distribution", wdist(PD), "0:1 3:2 6:1")


def check_shortening() -> None:
    section("Shortening")
    G6 = parse("100011 010101 001111")
    H6 = [h[1:] for h in H7]
    same("(6,3) generator", show(G6), "100011 010101 001111")
    require("the (6,3) code is the span of its generator", code(G6) == S6)
    same("(6,3) parity-check matrix, column 1 of H deleted", show(H6), "011100 101010 111001")
    require("(6,3) H times G^T = 0", all(v == 0 for row in mat_mul(H6, transpose(G6)) for v in row))
    require("columns of the (6,3) H are distinct and nonzero",
            len(set(transpose(H6))) == 6 and (0, 0, 0) not in transpose(H6))
    same("(6,3) codewords in message order", show(span(G6)),
         "000000 001111 010101 011010 100011 101100 110110 111001")
    same("(6,3) weight distribution", wdist(S6), "0:1 3:4 4:3")
    same("(6,3) code, (n,k,d)", params(S6), "(6,3,3)")
    G5 = parse("10101 01111")
    H5 = [h[2:] for h in H7]
    require("the (5,2) code is the span of its generator", code(G5) == S5)
    same("(5,2) parity-check matrix, columns 1 and 2 deleted", show(H5), "11100 01010 11001")
    require("(5,2) H times G^T = 0", all(v == 0 for row in mat_mul(H5, transpose(G5)) for v in row))
    same("(5,2) weight distribution", wdist(S5), "0:1 3:2 4:1")
    same("(5,2) code, (n,k,d)", params(S5), "(5,2,3)")
    require("(5,2) check bits 3 and 5 are equal", all(c[2] == c[4] for c in S5))
    checks7 = check_positions(G7)
    same("check positions of the (7,4) generator", " ".join(map(str, checks7)), "5 6 7")
    same("encoder XOR gates, (7,4)", gates(check_equations(G7, checks7)), 6)
    same("encoder XOR gates, (6,3) with m1 tied to 0", gates(check_equations(G7, checks7, {2, 3, 4})), 4)
    same("encoder XOR gates, (5,2) with m1 and m2 tied to 0", gates(check_equations(G7, checks7, {3, 4})), 1)
    c5, c6, c7 = map(set, check_equations(G7, checks7, {2, 3, 4}))
    require("(6,3) check c7 = c5 XOR m2", c7 == c5 ^ {2})
    same("encoder XOR gates, (6,3) with c7 built as c5 XOR m2", gates([c5, c6]) + 1, 3)


def check_expurgating() -> None:
    section("Expurgating")
    X = expurgate(C7)
    same("even-weight codewords of the (7,4) code", len(X), 8)
    same("expurgated code, weight distribution", wdist(X), "0:1 4:7")
    same("expurgated code, (n,k,d)", params(X), "(7,3,4)")
    require("the expurgated code equals the row space of H", X == SIMPLEX)
    require("a codeword has even weight exactly when m1 XOR m2 XOR m3 = 0",
            all((weight(c) % 2 == 0) == (m[0] ^ m[1] ^ m[2] == 0)
                for m, c in zip(product([0, 1], repeat=4), span(G7))))
    same("messages kept by expurgating", sum(1 for m in product([0, 1], repeat=4) if m[0] ^ m[1] ^ m[2] == 0), 8)
    odd = [r + 1 for r, g in enumerate(G7) if weight(g) % 2]
    same("message bits of the odd-weight rows of G", " ".join(map(str, odd)), "1 2 3")
    same("XOR gates for m3 = m1 XOR m2", gates([odd[:-1]]), 1)
    a, b = parse("1000110 0100011")
    same("w(a), w(b), w(a AND b), w(a XOR b)",
         f"{weight(a)} {weight(b)} {sum(x & y for x, y in zip(a, b))} {weight(xor(a, b))}", "3 3 1 4")
    words7 = list(product([0, 1], repeat=7))
    require("w(a XOR b) = w(a) + w(b) - 2 w(a AND b) for all words of length 7",
            all(weight(xor(p, q)) == weight(p) + weight(q) - 2 * sum(x & y for x, y in zip(p, q))
                for p in words7 for q in words7))


def check_augmenting() -> None:
    section("Augmenting")
    wmax = max(weight(c) for c in SIMPLEX)
    same("largest weight in the (7,3) code", wmax, 4)
    require("1111111 is not in the (7,3) code", (1,) * 7 not in SIMPLEX)
    A = augment(SIMPLEX)
    require("augmenting the (7,3) code gives the (7,4) code", A == C7)
    same("min(d, n - w_max) for the (7,3) code", min(dmin(SIMPLEX), 7 - wmax), 3)
    same("augmented code, minimum distance", dmin(A), 3)
    wmax5 = max(weight(c) for c in S5)
    same("largest weight in the (5,2) code", wmax5, 4)
    same("min(d, n - w_max) for the (5,2) code", min(dmin(S5), 5 - wmax5), 1)
    same("augmented (5,2) code, (n,k,d)", params(augment(S5)), "(5,3,1)")
    same("augmented (5,2) code, weight distribution", wdist(augment(S5)), "0:1 1:1 2:2 3:2 4:1 5:1")


def check_lengthening() -> None:
    section("Lengthening")
    L = lengthen(SIMPLEX)
    same("lengthened (7,3) code, (n,k,d)", params(L), "(8,4,4)")
    require("lengthening the (7,3) code gives the extended (8,4) code", L == E8)


def check_direct_sum() -> None:
    section("Direct sum")
    same("even-weight (4,3) code, weight distribution", wdist(EVEN4), "0:1 2:6 4:1")
    same("repetition (4,1) code, weight distribution", wdist(REP4), "0:1 4:1")
    same("direct sum, (n,k,d)", params(DIRECT8), "(8,4,2)")
    same("direct sum, weight distribution", wdist(DIRECT8), "0:1 2:6 4:2 6:6 8:1")
    a = [1, 0, 6, 0, 1]
    b = [1, 0, 0, 0, 1]
    poly = [sum(a[i] * b[k - i] for i in range(len(a)) if 0 <= k - i < len(b)) for k in range(9)]
    same("coefficients of (1 + 6x^2 + x^4)(1 + x^4)", " ".join(map(str, poly)), "1 0 6 0 2 0 6 0 1")
    require("code D is the direct sum of the (3,1) and (4,1) repetition codes", direct_sum(REP3, REP4) == D7)
    same("rate of the direct sum", rate(4, 8), "0.500")


def check_u_u_plus_v() -> None:
    section("The (u | u+v) construction")
    G = parse("10011001 01010101 00110011 00001111")
    same("(u | u+v) generator", show(G), "10011001 01010101 00110011 00001111")
    require("the generator spans the (u | u+v) code", code(G) == UUV8)
    same("(u | u+v) code, (n,k,d)", params(UUV8), "(8,4,4)")
    same("min(2 d1, d2)", min(2 * dmin(EVEN4), dmin(REP4)), 4)
    same("(u | u+v) weight distribution", wdist(UUV8), "0:1 4:14 8:1")
    require("its weight distribution equals the extended code's", wdist(UUV8) == wdist(E8))
    same("upper-bound word (u, u) with u = 1001", word((1, 0, 0, 1) * 2), "10011001")
    same("upper-bound word (0, v) with v = 1111", word((0,) * 4 + (1,) * 4), "00001111")
    v_zero = [u + u for u in EVEN4 if any(u)]
    v_nonzero = [u + xor(u, v) for u in EVEN4 for v in REP4 if any(v)]
    same("least weight with v = 0, then with v nonzero",
         f"{min(map(weight, v_zero))} {min(map(weight, v_nonzero))}", "4 4")
    order = [1, 2, 3, 5, 7, 6, 8, 4]
    require("positions read as 1 2 3 5 7 6 8 4 give the extended code",
            sorted(tuple(c[i - 1] for i in order) for c in UUV8) == E8)
    u, v = (0, 1, 1, 0), (1, 1, 1, 1)
    same("encoding u = 0110, v = 1111", word(u + xor(u, v)), "01101001")
    require("01101001 is a codeword", u + xor(u, v) in UUV8)


def check_product() -> None:
    section("Product codes")
    M = [(1, 0), (1, 1)]
    X = product_array(G32, G32, M)
    same("(3,2) x (3,2) array for data 10 11", show(X), "101 110 011")
    require("rows first and columns first agree on all 16 data arrays",
            all(rows_then_columns(G32, G32, [m[:2], m[2:]]) == columns_then_rows(G32, G32, [m[:2], m[2:]])
                for m in product([0, 1], repeat=4)))
    even3 = code(G32)
    require("every row and every column of every array is even",
            all(tuple(r) in even3 and tuple(c) in even3
                for m in product([0, 1], repeat=4)
                for A in [product_array(G32, G32, [m[:2], m[2:]])]
                for r, c in zip(A, transpose(A))))
    same("(3,2) x (3,2) product, (n,k,d)", params(PRODUCT9), "(9,4,4)")
    same("(3,2) x (3,2) product, weight distribution", wdist(PRODUCT9), "0:1 4:9 6:6")
    same("outer product of 110 with 110", show([tuple(a * b for b in (1, 1, 0)) for a in (1, 1, 0)]),
         "110 110 000")
    G54 = parse("10001 01001 00101 00011")
    G25 = product_generator(G54, G54)
    counts = big_wdist(G25)
    same("rectangular parity, (n,k) from the generator's rank", f"({len(G25[0])},{rank(G25)})", "(25,16)")
    require("rectangular parity, codewords listed once each", counts[0] == 1)
    same("rectangular parity, minimum distance", min(w for w in counts if w), 4)
    same("rectangular parity, codewords of weight 4", counts[4], 100)
    eq32 = check_equations(G32, check_positions(G32))
    same("product encoder XOR gates, row and column encoders reused",
         len(G32) * gates(eq32) + len(G32[0]) * gates(eq32), 5)
    G9 = product_generator(G32, G32)
    same("check positions of the (9,4) product generator", " ".join(map(str, check_positions(G9))), "3 6 7 8 9")
    same("product encoder XOR gates, one tree per check", gates(check_equations(G9, check_positions(G9))), 7)


def check_worked_examples() -> None:
    section("Worked examples")
    require("shortening the extended code at position 8 gives the (7,3) code", shorten(E8, 8) == SIMPLEX)
    same("extended code shortened at position 8, (n,k,d)", params(shorten(E8, 8)), "(7,3,4)")
    same("extended code punctured at position 8, (n,k,d)", params(puncture(E8, 8)), "(7,4,3)")
    PR = sorted(tuple(x for row in product_array(G43, G32, [m[0:2], m[2:4], m[4:6]]) for x in row)
                for m in product([0, 1], repeat=6))
    same("(4,3) x (3,2) product, (n,k,d)", params(PR), "(12,6,4)")
    same("(4,3) x (3,2) array for data 10 11 00", show(product_array(G43, G32, [(1, 0), (1, 1), (0, 0)])),
         "101 110 000 011")
    same("weight of that array", sum(map(sum, product_array(G43, G32, [(1, 0), (1, 1), (0, 0)]))), 6)
    same("rate of the (12,6) product", rate(6, 12), "0.500")
    require("111111 is not in the (6,3) code", (1,) * 6 not in S6)
    wmax6 = max(weight(c) for c in S6)
    same("largest weight in the (6,3) code", wmax6, 4)
    same("min(d, n - w_max) for the (6,3) code", min(dmin(S6), 6 - wmax6), 2)
    same("augmented (6,3) code, (n,k,d)", params(augment(S6)), "(6,4,2)")
    same("lengthened (6,3) code, (n,k,d)", params(lengthen(S6)), "(7,4,2)")
    require("110110 is a (6,3) codeword", parse("110110")[0] in S6)
    same("111111 XOR 110110", word(xor((1,) * 6, parse("110110")[0])), "001001")
    same("its weight", weight(xor((1,) * 6, parse("110110")[0])), 2)


def check_rate_and_distance() -> None:
    section("Rate and distance together")
    same("rates of the chart codes", " ".join(rate(dimension(C), len(C[0])) for _, C in CHART_CODES),
         "0.571 0.500 0.667 0.500 0.400 0.429 0.500 0.500 0.444")
    same("distances of the chart codes", " ".join(str(dmin(C)) for _, C in CHART_CODES), "3 4 2 3 3 4 2 4 4")
    same("(n,k) of the chart codes", " ".join(f"({len(C[0])},{dimension(C)})" for _, C in CHART_CODES),
         "(7,4) (8,4) (6,4) (6,3) (5,2) (7,3) (8,4) (8,4) (9,4)")


def check_exercises() -> None:
    section("Exercises")
    same("extended (6,3) code, (n,k,d)", params(extend(S6)), "(7,3,4)")
    same("extended (6,3) code, weight distribution", wdist(extend(S6)), "0:1 4:7")
    PC = puncture(PRODUCT9, 9)
    same("(9,4) product punctured at the corner, (n,k,d)", params(PC), "(8,4,3)")
    same("weight-4 product codewords through the corner", sum(1 for c in PRODUCT9 if weight(c) == 4 and c[8]), 4)
    same("extended code shortened at position 1, (n,k,d)", params(shorten(E8, 1)), "(7,3,4)")
    require("extend then shorten equals shorten then extend at position 1",
            shorten(extend(C7), 1) == extend(shorten(C7, 1)))
    dual = puncture(SIMPLEX, 1)
    require("the (6,3) shortened code is the dual of the punctured (7,3) code",
            all(sum(x & y for x, y in zip(p, q)) % 2 == 0 for p in S6 for q in dual)
            and dimension(S6) + dimension(dual) == 6)
    EVEN8 = code([tuple(1 if j in (i, 7) else 0 for j in range(8)) for i in range(7)])
    same("even-weight (8,7) code, (n,k,d)", params(EVEN8), "(8,7,2)")
    same("direct sum of (8,7,2) and (8,4,4), (n,k,d)", params(direct_sum(EVEN8, E8)), "(16,11,2)")
    same("(u | u+v) of (8,7,2) and (8,4,4), (n,k,d)", params(u_u_plus_v(EVEN8, E8)), "(16,11,4)")
    G49 = product_generator(G7, G7)
    counts = big_wdist(G49)
    same("(7,4) x (7,4) product, (n,k) from the generator's rank", f"({len(G49[0])},{rank(G49)})", "(49,16)")
    require("(7,4) x (7,4) product, codewords listed once each", counts[0] == 1)
    same("(7,4) x (7,4) product, minimum distance", min(w for w in counts if w), 9)
    same("(7,4) x (7,4) product, codewords of weight 9", counts[9], 49)
    same("rate of the (49,16) product", rate(rank(G49), len(G49[0])), "0.327")


def main() -> None:
    check_hamming_code()
    check_memory_word()
    check_extending()
    check_puncturing()
    check_shortening()
    check_expurgating()
    check_augmenting()
    check_lengthening()
    check_direct_sum()
    check_u_u_plus_v()
    check_product()
    check_rate_and_distance()
    check_worked_examples()
    check_exercises()
    print("\nAll chapter 24 numbers reproduced.")


if __name__ == "__main__":
    main()
