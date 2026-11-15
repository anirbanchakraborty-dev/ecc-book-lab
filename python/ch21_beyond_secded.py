"""Rebuild the codes of chapter 21 of Error Correction Coding, "Beyond SECDED", and recompute its numbers.

    python3 python/ch21_beyond_secded.py

It needs Python 3.10 or later and nothing outside the standard library. Each
check prints the value computed here beside the value the chapter states, and
the program stops with exit status 1 at the first pair that disagrees.

An element of GF(2^m) is held as a whole number whose bit j is the coefficient
of alpha^j, and it is printed with the highest coefficient first, so alpha^4 =
alpha + 1 in GF(16) prints as 0011. A column of a binary check matrix with r
rows is a whole number too, whose highest bit is row 0, so the printed string
is the column read from top to bottom. Positions are counted from 0. Every search
in the file runs in a fixed order, so it finds the same code on every run and
on the Octave route.
"""

import sys
from functools import lru_cache
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


def popcount(x: int) -> int:
    return bin(x).count("1")


def bits(v: int, width: int) -> str:
    return format(v, f"0{width}b")


def pct(part: int, whole: int) -> str:
    """A share in percent to two decimals, trailing zeros dropped."""
    return f"{100 * part / whole:.2f}".rstrip("0").rstrip(".")


def pct1(part: int, whole: int) -> str:
    return f"{100 * part / whole:.1f}"


# Binary polynomials and finite fields -----------------------------------------

def poly_deg(p: int) -> int:
    return p.bit_length() - 1


def poly_rem(a: int, b: int) -> int:
    db = poly_deg(b)
    while a and poly_deg(a) >= db:
        a ^= b << (poly_deg(a) - db)
    return a


def irreducible(p: int) -> bool:
    """No factor of degree 1 to deg(p) // 2, found by trial division."""
    m = poly_deg(p)
    for d in range(2, 1 << (m // 2 + 1)):
        if poly_rem(p, d) == 0:
            return False
    return True


def poly_mulmod(a: int, b: int, poly: int, m: int) -> int:
    """Shift-and-add product of two field elements, reduced by poly."""
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if (a >> m) & 1:
            a ^= poly
    return r


def power_table(poly: int, m: int) -> tuple:
    q1 = (1 << m) - 1
    exp = [0] * q1
    log = [-1] * (1 << m)
    v = 1
    for i in range(q1):
        exp[i] = v
        log[v] = i
        v = poly_mulmod(v, 2, poly, m)
    return exp, log, v


def gmul(a: int, b: int, exp: list, log: list) -> int:
    if a == 0 or b == 0:
        return 0
    return exp[(log[a] + log[b]) % len(exp)]


def gdiv(a: int, b: int, exp: list, log: list) -> int:
    if a == 0:
        return 0
    return exp[(log[a] - log[b]) % len(exp)]


def rank_gf2(cols: list) -> int:
    basis = {}
    rank = 0
    for c in cols:
        x = c
        while x:
            h = x.bit_length() - 1
            if h in basis:
                x ^= basis[h]
            else:
                basis[h] = x
                rank += 1
                break
    return rank


GF16_POLY = 0b10011        # x^4 + x + 1
GF128_POLY = 0b10000011    # x^7 + x + 1
GF16 = power_table(GF16_POLY, 4)
GF128 = power_table(GF128_POLY, 7)


def bch_columns(m: int, field: tuple, n: int) -> list:
    """Columns (beta, beta^3) for beta = alpha^0, ..., alpha^(n-1), beta in the top m rows."""
    exp = field[0]
    q1 = len(exp)
    return [(exp[i % q1] << m) | exp[(3 * i) % q1] for i in range(n)]


def ball(n: int, t: int) -> int:
    return sum(comb(n, i) for i in range(t + 1))


def least_check_bits(k: int, t: int, extended: bool) -> int:
    """Smallest r that the sphere-packing count allows for t-error correction of k data bits.

    With extended set, the code also has distance 2t + 2, and deleting one
    coordinate leaves 2^k words of length k + r - 1 at distance 2t + 1."""
    r = 1
    while True:
        shift = 1 if extended else 0
        if (1 << (r - shift)) >= ball(k + r - shift, t):
            return r
        r += 1


def hamming_columns() -> list:
    """The (71,64) Hamming code: the first 71 nonzero 7-bit words, all distinct."""
    return list(range(1, 72))


def hsiao_columns() -> list:
    """A (72,64) Hsiao code: the 8 unit words, all 56 words of weight 3 and the
    first 8 words of weight 5, all of odd weight and all distinct."""
    def weight(w):
        return [v for v in range(256) if popcount(v) == w]
    return weight(1) + weight(3) + weight(5)[:8]


def check_baselines() -> None:
    section("The Hamming SEC and Hsiao SEC-DED codes for 64 data bits, built")
    ham = hamming_columns()
    require("Hamming (71,64): 71 distinct nonzero columns", len(set(ham)) == 71 and 0 not in ham)
    same("Hamming (71,64): rank of the 7 x 71 check matrix", rank_gf2(ham), 7)
    hs = hsiao_columns()
    require("Hsiao (72,64): 72 distinct columns of odd weight", len(set(hs)) == 72 and all(popcount(c) % 2 for c in hs))
    same("Hsiao (72,64): rank of the 8 x 72 check matrix", rank_gf2(hs), 8)
    cols = set(hs)
    require("Hsiao (72,64): every double has a nonzero syndrome no single has",
            all(hs[i] ^ hs[j] and hs[i] ^ hs[j] not in cols for i in range(72) for j in range(i + 1, 72)))


# The opening problem ----------------------------------------------------------

def check_opening() -> None:
    section("One strike, two neighboring bits of a (72,64) word")
    same("double errors in a 72-bit word", comb(72, 2), 2556)
    same("adjacent pairs in a 72-bit word", sum(1 for i in range(72) for j in range(i + 1, 72) if j - i == 1), 71)
    same("check bits of the (72,64) SEC-DED code", rank_gf2(hsiao_columns()), 8)


# DEC-TED ----------------------------------------------------------------------

STATED_GF16 = ["0001", "0010", "0100", "1000", "0011", "0110", "1100", "1011",
               "0101", "1010", "0111", "1110", "1111", "1101", "1001"]


def check_sphere_counts() -> None:
    section("Sphere-packing minimum for 64 data bits")
    same("ball of radius 1 in 71 bits", ball(71, 1), 72)
    same("ball of radius 1 in 70 bits", ball(70, 1), 71)
    same("least check bits for SEC", least_check_bits(64, 1, False), 7)
    same("least check bits for SEC-DED", least_check_bits(64, 1, True), 8)
    same("ball of radius 2 in 75 bits", ball(75, 2), 2851)
    same("ball of radius 2 in 76 bits", ball(76, 2), 2927)
    require("2851 exceeds 2^11 = 2048, so 11 check bits fail", ball(75, 2) > (1 << 11))
    require("2927 fits in 2^12 = 4096, so 12 check bits pass", ball(76, 2) <= (1 << 12))
    same("least check bits for DEC", least_check_bits(64, 2, False), 12)
    same("least check bits for DEC-TED", least_check_bits(64, 2, True), 13)
    same("least check bits for DEC of 32 data bits", least_check_bits(32, 2, False), 10)
    same("ball of radius 2 in 41 bits", ball(41, 2), 862)
    same("ball of radius 2 in 42 bits", ball(42, 2), 904)


def check_gf16() -> None:
    section("GF(16) built from x^4 + x + 1")
    exp, log, wrap = GF16
    same("x^4 + x + 1 has no factor of degree 1 or 2", irreducible(GF16_POLY), True)
    for i in range(15):
        same(f"alpha^{i} as c3 c2 c1 c0", bits(exp[i], 4), STATED_GF16[i])
    same("alpha^15", bits(wrap, 4), "0001")
    same("distinct nonzero powers of alpha", len(set(exp)), 15)
    els = range(16)

    def mul(a, b):
        return poly_mulmod(a, b, GF16_POLY, 4)

    require("table product equals shift-and-add product, 256 pairs",
            all(gmul(a, b, exp, log) == mul(a, b) for a in els for b in els))
    require("every nonzero element has an inverse",
            all(any(mul(a, b) == 1 for b in els) for a in range(1, 16)))
    require("no two nonzero elements multiply to zero",
            all(mul(a, b) != 0 for a in range(1, 16) for b in range(1, 16)))
    require("multiplication is associative on all 4,096 triples",
            all(mul(mul(a, b), c) == mul(a, mul(b, c)) for a in els for b in els for c in els))
    require("multiplication distributes over addition, 4,096 triples",
            all(mul(a, b ^ c) == mul(a, b) ^ mul(a, c) for a in els for b in els for c in els))
    same("alpha^6 times alpha^13, as a power of alpha", log[gmul(exp[6], exp[13], exp, log)], 4)
    same("alpha^6 times alpha^13, shift and add", bits(mul(exp[6], exp[13]), 4), "0011")
    same("inverse of alpha^5, as a power of alpha", log[gdiv(1, exp[5], exp, log)], 10)


def check_dec_identities() -> None:
    section("The algebra behind the DEC proof, checked in GF(16)")
    exp, log, _ = GF16

    def m(a, b):
        return gmul(a, b, exp, log)

    def cube(a):
        return m(a, m(a, a))

    els = range(16)
    require("(a + b)^3 = a^3 + a^2 b + a b^2 + b^3, all 256 pairs",
            all(cube(a ^ b) == cube(a) ^ m(m(a, a), b) ^ m(a, m(b, b)) ^ cube(b)
                for a in els for b in els))
    require("a^3 + b^3 = s(s^2 + ab) with s = a + b, all 256 pairs",
            all(cube(a) ^ cube(b) == m(a ^ b, m(a ^ b, a ^ b) ^ m(a, b)) for a in els for b in els))
    require("a + b + c = 0 gives a^3 + b^3 + c^3 = abc, all 256 pairs",
            all(cube(a) ^ cube(b) ^ cube(a ^ b) == m(m(a, b), a ^ b) for a in els for b in els))
    require("x^2 + s x + p has at most two roots, all 256 pairs (s, p)",
            all(sum(1 for x in els if m(x, x) ^ m(s, x) ^ p == 0) <= 2 for s in els for p in els))


def codeword_weights(cols: list) -> list:
    """Weight distribution of the code whose check matrix has these columns."""
    n = len(cols)
    syn = [0] * (1 << n)
    dist = [0] * (n + 1)
    for v in range(1, 1 << n):
        low = (v & -v).bit_length() - 1
        syn[v] = syn[v & (v - 1)] ^ cols[low]
    for v in range(1 << n):
        if syn[v] == 0:
            dist[popcount(v)] += 1
    return dist


def check_bch_15_7() -> None:
    section("The (15,7) double-error-correcting code")
    exp, log, _ = GF16
    cols = bch_columns(4, GF16, 15)
    same("column 0, read top to bottom", bits(cols[0], 8), "00010001")
    same("column 1, read top to bottom", bits(cols[1], 8), "00101000")
    same("column 6, read top to bottom", bits(cols[6], 8), "11001000")
    same("column 14, read top to bottom", bits(cols[14], 8), "10011111")
    same("distinct bottom halves among the 15 columns", len({c & 15 for c in cols}), 5)
    same("rank of the 8 x 15 check matrix", rank_gf2(cols), 8)
    same("data bits k = 15 - 8", 15 - rank_gf2(cols), 7)
    dist = codeword_weights(cols)
    same("codewords", sum(dist), 128)
    same("minimum nonzero weight", min(w for w in range(1, 16) if dist[w]), 5)
    same("weight distribution A0 to A15", " ".join(map(str, dist)), "1 0 0 0 0 18 30 15 15 30 18 0 0 0 0 1")
    synd = {0} | set(cols) | {cols[i] ^ cols[j] for i in range(15) for j in range(i + 1, 15)}
    same("distinct syndromes of patterns of weight at most 2", len(synd), 121)
    same("patterns of weight at most 2 in 15 bits", ball(15, 2), 121)

    section("Locating two errors at positions 3 and 10")
    s = cols[3] ^ cols[10]
    s1, s3 = s >> 4, s & 15
    same("syndrome, top half S1", bits(s1, 4), "1111")
    same("syndrome, bottom half S3", bits(s3, 4), "1011")
    same("S1 as a power of alpha", log[s1], 12)
    same("S3 as a power of alpha", log[s3], 7)
    same("S3 / S1 as a power of alpha", log[gdiv(s3, s1, exp, log)], 10)
    same("S1^2 as a power of alpha", log[gmul(s1, s1, exp, log)], 9)
    p = gdiv(s3, s1, exp, log) ^ gmul(s1, s1, exp, log)
    same("ab = S3 / S1 + S1^2, as bits", bits(p, 4), "1101")
    same("ab as a power of alpha", log[p], 13)
    roots = [i for i in range(15) if gmul(exp[i], exp[i], exp, log) ^ gmul(s1, exp[i], exp, log) ^ p == 0]
    same("positions whose alpha^i solves x^2 + S1 x + ab", " ".join(map(str, roots)), "3 10")
    # The division-free form, multiplied through by S1, that each position tests at once.
    s1sq = gmul(s1, s1, exp, log)
    c0 = s3 ^ gmul(s1sq, s1, exp, log)
    free = [i for i in range(15)
            if gmul(s1, gmul(exp[i], exp[i], exp, log), exp, log) ^ gmul(s1sq, exp[i], exp, log) ^ c0 == 0]
    same("positions whose alpha^i solves S1 x^2 + S1^2 x + S3 + S1^3", " ".join(map(str, free)), "3 10")


def dec_ted_columns() -> list:
    """The (79,64) DEC-TED code: the first 78 columns (beta, beta^3) over GF(128),
    each with a parity bit 1 appended, then the unit column of that parity bit."""
    return [(c << 1) | 1 for c in bch_columns(7, GF128, 78)] + [1]


def check_bch_79_64() -> None:
    section("A DEC-TED code for 64 data bits from GF(128)")
    exp, log, wrap = GF128
    same("x^7 + x + 1 has no factor of degree 1 to 3", irreducible(GF128_POLY), True)
    same("distinct nonzero powers of alpha in GF(128)", len(set(exp)), 127)
    same("alpha^127", wrap, 1)
    cols = bch_columns(7, GF128, 127)
    same("rank of the 14 x 127 check matrix", rank_gf2(cols), 14)
    same("data bits of the (127,k) code", 127 - rank_gf2(cols), 113)
    short = cols[:78]
    same("rank of the first 78 columns", rank_gf2(short), 14)
    same("data bits of the shortened code", 78 - rank_gf2(short), 64)
    synd = {0} | set(short) | {short[i] ^ short[j] for i in range(78) for j in range(i + 1, 78)}
    same("distinct syndromes, weight at most 2, length 78", len(synd), 3082)
    ext = [(c << 1) | 1 for c in short] + [1]
    require("every column of the extended matrix has parity bit 1", all(c & 1 for c in ext))
    synd_ext = {0} | set(ext) | {ext[i] ^ ext[j] for i in range(79) for j in range(i + 1, 79)}
    same("distinct syndromes, weight at most 2, length 79", len(synd_ext), 3161)
    same("check bits of the DEC-TED code, the rank of its matrix", rank_gf2(ext), 15)
    same("the same matrix as dec_ted_columns() builds it", ext == dec_ted_columns(), True)
    same("overhead of the (79,64) code, percent", pct(15, 64), "23.44")
    same("check bits above the sphere-packing minimum of 13", rank_gf2(dec_ted_columns()) - least_check_bits(64, 2, True), 2)
    same("overhead of the (78,64) DEC code, percent", pct(14, 64), "21.88")


# DAEC -------------------------------------------------------------------------

@lru_cache(maxsize=None)
def daec_search(r: int, n: int):
    """Order n distinct odd-weight r-bit columns so that the n - 1 sums of
    neighbors are distinct. Depth-first, candidates in increasing order."""
    cands = [v for v in range(1, 1 << r) if popcount(v) % 2 == 1]
    m = len(cands)
    idx = [-1] * n
    used = [False] * m
    sum_used = [False] * (1 << r)
    depth = 0
    while True:
        if depth == n:
            return tuple(cands[i] for i in idx)
        j = idx[depth]
        if j >= 0:
            used[j] = False
            if depth > 0:
                sum_used[cands[idx[depth - 1]] ^ cands[j]] = False
        j += 1
        while j < m:
            if not used[j] and (depth == 0 or not sum_used[cands[idx[depth - 1]] ^ cands[j]]):
                break
            j += 1
        if j == m:
            idx[depth] = -1
            depth -= 1
            if depth < 0:
                return None
            continue
        idx[depth] = j
        used[j] = True
        if depth > 0:
            sum_used[cands[idx[depth - 1]] ^ cands[j]] = True
        depth += 1


def daec_outcomes(cols: list) -> tuple:
    """Run the decoder on every single and double error.

    Returns singles corrected, adjacent doubles corrected, nonadjacent doubles
    miscorrected and nonadjacent doubles flagged."""
    n = len(cols)
    single = {c: i for i, c in enumerate(cols)}
    adjacent = {cols[i] ^ cols[i + 1]: i for i in range(n - 1)}
    ok1 = sum(1 for i in range(n) if single.get(cols[i]) == i)
    ok2 = sum(1 for i in range(n - 1) if adjacent.get(cols[i] ^ cols[i + 1]) == i)
    mis = flag = 0
    for i in range(n):
        for j in range(i + 2, n):
            s = cols[i] ^ cols[j]
            if s in single or s in adjacent:
                mis += 1
            else:
                flag += 1
    return ok1, ok2, mis, flag


def daec_lower_bound(r: int, n: int) -> int:
    """Each nonzero even syndrome is the sum of at least (2n - 2^(r-1)) / 2 pairs."""
    pairs = max(0, 2 * n - (1 << (r - 1))) // 2
    return (n - 1) * max(0, pairs - 1)


def check_daec() -> None:
    section("SEC-DED-DAEC with odd-weight columns")
    same("odd-weight columns with 4 rows", sum(1 for v in range(16) if popcount(v) % 2), 8)
    same("nonzero even-weight syndromes with 4 rows", sum(1 for v in range(1, 16) if popcount(v) % 2 == 0), 7)
    total = 0
    for v in range(1, 16):
        if popcount(v) % 2 == 0:
            total ^= v
    same("sum of all seven nonzero even 4-bit words", bits(total, 4), "0000")
    same("a search for an (8,4) code finds one", daec_search(4, 8) is not None, False)
    c7 = daec_search(4, 7)
    same("columns of the (7,3) code", " ".join(bits(c, 4) for c in c7), "0001 0010 0100 1000 0111 1110 1011")
    sums = [c7[i] ^ c7[i + 1] for i in range(6)]
    same("neighbor sums of the (7,3) code", " ".join(bits(s, 4) for s in sums), "0011 0110 1100 1111 1001 0101")
    unused = [v for v in range(1, 16) if popcount(v) % 2 == 0 and v not in sums]
    same("even syndrome no neighbor sum uses", " ".join(bits(v, 4) for v in unused), "1010")
    same("rank of the (7,3) check matrix", rank_gf2(c7), 4)
    flagged = [f"{i},{j}" for i in range(7) for j in range(i + 2, 7) if c7[i] ^ c7[j] == 10]
    same("nonadjacent doubles the (7,3) code flags", " ".join(flagged), "0,6 1,3 2,5")
    same("syndrome of a double error at 1 and 4", bits(c7[1] ^ c7[4], 4), "0101")
    same("syndrome of a double error at 2 and 3", bits(c7[2] ^ c7[3], 4), "1100")
    same("nonadjacent pairs whose syndrome is also 1100",
         " ".join(f"{i},{j}" for i in range(7) for j in range(i + 2, 7) if c7[i] ^ c7[j] == 0b1100), "1,5 4,6")
    same("syndrome of a single error at 5", bits(c7[5], 4), "1110")

    stated = {(4, 7): (6, 15, 12, 3, 12), (5, 13): (12, 66, 51, 15, 48),
              (7, 39): (38, 703, 490, 213, 228), (8, 72): (71, 2485, 1465, 1020, 497)}
    for (r, n), (adj, nonadj, mis_s, flag_s, bound_s) in stated.items():
        cols = daec_search(r, n)
        label = f"({n},{n - r})"
        require(f"{label} search finds a code", cols is not None)
        require(f"{label} every column has odd weight", all(popcount(c) % 2 for c in cols))
        same(f"{label} rank of the check matrix", rank_gf2(cols), r)
        ok1, ok2, mis, flag = daec_outcomes(cols)
        same(f"{label} single errors corrected", ok1, n)
        same(f"{label} adjacent doubles corrected", ok2, adj)
        same(f"{label} nonadjacent doubles", mis + flag, nonadj)
        same(f"{label} nonadjacent doubles miscorrected", mis, mis_s)
        same(f"{label} nonadjacent doubles flagged", flag, flag_s)
        same(f"{label} lower bound on miscorrected doubles", daec_lower_bound(r, n), bound_s)
    cols72 = daec_search(8, 72)
    same("(72,64) share of nonadjacent doubles miscorrected", pct1(1465, 2485), "59.0")
    same("(72,64) share at the lower bound", pct1(497, 2485), "20.0")
    same("(72,64) unit columns among its 72 columns", sum(1 for c in cols72 if popcount(c) == 1), 8)
    same("(72,64) columns of weight 3", sum(1 for c in cols72 if popcount(c) == 3), 47)
    same("(72,64) columns of weight 5", sum(1 for c in cols72 if popcount(c) == 5), 17)
    adj = [cols72[i] ^ cols72[i + 1] for i in range(71)]
    counts = [sum(1 for i in range(72) for j in range(i + 1, 72) if cols72[i] ^ cols72[j] == s) for s in adj]
    same("(72,64) fewest pairs sharing a neighbor sum", min(counts), 8)
    same("(72,64) neighbor sums with that fewest count, of 71", counts.count(min(counts)), 4)
    same("(72,64) neighbor sums with more pairs than that", sum(1 for c in counts if c > min(counts)), 67)
    same("pairs on each even syndrome at least (144 - 128) / 2", (2 * 72 - 128) // 2, 8)


def interleave_extra_check_bits() -> int:
    """Lay four (72,64) words out bit by bit, cell p holding bit p // 4 of word p % 4,
    and count the cells holding check bits beyond the 4 x 8 the words carry anyway."""
    layout = [(p % 4, p // 4) for p in range(4 * 72)]
    check_cells = sum(1 for _, bit in layout if bit >= 64)
    return check_cells - 4 * (72 - 64)


def check_interleave() -> None:
    section("Four-way bit interleaving of four 72-bit words")
    cells = 4 * 72
    same("cells in the interleaved row", cells, 288)
    same("word holding cell 13", 13 % 4, 1)
    same("bit of that word in cell 13", 13 // 4, 3)
    worst = []
    for length in range(1, 10):
        top = 0
        for start in range(cells - length + 1):
            hits = [0, 0, 0, 0]
            for p in range(start, start + length):
                hits[p % 4] += 1
            top = max(top, max(hits))
        worst.append(top)
    same("most errors in one word, bursts of 1 to 9 cells", " ".join(map(str, worst)), "1 1 1 1 2 2 2 2 3")
    same("longest burst SEC-DED words always correct", max(L for L in range(1, 10) if worst[L - 1] <= 1), 4)
    same("longest burst SEC-DED words always detect", max(L for L in range(1, 10) if worst[L - 1] <= 2), 8)
    same("extra check bits the interleave costs", interleave_extra_check_bits(), 0)


# Single-symbol correction -----------------------------------------------------

# A binary (73,64) code found by a greedy search: one line per symbol, each column
# a 9-bit string read from top to bottom. Lines 1 to 16 are the data symbols, lines
# 17 and 18 the two 4-bit check symbols and line 19 the 1-bit check symbol.
S4EC_R9 = [
    "111010000 011100011 111111100 101111110",
    "000100110 011000001 110010100 100001101",
    "011111100 011011001 110110101 100000110",
    "101101000 010100011 111110100 001010010",
    "010111110 001111001 111000110 001000001",
    "110010110 011010010 101101111 110100101",
    "110011100 010001001 011111011 001010110",
    "111111101 001010011 101000101 111100000",
    "110011110 001100111 010110110 110111100",
    "011111110 101100001 101111010 010100010",
    "110111011 100011101 100110100 100001001",
    "001000011 101110001 000010111 000110110",
    "110111010 100000111 010001010 101101001",
    "100011010 110110000 001101100 110001011",
    "101000110 011011110 010110001 010000011",
    "000011110 100100010 111100001 001110110",
    "000000001 000000010 000000100 000001000",
    "000010000 000100000 001000000 010000000",
    "100000000",
]


def symbol_patterns(n: int, b: int) -> int:
    """Nonzero errors inside one symbol when n bits are cut into b-bit symbols in
    order, the last one shorter when b does not divide n."""
    return (n // b) * ((1 << b) - 1) + (1 << (n % b)) - 1


def least_binary_symbol_bits(k: int, b: int) -> int:
    """Least r with 1 + symbol_patterns(k + r, b) <= 2^r."""
    r = 1
    while 1 + symbol_patterns(k + r, b) > (1 << r):
        r += 1
    return r


def binary_symbol_groups() -> list:
    return [[int(c, 2) for c in line.split()] for line in S4EC_R9]


def group_patterns(group: list) -> list:
    """Syndromes of every nonzero error inside one symbol, in order of the error value."""
    out = []
    for e in range(1, 1 << len(group)):
        v = 0
        for j, c in enumerate(group):
            if (e >> j) & 1:
                v ^= c
        out.append(v)
    return out


def double_symbol_outcomes(patterns: list) -> tuple:
    """Errors in two different symbols: how many are miscorrected, flagged or unseen."""
    single = {v for p in patterns for v in p}
    mis = flag = unseen = 0
    for i in range(len(patterns)):
        for j in range(i + 1, len(patterns)):
            for a in patterns[i]:
                for b in patterns[j]:
                    v = a ^ b
                    if v == 0:
                        unseen += 1
                    elif v in single:
                        mis += 1
                    else:
                        flag += 1
    return mis, flag, unseen


def symbol_columns() -> list:
    """Columns of the (19,16) code over GF(16): (1, g, g^2) for g = alpha^0..alpha^14,
    then (0, 1, 1), then the three unit columns for the check symbols."""
    exp, log, _ = GF16
    cols = [(1, exp[i], gmul(exp[i], exp[i], exp, log)) for i in range(15)]
    cols.append((0, 1, 1))
    cols += [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    return cols


def check_symbol_code() -> None:
    section("Single 4-bit symbol correction for 64 data bits")
    exp, log, _ = GF16
    same("nonzero error patterns inside one 4-bit symbol", sum(1 for e in range(1 << 4) if e), 15)
    same("longest code with 2 check symbols, (16^2 - 1) / 15", (16 ** 2 - 1) // 15, 17)
    same("data bits that 17 symbols leave room for", (17 - 2) * 4, 60)
    same("syndromes 18 symbols need, 1 + 15 x 18", 1 + 15 * 18, 271)
    require("271 exceeds 16^2 = 256", 1 + 15 * 18 > 256)
    same("longest code with 3 check symbols, (16^3 - 1) / 15", (16 ** 3 - 1) // 15, 273)
    same("syndromes 19 symbols need, 1 + 15 x 19", 1 + 15 * 19, 286)
    cols = symbol_columns()
    same("symbols in the code", len(cols), 19)

    def scale(lam, col):
        return tuple(gmul(lam, x, exp, log) for x in col)

    require("no column is zero", all(any(c) for c in cols))
    require("no column is a nonzero multiple of another, 171 pairs",
            all(scale(lam, cols[i]) != cols[j]
                for i in range(19) for j in range(19) if i != j for lam in range(1, 16)))

    def binary_syndrome(i, v):
        a, b, c = scale(v, cols[i])
        return (a << 8) | (b << 4) | c

    synd = {binary_syndrome(i, v) for i in range(19) for v in range(1, 16)}
    same("distinct nonzero syndromes of single-symbol errors", len(synd - {0}), 285)
    binary_cols = [binary_syndrome(i, 1 << j) for i in range(19) for j in range(4)]
    same("rank of the 12 x 76 binary check matrix", rank_gf2(binary_cols), 12)
    same("binary length and data bits", f"{len(binary_cols)},{len(binary_cols) - 12}", "76,64")
    same("overhead of the (76,64) code, percent", pct(12, 64), "18.75")
    require("the same matrix as symbol_binary_columns() builds it", binary_cols == symbol_binary_columns())

    section("A binary code for one bad 4-bit symbol with fewer check bits")
    same("8 check bits: syndromes needed, 1 + 18 x 15", 1 + symbol_patterns(72, 4), 271)
    require("271 exceeds 2^8 = 256, so 8 check bits fail", 1 + symbol_patterns(72, 4) > (1 << 8))
    same("9 check bits: syndromes needed, 1 + 18 x 15 + 1", 1 + symbol_patterns(73, 4), 272)
    require("272 fits in 2^9 = 512, so 9 check bits pass the count", 1 + symbol_patterns(73, 4) <= (1 << 9))
    same("least check bits for a binary single-symbol code", least_binary_symbol_bits(64, 4), 9)
    groups = binary_symbol_groups()
    widths = [len(g) for g in groups]
    same("symbols of 4 bits and of 1 bit in the 73-bit word", f"{widths.count(4)} {widths.count(1)}", "18 1")
    flat = [c for g in groups for c in g]
    require("every column has 9 rows", all(c < (1 << 9) for c in flat))
    same("rank of the 9 x 73 binary check matrix", rank_gf2(flat), 9)
    same("binary length and data bits", f"{len(flat)},{len(flat) - rank_gf2(flat)}", "73,64")
    require("the last 9 columns are the 9 unit columns, so the word is systematic",
            sorted(flat[64:]) == [1 << i for i in range(9)])
    pats = [group_patterns(g) for g in groups]
    every = [v for p in pats for v in p]
    same("single-symbol errors in the (73,64) code", len(every), 271)
    same("distinct nonzero syndromes they give", len(set(every) - {0}), 271)
    same("overhead of the (73,64) code, percent", pct(9, 64), "14.06")
    same("4-bit chips the 73-bit and 76-bit words fill", f"{-(-73 // 4)} {-(-76 // 4)}", "19 19")
    same("4-bit chips the 9 check bits fill", -(-rank_gf2(flat) // 4), 3)
    same("(76,64) check bits above the binary minimum", rank_gf2(symbol_binary_columns()) - least_binary_symbol_bits(64, 4), 3)

    section("Errors in two different symbols")
    mis73, flag73, unseen73 = double_symbol_outcomes(pats)
    same("(73,64): double-symbol errors", mis73 + flag73 + unseen73, 34695)
    same("(73,64): miscorrected, flagged, unseen", f"{mis73} {flag73} {unseen73}", "17553 17142 0")
    same("(73,64): share miscorrected, percent", pct1(mis73, 34695), "50.6")
    gf_pats = [[binary_syndrome(i, v) for v in range(1, 16)] for i in range(19)]
    mis76, flag76, unseen76 = double_symbol_outcomes(gf_pats)
    same("(76,64): double-symbol errors", mis76 + flag76 + unseen76, 38475)
    same("(76,64): miscorrected, flagged, unseen", f"{mis76} {flag76} {unseen76}", "405 38070 0")
    same("(76,64): share miscorrected, percent", pct1(mis76, 38475), "1.1")

    section("Correcting the error value alpha^2 in data symbol 5")
    s1, s2, s3 = scale(exp[2], cols[5])
    same("syndrome S1, S2, S3 as bits", f"{bits(s1, 4)} {bits(s2, 4)} {bits(s3, 4)}", "0100 1011 1111")
    same("S1, S2, S3 as powers of alpha", f"{log[s1]} {log[s2]} {log[s3]}", "2 7 12")
    same("S2 / S1 as a power of alpha, the locator", log[gdiv(s2, s1, exp, log)], 5)
    same("S3 / S2 as a power of alpha, the same locator", log[gdiv(s3, s2, exp, log)], 5)
    same("bits of symbol 5", f"d{4 * 5} to d{4 * 5 + 3}", "d20 to d23")

    section("Byte symbols for 64 data bits")
    same("longest code with 2 byte check symbols, (256^2 - 1) / 255", (256 ** 2 - 1) // 255, 257)
    same("syndromes 10 byte symbols need, 1 + 255 x 10", 1 + 255 * 10, 2551)
    same("overhead of the (80,64) byte code, percent", pct(16, 64), "25")


# Unidirectional errors --------------------------------------------------------

def zeros(word: int, k: int) -> int:
    return k - popcount(word)


def berger_bits(k: int) -> int:
    """ceil(log2(k + 1)), the length of k written in binary."""
    r = 0
    while (1 << r) < k + 1:
        r += 1
    return r


def covering_pairs(code: list) -> int:
    """Ordered pairs of distinct codewords X, Y where every 1 of Y is a 1 of X."""
    return sum(1 for x in code for y in code if x != y and x & y == y)


def least_covering_distance(code: list):
    best = None
    for x in code:
        for y in code:
            if x != y and x & y == y:
                d = popcount(x ^ y)
                best = d if best is None or d < best else best
    return best


def check_berger() -> None:
    section("The Berger code")
    same("check bits for 8 data bits", berger_bits(8), 4)
    same("check bits for 64 data bits", berger_bits(64), 7)
    same("check bits for 32 data bits", berger_bits(32), 6)
    same("check bits for 7 data bits", berger_bits(7), 3)
    same("overhead for 64 data bits, percent", pct(7, 64), "10.94")
    same("overhead for 32 data bits, percent", pct(6, 32), "18.75")
    d = 0b10110010
    same("zeros in 10110010", zeros(d, 8), 4)
    same("check of 10110010", bits(zeros(d, 8), 4), "0100")
    e = 0b10100010
    same("check recomputed after 10110010 becomes 10100010", bits(zeros(e, 8), 4), "0101")
    e = 0b11110010
    same("check recomputed after 10110010 becomes 11110010", bits(zeros(e, 8), 4), "0011")
    code = [(w << 4) | zeros(w, 8) for w in range(256)]
    same("codeword pairs where one covers the other, k = 8", covering_pairs(code), 0)
    ones = [(w << 4) | popcount(w) for w in range(256)]
    same("the same count when the check counts ones", covering_pairs(ones) > 0, True)
    x, y = (0b11011010 << 4) | 0b0101, (0b11010010 << 4) | 0b0100
    require("11011010 0101 and 11010010 0100 are both ones-count words", x in ones and y in ones)
    same("bits the 1-to-0 error flips", popcount(x ^ y), 2)
    same("least check bits for all-unidirectional detection, k = 64", 7, berger_bits(64))


def f_mod(r: int):
    return lambda z: z % (1 << r)


def f_bl4(z: int) -> int:
    return (z % 8) + 4


def f_bl1(r: int):
    return lambda z: (z % (1 << (r - 1))) + (1 << (r - 2))


F_MAP = [0b0011, 0b0101, 0b0110, 0b1001, 0b1010, 0b1100]


def f_bl2(r: int):
    low = r - 4

    def f(z: int) -> int:
        v = z % (6 * (1 << low))
        return (F_MAP[v >> low] << low) | (v & ((1 << low) - 1))
    return f


def detects_up_to(f, k: int):
    """Largest t such that every unidirectional error of weight at most t is
    detected, from the zero counts alone. None means every weight is detected."""
    best = None
    for z in range(k + 1):
        for j in range(1, k - z + 1):
            cz, cj = f(z), f(z + j)
            if cz & cj == cj:
                d = j + popcount(cz ^ cj)
                best = d if best is None or d < best else best
    return None if best is None else best - 1


def check_bose_lin() -> None:
    section("Counting zeros modulo 2^r")
    for r in range(2, 7):
        same(f"zeros mod 2^{r} detects up to, k = 64", detects_up_to(f_mod(r), 64), r)
    same("zeros mod 2^7, k = 64, is the Berger code: every weight", detects_up_to(f_mod(7), 64) is None, True)
    d = 0b11011010
    same("zeros in 11011010", zeros(d, 8), 3)
    same("its 2-bit check", bits(zeros(d, 8) % 4, 2), "11")
    same("zeros after 11011010 becomes 10011010", zeros(0b10011010, 8), 4)
    same("the check that word needs", bits(zeros(0b10011010, 8) % 4, 2), "00")

    section("Bose-Lin codes")
    same("(zeros mod 8) + 4 detects up to, k = 64", detects_up_to(f_bl4, 64), 6)
    same("check of a word with 2 zeros", bits(f_bl4(2), 4), "0110")
    same("check of a word with 8 zeros", bits(f_bl4(8), 4), "0100")
    stated1 = {5: 11, 6: 20, 7: 37}
    stated2 = {5: 11, 6: 22, 7: 43}
    for r in range(5, 8):
        same(f"offset count, r = {r}, detects up to", detects_up_to(f_bl1(r), 64), stated1[r])
        same(f"offset count, r = {r}, 2^(r-2) + r - 2", (1 << (r - 2)) + r - 2, stated1[r])
        same(f"second method, r = {r}, detects up to", detects_up_to(f_bl2(r), 64), stated2[r])
        same(f"second method, r = {r}, 5 x 2^(r-4) + r - 4", 5 * (1 << (r - 4)) + r - 4, stated2[r])
    same("offset 2^(r-1) instead, r = 5, detects only", detects_up_to(lambda z: (z % 16) + 16, 64), 4)
    same("Berger, k = 64, detects every weight", detects_up_to(lambda z: z, 64) is None, True)

    section("Every codeword pair, k = 8, against the zero-count rule")
    for name, f, r in (("zeros mod 4", f_mod(2), 2), ("zeros mod 8", f_mod(3), 3),
                       ("(zeros mod 8) + 4", f_bl4, 4), ("Berger", lambda z: z, 4)):
        code = [(w << r) | f(zeros(w, 8)) for w in range(256)]
        brute = least_covering_distance(code)
        brute_t = "every weight" if brute is None else str(brute - 1)
        rule = detects_up_to(f, 8)
        rule_t = "every weight" if rule is None else str(rule)
        same(f"{name}: all 65,536 pairs agree with the rule", brute_t, rule_t)


# Choosing a code --------------------------------------------------------------

def least_daec_bits(k: int) -> int:
    """Least r with 1 + n + (n - 1) <= 2^r for n = k + r: zero, singles, adjacent doubles."""
    r = 1
    while 2 * (k + r) > (1 << r):
        r += 1
    return r


def least_symbol_checks(k: int, b: int) -> int:
    """Least check symbols c with k / b + c symbols within (q^c - 1) / (q - 1), q = 2^b."""
    q = 1 << b
    c = 1
    while k // b + c > (q ** c - 1) // (q - 1):
        c += 1
    return c


def symbol_binary_columns() -> list:
    """The 76 binary columns of the (76,64) code over GF(16), four per symbol."""
    exp, log, _ = GF16
    out = []
    for col in symbol_columns():
        for j in range(4):
            a, b, c = (gmul(1 << j, x, exp, log) for x in col)
            out.append((a << 8) | (b << 4) | c)
    return out


def choice_rows() -> list:
    """Each row's check bits are the rank of a check matrix this file builds, or the
    width of a check it computes, so none of them is typed in."""
    bose_lin_width = max(f_bl4(z) for z in range(65)).bit_length()
    return [
        ("one random bit", "Hamming SEC (71,64)", rank_gf2(hamming_columns())),
        ("one random bit, two detected", "Hsiao SEC-DED (72,64)", rank_gf2(hsiao_columns())),
        ("two random bits, three detected", "DEC-TED (79,64)", rank_gf2(dec_ted_columns())),
        ("two adjacent bits", "SEC-DED-DAEC (72,64)", rank_gf2(daec_search(8, 72))),
        ("one 4-bit symbol", "S4EC (76,64)", rank_gf2(symbol_binary_columns())),
        ("unidirectional, any weight", "Berger (71,64)", berger_bits(64)),
        ("unidirectional, at most 6 bits", "Bose-Lin (68,64)", bose_lin_width),
    ]


def check_choice_table() -> None:
    section("Choosing a code for 64 data bits")
    stated = ["7 10.94", "8 12.5", "15 23.44", "8 12.5", "12 18.75", "7 10.94", "4 6.25"]
    for (model, code, r), s in zip(choice_rows(), stated):
        same(f"{model}: {code}", f"{r} {pct(r, 64)}", s)
    same("sphere-packing minimum, one random bit", least_check_bits(64, 1, False), 7)
    same("sphere-packing minimum, two random bits, three detected", least_check_bits(64, 2, True), 13)
    same("Bose-Lin with 4 check bits detects up to", detects_up_to(f_bl4, 64), 6)
    same("SEC-DED-DAEC counting minimum, 1 + n + (n - 1) <= 2^r", least_daec_bits(64), 8)
    same("counting at r = 7: 1 + 71 + 70", 1 + 71 + 70, 142)
    same("single 4-bit symbol: least check symbols", least_symbol_checks(64, 4), 3)
    same("offset count, r = 4, is the Bose-Lin (zeros mod 8) + 4", detects_up_to(f_bl1(4), 64), 6)
    same("4-way interleave: extra check bits", interleave_extra_check_bits(), 0)
    lengths = [len(hamming_columns()), len(hsiao_columns()), len(dec_ted_columns()),
               len(daec_search(8, 72)), len(symbol_binary_columns())]
    same("lengths of the five built matrices", " ".join(map(str, lengths)), "71 72 79 72 76")
    same("data bits each leaves, length minus rank",
         " ".join(str(n - r) for n, (_, _, r) in zip(lengths, choice_rows())), "64 64 64 64 64")
    same("single 4-bit symbol, binary counting minimum", least_binary_symbol_bits(64, 4), 9)
    same("single 4-bit symbol, whole GF(16) check symbols", 4 * least_symbol_checks(64, 4), 12)


def main() -> None:
    check_opening()
    check_sphere_counts()
    check_baselines()
    check_gf16()
    check_dec_identities()
    check_bch_15_7()
    check_bch_79_64()
    check_daec()
    check_interleave()
    check_symbol_code()
    check_berger()
    check_bose_lin()
    check_choice_table()
    print("\nAll chapter 21 numbers reproduced.")


if __name__ == "__main__":
    main()
