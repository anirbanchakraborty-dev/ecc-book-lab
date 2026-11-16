"""Reproduce every number in chapter 25 of Error Correction Coding, "Golay Codes".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch25_golay_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in the chapter is a whole number, a binary or ternary word or a
matrix, so every comparison is exact. Both codes are checked exhaustively: all
4,096 codewords of the extended binary Golay code and all 729 codewords of the
extended ternary Golay code are enumerated, every 5-subset of the 24
coordinates is shown to lie in exactly one octad, every weight-4 word of
length 23 is shown to lie within distance 3 of exactly one codeword of the
perfect code, and the decoder is run on all 2,324 error patterns of weight 1
to 3 and all 10,626 of weight 4.

A binary word of length 24 is held as an integer whose most significant bit
is coordinate 0, so format(x, "024b") prints it left to right. The left half,
coordinates 0 to 11, is the top 12 bits, and the right half, coordinates 12
to 23, is the bottom 12 bits. Row i of B is held the same way as a 12-bit
integer. The module also holds the functions that ch25_charts.py draws, so
the figures and this check share one source.
"""

import sys
from itertools import combinations, product
from math import comb

SQUARES0 = (0, 1, 3, 4, 5, 9)   # 0 and the nonzero squares modulo 11
CHI5 = (0, 1, 2, 2, 1)          # the quadratic character modulo 5, with -1 written as 2


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


def bits(x: int, n: int) -> str:
    return format(x, f"0{n}b")


def pop(x: int) -> int:
    return bin(x).count("1")


def support(x: int, n: int = 24) -> list:
    return [k for k in range(n) if (x >> (n - 1 - k)) & 1]


def mask(positions, n: int = 24) -> int:
    m = 0
    for p in positions:
        m |= 1 << (n - 1 - p)
    return m


def distribution(weights) -> str:
    counts = {}
    for w in weights:
        counts[w] = counts.get(w, 0) + 1
    return " ".join(f"{w}:{counts[w]}" for w in sorted(counts))


# ---------------------------------------------------------------- the binary code

def b_matrix() -> list:
    """B as a 12-by-12 list of bits: squares pattern, then the border."""
    rows = []
    for i in range(11):
        rows.append([1 if (i + j) % 11 in SQUARES0 else 0 for j in range(11)] + [1])
    rows.append([1] * 11 + [0])
    return rows


B = b_matrix()
BROW = [int("".join(str(x) for x in row), 2) for row in B]
UNIT = [1 << (11 - i) for i in range(12)]


def times_b(v: int) -> int:
    """The 12-bit row vector v times B over GF(2)."""
    out = 0
    for k in range(12):
        if (v >> (11 - k)) & 1:
            out ^= BROW[k]
    return out


def encode(a: int) -> int:
    return (a << 12) | times_b(a)


CODE24 = [encode(a) for a in range(4096)]
GROWS = [encode(UNIT[i]) for i in range(12)]


def puncture_last(c: int) -> int:
    return c >> 1


CODE23 = [puncture_last(c) for c in CODE24]


def overlap11(i: int, j: int) -> int:
    return sum(B[i][k] & B[j][k] for k in range(11))


def octads() -> list:
    return [c for c in CODE24 if pop(c) == 8]


# ---------------------------------------------------------------- the decoder

def syndrome(r: int) -> int:
    return (r >> 12) ^ times_b(r & 0xFFF)


def decode(r: int):
    """Return (case, error) for a received word r; case 5 means flagged."""
    s = syndrome(r)
    if pop(s) <= 3:
        return 1, s << 12
    for i in range(12):
        if pop(s ^ BROW[i]) <= 2:
            return 2, ((s ^ BROW[i]) << 12) | UNIT[i]
    t = times_b(s)
    if pop(t) <= 3:
        return 3, t
    for i in range(12):
        if pop(t ^ BROW[i]) <= 2:
            return 4, (UNIT[i] << 12) | (t ^ BROW[i])
    return 5, None


def decoder_cases() -> dict:
    """For each error weight 1 to 4, how many patterns each case resolves."""
    table = {}
    k = 0
    for w in range(1, 5):
        counts = [0, 0, 0, 0, 0]
        for pos in combinations(range(24), w):
            e = mask(pos)
            c = CODE24[k % 4096]
            k += 1
            case, found = decode(c ^ e)
            if case < 5 and found != e:
                counts[4] += 1000000
            counts[case - 1] += 1
        table[w] = counts
    return table


# ---------------------------------------------------------------- the ternary code

def b3_matrix() -> list:
    rows = [[0, 1, 1, 1, 1, 1]]
    for i in range(1, 6):
        rows.append([1] + [CHI5[(j - i) % 5] for j in range(1, 6)])
    return rows


B3 = b3_matrix()


def encode3(a) -> list:
    right = [sum(a[k] * B3[k][j] for k in range(6)) % 3 for j in range(6)]
    return list(a) + right


CODE12 = [encode3(a) for a in product(range(3), repeat=6)]
CODE11 = [c[:11] for c in CODE12]


def wt3(c) -> int:
    return sum(1 for x in c if x)


# ---------------------------------------------------------------- the hardware

def chain_gates(n: int, m: int) -> int:
    """Two-input AND and OR gates in a thermometer chain that flags weight >= m."""
    gates = 0
    for k in range(1, n):
        gates += 1
        for j in range(2, min(k + 1, m) + 1):
            gates += 2 if j <= k else 1
    return gates


def chain_flag(x: int, n: int, m: int) -> int:
    t = [0] * (m + 1)
    for k in range(n):
        bit = (x >> (n - 1 - k)) & 1
        for j in range(m, 1, -1):
            t[j] = t[j] | (t[j - 1] & bit)
        t[1] = t[1] | bit
    return t[m]


def column_weights() -> list:
    return [sum(B[k][j] for k in range(12)) for j in range(12)]


def gate_counts() -> dict:
    cols = column_weights()
    syn = sum(cols)
    deep = max((1 + c - 1).bit_length() for c in cols)
    sb = sum(c - 1 for c in cols)
    sb_deep = max((c - 1).bit_length() for c in cols)
    inv = sum(pop(r) for r in BROW)
    return {
        "syndrome_xor": syn, "syndrome_depth": deep,
        "sb_xor": sb, "sb_depth": sb_deep,
        "fold_inverters": inv, "literal_xor": 12 * 12,
        "chain3": chain_gates(12, 3), "chain4": chain_gates(12, 4),
        "shared_inverters": pop(or_rows()), "flag_inversions": 2 * (len(BROW) + 1),
        "series_depth": deep + sb_deep, "correction_xor": 2 * len(BROW),
    }


def or_rows() -> int:
    """The bits of s that some row of B complements, one shared inverter each."""
    out = 0
    for r in BROW:
        out |= r
    return out


# ---------------------------------------------------------------- the checks

STATED_B = ["110111000101", "101110001011", "011100010111", "111000101101",
            "110001011011", "100010110111", "000101101111", "001011011101",
            "010110111001", "101101110001", "011011100011", "111111111110"]
STATED_MEET = {1: "0 3 4", 2: "1 3 9", 3: "0 1 9", 4: "0 1 5", 5: "0 4 9"}


def check_b() -> None:
    section("The matrix B")
    for i in range(12):
        same(f"row {i} of B", bits(BROW[i], 12), STATED_B[i])
    require("B is symmetric", all(B[i][j] == B[j][i] for i in range(12) for j in range(12)))
    require("each of rows 1 to 10 is the row above shifted left",
            all(B[i][:11] == B[i - 1][1:11] + B[i - 1][:1] for i in range(1, 11)))
    same("weights of rows 0 to 10", join(sorted({pop(BROW[i]) for i in range(11)})), "7")
    same("weight of row 11", pop(BROW[11]), 11)
    same("column weights of B", join(column_weights()), "7 7 7 7 7 7 7 7 7 7 7 11")

    section("Row overlaps by shift")
    S = set(SQUARES0)
    for t in range(1, 6):
        meet = sorted(S & {(x - t) % 11 for x in S})
        same(f"S intersected with S - {t}", join(meet), STATED_MEET[t])
        same(f"overlap of rows 0 and {t} in columns 0 to 10", overlap11(0, t), 3)
        same(f"overlap of rows 0 and {11 - t} in columns 0 to 10", overlap11(0, 11 - t), 3)
    pairs = [(i, j) for i in range(11) for j in range(i + 1, 11)]
    same("pairs among rows 0 to 10", len(pairs), 55)
    same("pairs whose shift is t or -t, for t = 1 to 5",
         join(sum(1 for i, j in pairs if (j - i) % 11 in (t, 11 - t)) for t in range(1, 6)),
         "11 11 11 11 11")
    require("overlap of any two of rows 0 to 10 is 3", all(overlap11(i, j) == 3 for i, j in pairs))
    same("overlap of row 11 with each of rows 0 to 10",
         join(sorted({pop(BROW[11] & BROW[i]) for i in range(11)})), "6")
    prod = [[sum(B[i][k] * B[k][j] for k in range(12)) % 2 for j in range(12)] for i in range(12)]
    require("B times B is the identity over GF(2)",
            all(prod[i][j] == (1 if i == j else 0) for i in range(12) for j in range(12)))
    dists = sorted({pop(BROW[i] ^ BROW[j]) for i in range(12) for j in range(i + 1, 12)})
    same("distances between two different rows of B", join(dists), "6")


def check_g24() -> None:
    section("The extended binary Golay code G24")
    require("G times G transpose is zero",
            all(pop(GROWS[i] & GROWS[j]) % 2 == 0 for i in range(12) for j in range(12)))
    same("weights of the rows of G", join(sorted({pop(g) for g in GROWS})), "8 12")
    same("codewords enumerated", len(set(CODE24)), 4096)
    require("every codeword weight is divisible by 4", all(pop(c) % 4 == 0 for c in CODE24))
    same("codewords of weight 4", sum(1 for c in CODE24 if pop(c) == 4), 0)
    same("minimum distance", min(pop(c) for c in CODE24 if c), 8)
    total = 0
    for g in GROWS:
        total ^= g
    same("sum of the 12 rows of G", bits(total, 24), "1" * 24)
    require("the all-ones word is a codeword", (1 << 24) - 1 in set(CODE24))
    same("weight distribution of G24", distribution(pop(c) for c in CODE24),
         "0:1 8:759 12:2576 16:759 24:1")
    same("A_12 = 4096 - 2 - 2 * 759", 4096 - 2 - 2 * 759, 2576)
    require("r_L + r_R B equals r_L + r_R B transpose for every r_R",
            all(times_b(v) == sum_cols(v) for v in range(4096)))
    require("every codeword has syndrome zero", all(syndrome(c) == 0 for c in CODE24))


def sum_cols(v: int) -> int:
    """v times B transpose, built from the columns of B."""
    out = 0
    for j in range(12):
        bit = 0
        for k in range(12):
            bit ^= ((v >> (11 - k)) & 1) & B[j][k]
        out |= bit << (11 - j)
    return out


def check_octads() -> None:
    section("Octads and the Steiner system S(5,8,24)")
    octs = octads()
    same("C(24,5)", comb(24, 5), 42504)
    same("C(8,5)", comb(8, 5), 56)
    same("A_8 = C(24,5) / C(8,5)", comb(24, 5) // comb(8, 5), 759)
    same("octads found", len(octs), 759)
    for size, stated in ((5, 1), (4, 5), (3, 21), (2, 77), (1, 253)):
        counts = {}
        for o in octs:
            for sub in combinations(support(o), size):
                counts[sub] = counts.get(sub, 0) + 1
        same(f"{size}-subsets reached by the octads", len(counts), comb(24, size))
        same(f"octads through each {size}-subset", join(sorted(set(counts.values()))), str(stated))
    target = mask([0, 1, 2, 3, 4])
    found = [support(o) for o in octs if o & target == target]
    same("octads containing 0 1 2 3 4", len(found), 1)
    same("the octad containing 0 1 2 3 4", join(found[0]), "0 1 2 3 4 6 14 17")
    partial, x = [], 0
    for i in range(5):
        x ^= BROW[i]
        if i:
            partial.append(bits(x, 12))
    same("partial sums b_0 + b_1, ..., b_0 + ... + b_4", " ".join(partial),
         "011001001110 000101011001 111101110100 001100101111")
    same("weight of x = b_0 + ... + b_4", pop(x), 7)
    same("row 6 of B", bits(BROW[6], 12), "000101101111")
    same("x + b_6", bits(x ^ BROW[6], 12), "001001000000")
    octad = encode(int("111110100000", 2))
    same("octad of the codeword with message 111110100000", join(support(octad)), "0 1 2 3 4 6 14 17")


def check_g23() -> None:
    section("The perfect binary Golay code G23")
    same("distinct words after puncturing", len(set(CODE23)), 4096)
    same("minimum distance of G23", min(pop(c) for c in CODE23 if c), 7)
    v = sum(comb(23, i) for i in range(4))
    same("V(23,3) = 1 + 23 + 253 + 1771", v, 2048)
    same("2^12 * V(23,3)", 4096 * v, 2 ** 23)
    same("weight distribution of G23", distribution(pop(c) for c in CODE23),
         "0:1 7:253 8:506 11:1288 12:1288 15:506 16:253 23:1")
    same("C(23,4)", comb(23, 4), 8855)
    same("C(7,4)", comb(7, 4), 35)
    same("A_7 = C(23,4) / C(7,4)", comb(23, 4) // comb(7, 4), 253)
    octs = octads()
    same("octads with a 1 in coordinate 23", sum(1 for o in octs if o & 1), 253)
    same("octads with a 0 in coordinate 23", sum(1 for o in octs if not o & 1), 506)
    near = [c for c in CODE23 if pop(c) <= 7]
    good = 0
    for pos in combinations(range(23), 4):
        w = mask(pos, 23)
        close = [c for c in near if pop(c ^ w) <= 3]
        if len(close) == 1 and pop(close[0]) == 7 and close[0] & w == w:
            good += 1
    same("weight-4 words with one codeword within 3, of weight 7", good, 8855)
    rebuilt = [(c << 1) | (pop(c) % 2) for c in CODE23]
    require("extending G23 gives back G24", sorted(rebuilt) == sorted(CODE24))


def check_ternary() -> None:
    section("The ternary Golay codes")
    stated = ["011111", "101221", "110122", "121012", "122101", "112210"]
    for i in range(6):
        same(f"row {i} of B3", "".join(str(x) for x in B3[i]), stated[i])
    require("B3 is symmetric", all(B3[i][j] == B3[j][i] for i in range(6) for j in range(6)))
    gram = [[sum(B3[i][k] * B3[j][k] for k in range(6)) % 3 for j in range(6)] for i in range(6)]
    require("B3 times B3 transpose is -I over GF(3)",
            all(gram[i][j] == (2 if i == j else 0) for i in range(6) for j in range(6)))
    g = [encode3([1 if k == i else 0 for k in range(6)]) for i in range(6)]
    require("G times G transpose is zero over GF(3)",
            all(sum(g[i][k] * g[j][k] for k in range(12)) % 3 == 0 for i in range(6) for j in range(6)))
    same("weights of the rows of G", join(sorted({wt3(r) for r in g})), "6")
    same("row 1 of B3 with itself, before reducing mod 3", sum(x * x for x in B3[1]), 11)
    same("row 1 of B3 with itself, mod 3", sum(x * x for x in B3[1]) % 3, 2)
    same("row 1 with row 2 of B3, before reducing mod 3", sum(x * y for x, y in zip(B3[1], B3[2])), 9)
    same("row 1 with row 2 of B3, mod 3", sum(x * y for x, y in zip(B3[1], B3[2])) % 3, 0)
    chi = [0, 1, -1, -1, 1]
    for t, stated in ((1, "0 -1 1 -1 0"), (2, "0 -1 -1 0 1")):
        terms = [chi[(x + t) % 5] * chi[x] for x in range(5)]
        same(f"terms chi(x + {t}) chi(x), x = 0 to 4", join(terms), stated)
        same(f"their sum for t = {t}", sum(terms), -1)
    same("row 0 of G3", "".join(str(x) for x in g[0]), "100000011111")
    same("row 1 of G3", "".join(str(x) for x in g[1]), "010000101221")
    same("digit sums of the rows of G3, mod 3", join(sum(r) % 3 for r in g), "0 2 2 2 2 2")
    same("digit sum of row 1 of G3", sum(g[1]), 8)
    same("check digit -(c0 + c7 + c8 + c9 + c10) for row 1", -(g[1][0] + sum(g[1][7:11])) % 3, 1)
    require("c11 = -(c0 + c7 + c8 + c9 + c10) on all 729 codewords",
            all(c[11] == -(c[0] + sum(c[7:11])) % 3 for c in CODE12))
    zero_sum = {tuple(c[:11] + [-sum(c[:11]) % 3]) for c in CODE12}
    require("appending the zero-sum digit to G11 does not give G12", zero_sum != {tuple(c) for c in CODE12})
    same("minimum distance of G11 with the zero-sum digit", min(wt3(c) for c in zero_sum if any(c)), 5)
    same("codewords enumerated", len({tuple(c) for c in CODE12}), 729)
    require("every codeword weight is divisible by 3", all(wt3(c) % 3 == 0 for c in CODE12))
    same("codewords of weight 3", sum(1 for c in CODE12 if wt3(c) == 3), 0)
    same("minimum distance of G12", min(wt3(c) for c in CODE12 if any(c)), 6)
    same("weight distribution of G12", distribution(wt3(c) for c in CODE12), "0:1 6:264 9:440 12:24")
    zeros = sorted({sum(1 for c in CODE12 if c[i] == 0) for i in range(12)})
    same("codewords that vanish at a given coordinate", join(zeros), "243")
    six = sorted({sum(1 for c in CODE12 if wt3(c) == 6 and c[i]) for i in range(12)})
    same("weight-6 codewords nonzero at a given coordinate", join(six), "132")
    same("C(11,3) * 2^3", comb(11, 3) * 8, 1320)
    same("C(11,3) * 2^3 / C(5,3)", comb(11, 3) * 8 // comb(5, 3), 132)
    same("A_6 = 12 * 132 / 6", 12 * 132 // 6, 264)
    same("zero entries at a coordinate: 243 - 1 - 132", 243 - 1 - 132, 110)
    same("A_9 = 12 * 110 / 3", 12 * 110 // 3, 440)
    same("A_12 = 729 - 1 - 264 - 440", 729 - 1 - 264 - 440, 24)
    v3 = 1 + 11 * 2 + comb(11, 2) * 4
    same("V_3(11,2) = 1 + 22 + 220", v3, 243)
    same("3^6 * V_3(11,2)", 729 * v3, 3 ** 11)
    same("distinct words after puncturing", len({tuple(c) for c in CODE11}), 729)
    same("minimum distance of G11", min(wt3(c) for c in CODE11 if any(c)), 5)
    same("weight distribution of G11", distribution(wt3(c) for c in CODE11),
         "0:1 5:132 6:132 8:330 9:110 11:24")
    same("weight-9 codewords of G12 with a 0 in coordinate 11",
         sum(1 for c in CODE12 if wt3(c) == 9 and c[11] == 0), 110)
    same("A_8 of G11 = 440 - 110", 440 - 110, 330)


def check_decoder() -> None:
    section("Decoding G24")
    same("patterns of weight 1, 2, 3", join(comb(24, w) for w in (1, 2, 3)), "24 276 2024")
    same("patterns of weight 1 to 3 in all", sum(comb(24, w) for w in (1, 2, 3)), 2324)
    same("patterns of weight 4", comb(24, 4), 10626)
    table = decoder_cases()
    same("weight 1 resolved by cases 1 to 4, then flagged", join(table[1]), "12 12 0 0 0")
    same("weight 2 resolved by cases 1 to 4, then flagged", join(table[2]), "66 144 66 0 0")
    same("weight 3 resolved by cases 1 to 4, then flagged", join(table[3]), "220 792 220 792 0")
    same("weight 4 resolved by cases 1 to 4, then flagged", join(table[4]), "0 0 0 0 10626")
    same("patterns of weight 1 to 3 corrected", sum(sum(table[w][:4]) for w in (1, 2, 3)), 2324)
    same("patterns of weight 1 to 3 by case",
         join(sum(table[w][k] for w in (1, 2, 3)) for k in range(4)), "298 948 286 792")


EXAMPLES = [
    ("case 1", "101000000000", (3, 7, 10)),
    ("case 2", "010010000001", (1, 6, 16)),
    ("case 3", "111100000000", (12, 17, 22)),
    ("case 4", "000000110001", (3, 14, 20)),
    ("weight 4", "100000000000", (0, 5, 13, 18)),
]
STATED_EXAMPLES = {
    "case 1": ["101000000000101011010010", "101100010010", "101011010010", "101000000000",
               "000100010010", 3, "8 6 4 10 6 6 6 8 6 6 8 8", "101000010011", 5,
               "8 4 4 6 4 4 8 6 8 4 6 8", 1, "3 7 10"],
    "case 2": ["010010000001100000101110", "000010100001", "100010101110", "100011011010",
               "100001111011", 7, "8 6 8 6 2 4 4 6 6 4 6 6", "101001100100", 5,
               "6 8 8 4 8 6 6 6 10 4 6 6", 2, "1 6 16"],
    "case 3": ["111100000000111101110100", "111100000000", "011100110110", "110010010001",
               "001110010001", 5, "6 4 4 8 8 6 8 4 4 4 6 8", "100001000010", 3,
               "6 6 8 8 4 6 6 8 10 6 6 8", 3, "12 17 22"],
    "case 4": ["000000110001110001001100", "000100110001", "111001000100", "001010011111",
               "001110101110", 7, "8 4 6 6 10 6 4 6 6 8 6 4", "110000100101", 5,
               "4 8 6 2 6 4 6 8 6 6 6 8", 4, "3 14 20"],
    "weight 4": ["100000000000110111000101", "000001000000", "100111100101", "001011100100",
                 "001010100100", 4, "7 7 7 5 11 5 7 5 7 7 5 7", "000101010010", 4,
                 "7 7 5 11 5 7 5 7 7 5 7 7", 5, "none"],
}


def check_examples() -> None:
    section("Worked decoding examples")
    for name, message, errors in EXAMPLES:
        want = STATED_EXAMPLES[name]
        a = int(message, 2)
        c = encode(a)
        r = c ^ mask(errors)
        s = syndrome(r)
        t = times_b(s)
        case, found = decode(r)
        same(f"{name}: codeword for message {message}", bits(c, 24), want[0])
        same(f"{name}: received r_L", bits(r >> 12, 12), want[1])
        same(f"{name}: received r_R", bits(r & 0xFFF, 12), want[2])
        same(f"{name}: r_R B", bits(times_b(r & 0xFFF), 12), want[3])
        same(f"{name}: syndrome s", bits(s, 12), want[4])
        same(f"{name}: weight of s", pop(s), want[5])
        same(f"{name}: weights of s + b_i, i = 0 to 11", join(pop(s ^ BROW[i]) for i in range(12)), want[6])
        same(f"{name}: s B", bits(t, 12), want[7])
        same(f"{name}: weight of s B", pop(t), want[8])
        same(f"{name}: weights of s B + b_i, i = 0 to 11", join(pop(t ^ BROW[i]) for i in range(12)), want[9])
        same(f"{name}: case that fires (5 = flagged)", case, want[10])
        same(f"{name}: error positions found", join(support(found)) if found is not None else "none", want[11])
        if found is not None:
            same(f"{name}: message recovered", bits((r ^ found) >> 12, 12), message)


def check_exercises() -> None:
    section("Exercises")
    same("row 3 of B", bits(BROW[3], 12), "111000101101")
    same("overlap of rows 0 and 3, all 12 columns", pop(BROW[0] & BROW[3]), 4)
    g01 = GROWS[0] ^ GROWS[1]
    same("weight of g_0 + g_1", pop(g01), 8)
    same("support of g_0 + g_1", join(support(g01)), "0 1 13 14 17 20 21 22")
    same("octads through 1, 2, 3, 4 given points",
         join(comb(24 - i, 5 - i) // comb(8 - i, 5 - i) for i in range(1, 5)), "253 77 21 5")
    r = mask([5, 23])
    s = syndrome(r)
    same("syndrome of an error at coordinates 5 and 23", bits(s, 12), "111110111110")
    case, found = decode(r)
    same("case and positions found for that error", f"{case} {join(support(found))}", "2 5 23")
    c3 = encode3([1, 1, 0, 0, 0, 0])
    same("ternary codeword for message 110000", "".join(str(x) for x in c3), "110000112002")
    same("its weight", wt3(c3), 6)
    stats = gate_counts()
    same("XOR gates in the syndrome trees", stats["syndrome_xor"], 88)
    same("XOR gates in the s B network", stats["sb_xor"], 76)
    same("levels from r to s B, the two networks in series", stats["series_depth"], 8)


def check_hardware() -> None:
    section("Hardware gate counts")
    stats = gate_counts()
    same("two-input XOR gates for the syndrome", stats["syndrome_xor"], 88)
    same("levels in the deepest syndrome tree", stats["syndrome_depth"], 4)
    same("two-input XOR gates for s B", stats["sb_xor"], 76)
    same("levels in the deepest s B tree", stats["sb_depth"], 4)
    same("XOR gates in one bank of 12 tests, built literally", stats["literal_xor"], 144)
    same("inverters in one bank once the constants are folded", stats["fold_inverters"], 88)
    same("inverters in one bank with the inverters shared", stats["shared_inverters"], 12)
    same("AND and OR gates flagging weight at least 3 of 12", stats["chain3"], 51)
    same("AND and OR gates flagging weight at least 4 of 12", stats["chain4"], 68)
    require("both chains agree with a popcount on all 4096 inputs",
            all(chain_flag(x, 12, 3) == (pop(x) >= 3) and chain_flag(x, 12, 4) == (pop(x) >= 4)
                for x in range(4096)))
    same("AND and OR gates in the 24 weight-at-most-2 tests", 24 * stats["chain3"], 1224)
    same("AND and OR gates in the 2 weight-at-most-3 tests", 2 * stats["chain4"], 136)
    same("AND and OR gates in one bank of 13 tests", 12 * stats["chain3"] + stats["chain4"], 680)
    same("two-input XOR gates in both networks", stats["syndrome_xor"] + stats["sb_xor"], 164)
    same("flag inversions in both banks, left to the select", stats["flag_inversions"], 26)
    same("two-input XOR gates forming r + e", stats["correction_xor"], 24)


def main() -> None:
    print("Chapter 25, Golay Codes")
    check_b()
    check_g24()
    check_octads()
    check_g23()
    check_ternary()
    check_decoder()
    check_examples()
    check_exercises()
    check_hardware()
    print("\nAll chapter 25 numbers reproduced.")


if __name__ == "__main__":
    main()
