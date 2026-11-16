"""Reproduce every number in chapter 26 of Error Correction Coding, "Reed-Muller Codes".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch26_reed_muller_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in the chapter is a whole number, a binary word, a list of
integers or a probability printed to six places, so every comparison is exact.

The points of GF(2)^m are numbered j = 0, 1, ..., 2^m - 1, and variable x_i is
bit i - 1 of j, so x_1 is the least significant bit and x_m the most
significant. Position j of a codeword is the value of the polynomial at point
j, and a word is printed with position 0 first. Every code with at most 2,048
codewords is enumerated whole, and so is RM(2,5) with 65,536. RM(3,5) and
RM(4,5) get their weights from their duals by the MacWilliams identity, and the
distances of RM(2,7), RM(3,6) and RM(1,10) come from the proved formula. Reed's majority-logic decoder is run on every
error pattern up to its guaranteed radius, on every codeword of RM(1,3),
RM(1,4) and RM(2,4), and on the all-zero codeword of RM(1,5), where all
4,514,873 patterns of weight 0 to 7 are decoded. The decoder's verdict depends
only on the error pattern, which the chapter proves, so the all-zero codeword
stands for all 64. The module also holds the functions that ch26_charts.py
draws, so the chart and this check share one source.
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


def word(v) -> str:
    return "".join(str(x) for x in v)


def vec(v) -> str:
    return " ".join(str(x) for x in v)


# ---------------------------------------------------------------- the codes

def monomials(r: int, m: int) -> list:
    """Every monomial of degree at most r in x_1..x_m, by degree, then in dictionary order."""
    return [S for d in range(r + 1) for S in combinations(range(1, m + 1), d)]


def mono_name(S) -> str:
    return "1" if not S else "".join(f"x{i}" for i in S)


def truth_table(S, m: int) -> list:
    return [1 if all((j >> (i - 1)) & 1 for i in S) else 0 for j in range(1 << m)]


def generator(r: int, m: int) -> list:
    return [truth_table(S, m) for S in monomials(r, m)]


def add(a, b) -> list:
    return [x ^ y for x, y in zip(a, b)]


def encode(coeffs, r: int, m: int) -> list:
    c = [0] * (1 << m)
    for S, a in zip(monomials(r, m), coeffs):
        if a:
            c = add(c, truth_table(S, m))
    return c


def to_int(v) -> int:
    return sum(b << j for j, b in enumerate(v))


def from_int(x: int, n: int) -> list:
    return [(x >> j) & 1 for j in range(n)]


def span_ints(rows) -> list:
    """Every codeword spanned by the rows, as integers, in Gray-code order."""
    g = [to_int(row) for row in rows]
    words, c = [0], 0
    for t in range(1, 1 << len(g)):
        c ^= g[(t & -t).bit_length() - 1]
        words.append(c)
    return words


def rm_code(r: int, m: int) -> list:
    return span_ints(generator(r, m))


def k_rm(r: int, m: int) -> int:
    return sum(comb(m, i) for i in range(r + 1))


def weight_distribution(words) -> dict:
    d = {}
    for c in words:
        w = c.bit_count()
        d[w] = d.get(w, 0) + 1
    return dict(sorted(d.items()))


def dist_text(d: dict) -> str:
    return " ".join(f"{w}:{n}" for w, n in d.items())


def min_weight(words) -> int:
    return min(c.bit_count() for c in words if c)


def mobius_stages(f, m: int) -> list:
    """The Moebius transform, one variable per stage: a[j] ^= a[j - 2^(i-1)] where bit i-1 of j is 1."""
    a, stages = list(f), []
    for i in range(m):
        a = [a[j] ^ a[j ^ (1 << i)] if (j >> i) & 1 else a[j] for j in range(len(a))]
        stages.append(a)
    return stages


def anf_terms(coeffs, m: int) -> str:
    names = [mono_name(tuple(i for i in range(1, m + 1) if (s >> (i - 1)) & 1)) for s in range(1 << m)]
    order = sorted((s for s in range(1 << m) if coeffs[s]), key=lambda s: (s.bit_count(), [i for i in range(m) if (s >> i) & 1]))
    return " + ".join(names[s] for s in order) if order else "0"


# ------------------------------------------------- Reed's majority-logic decoder

def subcube(S, b: int, m: int) -> list:
    """Positions whose variables outside S equal the bits of b, the variables in S running free."""
    fixed = [i for i in range(1, m + 1) if i not in S]
    return [j for j in range(1 << m)
            if all(((j >> (i - 1)) & 1) == ((b >> t) & 1) for t, i in enumerate(fixed))]


def reed_decode(y, r: int, m: int):
    """Returns (coefficients by monomial, or None on a tied vote, and a trace of every vote)."""
    y, coeffs, trace = list(y), {}, []
    for deg in range(r, -1, -1):
        part = [0] * (1 << m)
        for S in combinations(range(1, m + 1), deg):
            sums = [sum(y[p] for p in subcube(S, b, m)) % 2 for b in range(1 << (m - deg))]
            ones = sum(sums)
            trace.append((S, word(sums)))
            if 2 * ones == len(sums):
                return None, trace
            coeffs[S] = 1 if 2 * ones > len(sums) else 0
            if coeffs[S]:
                part = add(part, truth_table(S, m))
        y = add(y, part)
        trace.append(("after", word(y)))
    return coeffs, trace


def reed_plan(r: int, m: int) -> list:
    """For each degree from r down to 0, the folds, mask, monomial and vote size of every coefficient."""
    n = 1 << m
    zero = [sum(1 << j for j in range(n) if not (j >> (i - 1)) & 1) for i in range(1, m + 1)]
    plan = []
    for deg in range(r, -1, -1):
        level = []
        for S in combinations(range(1, m + 1), deg):
            mask = (1 << n) - 1
            for i in S:
                mask &= zero[i - 1]
            level.append((tuple(1 << (i - 1) for i in S), mask, to_int(truth_table(S, m)), 1 << (m - deg)))
        plan.append(level)
    return plan


def reed_decode_int(y: int, plan) -> int:
    """The same decoder on a word held as an integer. Returns the decoded codeword, or -1 on a tie."""
    decoded = 0
    for level in plan:
        part = 0
        for shifts, mask, mono, nsums in level:
            z = y
            for s in shifts:
                z ^= z >> s
            ones = (z & mask).bit_count()
            if 2 * ones == nsums:
                return -1
            if 2 * ones > nsums:
                part ^= mono
        y ^= part
        decoded ^= part
    return decoded


def patterns(n: int, w: int):
    """Every n-bit integer of weight w, in increasing order (Gosper's step)."""
    if w == 0:
        yield 0
        return
    v, limit = (1 << w) - 1, 1 << n
    while v < limit:
        yield v
        c = v & -v
        r = v + c
        v = (((r ^ v) >> 2) // c) | r


def exhaustive_reed(r: int, m: int, every_codeword: bool):
    """(patterns tried, decodings tried, failures) for every error pattern up to the decoder's radius."""
    n, t = 1 << m, (1 << (m - r - 1)) - 1
    plan = reed_plan(r, m)
    words = rm_code(r, m) if every_codeword else [0]
    tried = failures = 0
    npat = sum(comb(n, w) for w in range(t + 1))
    for w in range(t + 1):
        for e in patterns(n, w):
            for c in words:
                tried += 1
                if reed_decode_int(c ^ e, plan) != c:
                    failures += 1
    return npat, tried, failures


# ------------------------------------------------- the fast Hadamard transform

def fht_stages(Y) -> list:
    """Butterflies (a + b, a - b) on positions j and j + h, for h = 2^(m-1) down to 1."""
    y, n, stages = list(Y), len(Y), []
    h = n // 2
    while h >= 1:
        z = list(y)
        for j in range(n):
            if not j & h:
                z[j], z[j + h] = y[j] + y[j + h], y[j] - y[j + h]
        y = z
        stages.append(y)
        h //= 2
    return stages


def signs(y) -> list:
    return [1 - 2 * b for b in y]


def direct_correlation(Y, m: int) -> list:
    n = 1 << m
    return [sum(Y[j] * (1 - 2 * ((a & j).bit_count() % 2)) for j in range(n)) for a in range(n)]


def hadamard(m: int) -> list:
    n = 1 << m
    return [[1 - 2 * ((a & j).bit_count() % 2) for j in range(n)] for a in range(n)]


def kron(A, B) -> list:
    return [[a * b for a in ra for b in rb] for ra in A for rb in B]


def fht_additions(m: int) -> int:
    return m * (1 << m)


def direct_additions(m: int) -> int:
    return (1 << m) * ((1 << m) - 1)


def fht_decode(y, m: int):
    """(index of the largest |T|, a_0, the transform). Ties go to the lowest index."""
    T = fht_stages(signs(y))[-1]
    best = max(range(1 << m), key=lambda a: (abs(T[a]), -a))
    return best, (1 if T[best] < 0 else 0), T


def rm1_codeword(a0: int, a: int, m: int) -> list:
    return [a0 ^ ((a & j).bit_count() % 2) for j in range(1 << m)]


def word_error(n: int, t: int, p: float) -> float:
    return sum(comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(t + 1, n + 1))


def macwilliams(dist: dict, n: int) -> list:
    """Weight distribution of the dual code, by the MacWilliams identity with Krawtchouk polynomials."""
    size = sum(dist.values())
    B = []
    for j in range(n + 1):
        total = sum(a * sum((-1) ** l * comb(i, l) * comb(n - i, j - l) for l in range(j + 1)) for i, a in dist.items())
        if total % size:
            raise ValueError("MacWilliams sum not divisible by the code size")
        B.append(total // size)
    return B


# ---------------------------------------------------------- hardware models

def serial_encoder_rm13():
    """The serial RM(1,3) encoder as a netlist: gates (kind, input, input, output), evaluated once per clock."""
    return [("AND", "a1", "x1", "p1"), ("AND", "a2", "x2", "p2"), ("AND", "a3", "x3", "p3"),
            ("XOR", "a0", "p1", "s1"), ("XOR", "p2", "p3", "s2"), ("XOR", "s1", "s2", "c")]


def run_netlist(netlist, signals: dict) -> dict:
    sig = dict(signals)
    for kind, a, b, out in netlist:
        sig[out] = sig[a] & sig[b] if kind == "AND" else sig[a] ^ sig[b]
    return sig


def ripple(a: int, b: int, w: int, subtract: bool):
    """A w-bit ripple-carry adder or subtracter of full-adder cells on w-bit two's-complement patterns.
    The subtracter is the adder fed the inverted b and a carry-in of 1. Returns (w-bit result, carry out, cells)."""
    if subtract:
        b ^= (1 << w) - 1
    carry, out = int(subtract), 0
    for i in range(w):
        x, y = (a >> i) & 1, (b >> i) & 1
        out |= (x ^ y ^ carry) << i
        carry = (x & y) | (carry & (x ^ y))
    return out, carry, w


def to_bits(v: int, w: int) -> int:
    return v & ((1 << w) - 1)


def from_bits(p: int, w: int) -> int:
    return p - (1 << w) if (p >> (w - 1)) & 1 else p


def width_for(bound: int) -> int:
    """Fewest two's-complement bits that hold every whole number from -bound to bound."""
    w = 1
    while (1 << (w - 1)) - 1 < bound:
        w += 1
    return w


def butterfly_array(Y, carry_as_new_bit: bool = False):
    """The fast transform built from fixed-width two's-complement adders and subtracters.
    Each stage sign-extends its inputs to the width that holds the stage's range and drops the
    top carry; with carry_as_new_bit it adds in the old width and keeps the carry as the new top bit.
    Returns (correlations, widths after the sign map and each stage, adders and subtracters, cells)."""
    n, bound = len(Y), 1
    w = width_for(bound)
    vals = [to_bits(v, w) for v in Y]
    widths, units, cells = [w], 0, 0
    h = n // 2
    while h >= 1:
        bound *= 2
        wn = width_for(bound)
        z = list(vals)
        for j in range(n):
            if not j & h:
                a, b = vals[j], vals[j + h]
                if carry_as_new_bit:
                    s, cs, c1 = ripple(a, b, w, False)
                    d, cd, c2 = ripple(a, b, w, True)
                    z[j], z[j + h] = s | (cs << w), d | (cd << w)
                else:
                    a, b = to_bits(from_bits(a, w), wn), to_bits(from_bits(b, w), wn)
                    s, _, c1 = ripple(a, b, wn, False)
                    d, _, c2 = ripple(a, b, wn, True)
                    z[j], z[j + h] = s, d
                units += 2
                cells += c1 + c2
        vals, w = z, wn
        widths.append(w)
        h //= 2
    return [from_bits(v, w) for v in vals], widths, units, cells


def max_finder(T):
    """A tree of two-input comparators on |T| that carries the index and sign of the winner.
    Ties keep the left input, so the lowest index wins. Returns (index, a0, comparators, levels)."""
    layer = [(abs(t), i, 1 if t < 0 else 0) for i, t in enumerate(T)]
    comparators = levels = 0
    while len(layer) > 1:
        layer = [layer[p] if layer[p][0] >= layer[p + 1][0] else layer[p + 1] for p in range(0, len(layer), 2)]
        comparators += len(layer)
        levels += 1
    return layer[0][1], layer[0][2], comparators, levels


# ----------------------------------------------------------------- the checks

def check_truth_tables():
    section("Boolean functions as codewords, m = 3")
    same("point j = 3 as (x1, x2, x3)", vec((3 >> (i - 1)) & 1 for i in (1, 2, 3)), "1 1 0")
    same("point j = 6 as (x1, x2, x3)", vec((6 >> (i - 1)) & 1 for i in (1, 2, 3)), "0 1 1")
    tables = {S: word(truth_table(S, 3)) for S in monomials(3, 3)}
    for S, stated in [((), "11111111"), ((1,), "01010101"), ((2,), "00110011"), ((3,), "00001111"),
                      ((1, 2), "00010001"), ((1, 3), "00000101"), ((2, 3), "00000011"), ((1, 2, 3), "00000001")]:
        same(f"truth table of {mono_name(S)}", tables[S], stated)
    same("Boolean functions of 3 variables", 2 ** 8, 256)
    same("monomials in 3 variables", len(monomials(3, 3)), 8)
    same("codeword of 1 + x2 + x3 in RM(1,3)", word(encode([1, 0, 1, 1], 1, 3)), "11000011")
    same("weight of that codeword", sum(encode([1, 0, 1, 1], 1, 3)), 4)
    G = generator(1, 3)
    same("RM(1,3) generator, rows 1, x1, x2, x3", " ".join(word(r) for r in G), "11111111 01010101 00110011 00001111")
    same("codewords of RM(1,3)", len(rm_code(1, 3)), 16)
    # every function is the sum of the indicators of its 1 points, so the monomials span everything
    span3 = set(rm_code(3, 3))
    same("distinct words spanned by the 8 monomials", len(span3), 256)


def check_mobius():
    section("The Moebius transform")
    f = truth_table((1, 2), 3)
    f = add(f, truth_table((3,), 3))
    same("truth table of x3 + x1x2", word(f), "00011110")
    st = mobius_stages(f, 3)
    same("after the x1 stage", word(st[0]), "00011011")
    same("after the x2 stage", word(st[1]), "00011001")
    same("after the x3 stage, the coefficients", word(st[2]), "00011000")
    same("polynomial read from the coefficients", anf_terms(st[2], 3), "x3 + x1x2")
    # the transform undoes itself, and it agrees with the subset sums for every function of 3 variables
    ok_inv = ok_sum = True
    for x in range(256):
        t = from_int(x, 8)
        a = mobius_stages(t, 3)[-1]
        ok_inv &= mobius_stages(a, 3)[-1] == t
        ok_sum &= all(a[s] == sum(t[u] for u in range(8) if u & s == u) % 2 for s in range(8))
    require("transform applied twice returns all 256 tables", ok_inv)
    require("coefficient a_S equals the sum of f over subsets of S", ok_sum)
    same("XOR gates in the transform for m = 3", 3 * 4, 12)


def check_dimensions():
    section("Dimension and the special cases")
    for r, m, stated in [(1, 3, 4), (2, 3, 7), (1, 4, 5), (2, 4, 11), (1, 5, 6), (2, 5, 16), (3, 5, 26)]:
        same(f"k of RM({r},{m})", k_rm(r, m), stated)
    for r, m in [(r, m) for m in range(1, 6) for r in range(m + 1) if k_rm(r, m) <= 11]:
        require(f"RM({r},{m}) has 2^k distinct codewords, k = {k_rm(r, m)}", len(set(rm_code(r, m))) == 1 << k_rm(r, m))
    rows = {m: " ".join(str(k_rm(r, m)) for r in range(m + 1)) for m in range(1, 8)}
    for m, stated in [(1, "1 2"), (2, "1 3 4"), (3, "1 4 7 8"), (4, "1 5 11 15 16"), (5, "1 6 16 26 31 32"),
                      (6, "1 7 22 42 57 63 64"), (7, "1 8 29 64 99 120 127 128")]:
        same(f"k for r = 0..{m}, m = {m}", rows[m], stated)
    for m in range(1, 6):
        require(f"RM(0,{m}) is the repetition code", sorted(rm_code(0, m)) == [0, (1 << (1 << m)) - 1])
    even8 = sorted(c for c in range(256) if c.bit_count() % 2 == 0)
    require("RM(2,3) is the even-weight code of length 8", sorted(rm_code(2, 3)) == even8)
    require("RM(3,3) is every word of length 8", sorted(rm_code(3, 3)) == list(range(256)))
    for m in range(1, 8):
        require(f"k(m-1,{m}) = 2^m - 1", k_rm(m - 1, m) == (1 << m) - 1)


def check_recursion():
    section("The (u | u + v) recursion")
    for r, m in [(1, 2), (1, 3), (2, 3), (1, 4), (2, 4), (1, 5)]:
        n2 = 1 << (m - 1)
        U = rm_code(r, m - 1)
        V = rm_code(r - 1, m - 1)
        built = {u | ((u ^ v) << n2) for u in U for v in V}
        require(f"RM({r},{m}) = (u | u + v), u in RM({r},{m - 1}), v in RM({r - 1},{m - 1})", built == set(rm_code(r, m)))
    for m in range(2, 8):
        require(f"k(r,{m}) = k(r,{m - 1}) + k(r-1,{m - 1}) for 1 <= r <= {m - 1}",
                all(k_rm(r, m) == k_rm(r, m - 1) + k_rm(r - 1, m - 1) for r in range(1, m)))
    G13, G03 = generator(1, 3), generator(0, 3)
    G = [row + row for row in G13] + [[0] * 8 + G03[0]]
    same("rows of the recursive RM(1,4) generator", len(G), 5)
    same("its last row, x4", word(G[4]), "0000000011111111")
    require("it spans RM(1,4)", set(span_ints(G)) == set(rm_code(1, 4)))
    g = encode([1, 1, 0, 0], 1, 3)
    v = encode([1], 0, 3)
    c = g + add(g, v)
    same("u = 1 + x1 on m = 3", word(g), "10101010")
    same("codeword (u | u + v) with v = 1", word(c), "1010101001010101")
    same("the same word from 1 + x1 + x4", word(encode([1, 1, 0, 0, 1], 1, 4)), "1010101001010101")
    same("its weight", sum(c), 8)


def check_distance():
    section("Minimum distance")
    checked = 0
    for m in range(1, 8):
        for r in range(m + 1):
            if k_rm(r, m) <= 11:
                require(f"d of RM({r},{m}) = 2^{m - r} = {1 << (m - r)}", min_weight(rm_code(r, m)) == 1 << (m - r))
                checked += 1
    same("codes enumerated whole", checked, 18)
    # RM(2,5) is small enough to list; RM(3,5) and RM(4,5) get their weights from their duals
    d_of = {(r, m): min_weight(rm_code(r, m)) for r, m in [(1, 5), (2, 4), (1, 7), (2, 5)]}
    same("RM(2,5): words listed", len(rm_code(2, 5)), 65536)
    for r, m, stated in [(3, 5, "67108864 4"), (4, 5, "2147483648 2")]:
        B = macwilliams(weight_distribution(rm_code(m - r - 1, m)), 1 << m)
        d_of[(r, m)] = next(j for j in range(1, len(B)) if B[j])
        same(f"RM({r},{m}) by MacWilliams from RM({m - r - 1},{m}): words, and d", f"{sum(B)} {d_of[(r, m)]}", stated)
    for r, m, stated in [(1, 5, "32 6 16"), (2, 5, "32 16 8"), (3, 5, "32 26 4"), (4, 5, "32 31 2"),
                         (2, 4, "16 11 4"), (1, 7, "128 8 64")]:
        same(f"(n, k, d) of RM({r},{m})", vec([1 << m, k_rm(r, m), d_of[(r, m)]]), stated)
    same("(n, k, d) of RM(2,7), d from the formula", vec([1 << 7, k_rm(2, 7), 1 << (7 - 2)]), "128 29 32")
    for r, m, stated in [(1, 3, 1), (1, 4, 3), (2, 4, 1), (1, 5, 7), (2, 5, 3)]:
        d = d_of.get((r, m)) or min_weight(rm_code(r, m))
        same(f"guaranteed radius of RM({r},{m}), (d - 1) / 2 rounded down", (d - 1) // 2, stated)
    same("guaranteed radius of RM(3,6), from the formula", (1 << (6 - 3 - 1)) - 1, 3)


def check_first_order():
    section("First-order weights")
    for m in range(1, 8):
        stated = {0: 1, 1 << (m - 1): (1 << (m + 1)) - 2, 1 << m: 1}
        if m == 1:
            stated = {0: 1, 1: 2, 2: 1}
        require(f"RM(1,{m}): 1 at 0, 2^(m+1) - 2 at 2^(m-1), 1 at 2^m", weight_distribution(rm_code(1, m)) == stated)
    same("weight distribution of RM(1,3)", dist_text(weight_distribution(rm_code(1, 3))), "0:1 4:14 8:1")
    same("weight distribution of RM(1,4)", dist_text(weight_distribution(rm_code(1, 4))), "0:1 8:30 16:1")
    same("weight distribution of RM(1,5)", dist_text(weight_distribution(rm_code(1, 5))), "0:1 16:62 32:1")
    same("words of weight 32 in RM(1,6)", weight_distribution(rm_code(1, 6))[32], 126)
    # the extended Hamming code of length 8, built the way Hamming Codes and Perfect Codes builds it
    ext = []
    for x in range(256):
        c = from_int(x, 8)
        s = 0
        for j in range(1, 8):
            if c[j]:
                s ^= j
        if s == 0 and c[0] == sum(c[1:]) % 2:
            ext.append(x)
    same("weight distribution of the extended (8,4) Hamming code", dist_text(weight_distribution(ext)), "0:1 4:14 8:1")
    require("RM(1,3) is that code, position for position", sorted(ext) == sorted(rm_code(1, 3)))
    ok = all(sum(row[j] for row in generator(1, 3)) >= 1 for j in range(8))
    G = generator(1, 3)
    cols = [(G[0][j], G[1][j], G[2][j], G[3][j]) for j in range(8)]
    require("column j of the RM(1,3) generator is (1, bits of j)", ok and all(cols[j] == (1, j & 1, (j >> 1) & 1, (j >> 2) & 1) for j in range(8)))


def check_duality():
    section("Duality")
    for m in range(1, 7):
        for r in range(m):
            A, B = generator(r, m), generator(m - r - 1, m)
            ok = all(sum(a & b for a, b in zip(ra, rb)) % 2 == 0 for ra in A for rb in B)
            require(f"RM({r},{m}) is orthogonal to RM({m - r - 1},{m}), k sum = {k_rm(r, m) + k_rm(m - r - 1, m)}",
                    ok and k_rm(r, m) + k_rm(m - r - 1, m) == 1 << m)
    for m in range(1, 6):
        require(f"every monomial of degree below {m} has even weight, m = {m}",
                all(sum(truth_table(S, m)) % 2 == 0 for S in monomials(m - 1, m)))
    same("weight of x1x2 at m = 4", sum(truth_table((1, 2), 4)), 4)
    # RM(m-2,m) is the extended Hamming code: even weight, and the positions holding 1s XOR to zero
    for m in (3, 4):
        def ext_ok(c):
            s = 0
            for j in range(1 << m):
                if (c >> j) & 1:
                    s ^= j
            return s == 0 and c.bit_count() % 2 == 0
        words = rm_code(m - 2, m)
        require(f"every word of RM({m - 2},{m}) passes both extended Hamming checks", all(ext_ok(c) for c in words))
        same(f"words of length {1 << m} passing both checks, all {1 << (1 << m)} listed",
             sum(1 for c in range(1 << (1 << m)) if ext_ok(c)), len(words))
    same("weight distribution of RM(2,4)", dist_text(weight_distribution(rm_code(2, 4))), "0:1 4:140 6:448 8:870 10:448 12:140 16:1")
    # the 140 words of weight 4, counted by the recursion from u in RM(2,3) and v in RM(1,3)
    by_v = {0: 0, 4: 0, 8: 0}
    for u in rm_code(2, 3):
        for v in rm_code(1, 3):
            if u.bit_count() + (u ^ v).bit_count() == 4:
                by_v[v.bit_count()] += 1
    same("weight-4 words of RM(2,4) with v = 0, with wt(v) = 4, with v all ones", vec(by_v.values()), "28 112 0")
    inside = {sum(1 for u in rm_code(2, 3) if not u & ~v) for v in rm_code(1, 3) if v.bit_count() == 4}
    same("even-weight u inside the 1s of each v of weight 4", vec(sorted(inside)), "8")
    same("x1x2 at m = 3, then at m = 4", f"{word(truth_table((1, 2), 3))} {word(truth_table((1, 2), 4))}",
         "00010001 0001000100010001")
    xor0 = sum(1 for Q in combinations(range(16), 4) if Q[0] ^ Q[1] ^ Q[2] ^ Q[3] == 0)
    same("sets of 4 of the 16 points whose numbers XOR to 0", xor0, 140)
    w4 = [c for c in rm_code(2, 4) if c.bit_count() == 4]
    def xor_of(c):
        s = 0
        for j in range(16):
            if (c >> j) & 1:
                s ^= j
        return s
    require("the 140 weight-4 words of RM(2,4) are exactly those sets", len(w4) == 140 and all(xor_of(c) == 0 for c in w4))
    same("dimension of the dual of RM(1,5)", 32 - k_rm(1, 5), 26)


def check_reed():
    section("Reed's majority-logic decoder")
    # the subcube sums: over the subcube with free set S, a monomial of degree at most |S| sums to 1 exactly when it is x_S
    ok, cases = True, 0
    for m in range(1, 6):
        for s in range(m + 1):
            for S in combinations(range(1, m + 1), s):
                for b in range(1 << (m - s)):
                    Q = subcube(S, b, m)
                    for T in monomials(s, m):
                        tt = truth_table(T, m)
                        ok &= (sum(tt[p] for p in Q) % 2 == 1) == (T == S)
                        cases += 1
    require("subcube sum of x_T is 1 exactly when T = S, |T| <= |S|", ok)
    same("cases checked for m = 1 to 5", cases, 3880)
    same("check sums per degree-1 monomial of RM(1,3)", 1 << 2, 4)
    same("positions in the x1 check sums of RM(1,3)", " ".join("".join(str(p) for p in subcube((1,), b, 3)) for b in range(4)), "01 23 45 67")
    same("positions in the x2 check sums", " ".join("".join(str(p) for p in subcube((2,), b, 3)) for b in range(4)), "02 13 46 57")
    same("positions in the x3 check sums", " ".join("".join(str(p) for p in subcube((3,), b, 3)) for b in range(4)), "04 15 26 37")
    y = encode([1, 0, 1, 1], 1, 3)
    y[5] ^= 1
    same("received word, error in position 5", word(y), "11000111")
    coeffs, trace = reed_decode(y, 1, 3)
    same("check sums for x1", trace[0][1], "0010")
    same("check sums for x2", trace[1][1], "1110")
    same("check sums for x3", trace[2][1], "1011")
    same("word after removing x2 + x3", trace[3][1], "11111011")
    same("votes for a0, the word itself", trace[4][1], "11111011")
    same("decoded (a0, a1, a2, a3)", vec(coeffs[S] for S in monomials(1, 3)), "1 0 1 1")
    # RM(2,4): f = 1 + x2 + x4 + x1x2 + x2x3 + x3x4, one error in position 10
    f24 = [1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1]
    y = encode(f24, 2, 4)
    same("RM(2,4) codeword of 1 + x2 + x4 + x1x2 + x2x3 + x3x4", word(y), "1101111000101110")
    y[10] ^= 1
    same("RM(2,4) received, error in position 10", word(y), "1101111000001110")
    coeffs, trace = reed_decode(y, 2, 4)
    same("RM(2,4) check sums, x1x2 x1x3 x1x4 x2x3 x2x4 x3x4", " ".join(t[1] for t in trace[0:6]),
         "1101 0001 0100 1101 1000 1101")
    same("RM(2,4) degree-2 part, x1x2 + x2x3 + x3x4",
         word(add(add(truth_table((1, 2), 4), truth_table((2, 3), 4)), truth_table((3, 4), 4))), "0001001000011101")
    same("RM(2,4) word after removing the degree-2 part", trace[6][1], "1100110000010011")
    same("RM(2,4) check sums, x1 x2 x3 x4", " ".join(t[1] for t in trace[7:11]),
         "00000100 11110111 00000010 11011111")
    same("RM(2,4) word after removing x2 + x4", trace[11][1], "1111111111011111")
    same("RM(2,4) left after removing a0 = 1, the error", trace[13][1], "0000000000100000")
    same("RM(2,4) decoded coefficients in monomial order", vec(coeffs[S] for S in monomials(2, 4)),
         "1 0 1 0 1 1 0 0 1 0 1")
    # the int decoder and the list decoder agree on all 256 received words of length 8
    plan = reed_plan(1, 3)
    agree = True
    for x in range(256):
        c, _ = reed_decode(from_int(x, 8), 1, 3)
        a = reed_decode_int(x, plan)
        agree &= (a == -1) if c is None else (a == to_int(encode([c[S] for S in monomials(1, 3)], 1, 3)))
    require("bit-mask decoder matches the step-by-step one, all 256 words", agree)
    for r, m, sp, st, sf in [(1, 3, 9, 144, 0), (1, 4, 697, 22304, 0), (2, 4, 17, 34816, 0)]:
        npat, tried, failures = exhaustive_reed(r, m, True)
        same(f"RM({r},{m}): error patterns up to the radius", npat, sp)
        same(f"RM({r},{m}): decodings, every codeword times every pattern", tried, st)
        same(f"RM({r},{m}): decodings that went wrong", failures, sf)
    same("weight-7 patterns of length 32", comb(32, 7), 3365856)
    npat, tried, failures = exhaustive_reed(1, 5, False)
    same("RM(1,5): patterns of weight 0 to 7, all decoded", tried, 4514873)
    same("RM(1,5): decodings that went wrong", failures, 0)
    # the verdict depends only on the error: spread patterns of weight 0 to 9, past the radius
    # included, succeed or fail on each of the 64 codewords exactly as they do on the all-zero word
    plan = reed_plan(1, 5)
    words = rm_code(1, 5)
    sample = sorted({sum(1 << ((s + step * t) % 32) for t in range(w)) for step in (3, 5, 7) for w in range(10) for s in range(32)})
    same("distinct spread patterns, strides 3 5 7, weights 0 to 9, starts 0 to 31", len(sample), 801)
    on_zero = [reed_decode_int(e, plan) for e in sample]
    same("on the all-zero word: ties, wrong codewords, successes",
         vec([sum(1 for d in on_zero if d == -1), sum(1 for d in on_zero if d > 0), sum(1 for d in on_zero if d == 0)]), "96 96 609")
    def lemma(c, e, d):
        got = reed_decode_int(c ^ e, plan)
        return got == -1 if d == -1 else got == c ^ d
    require("D(c + e) = c + D(e) on all 64 codewords, a tie on one side exactly when on the other",
            all(lemma(c, e, d) for e, d in zip(sample, on_zero) for c in words))
    z = [0] * 8
    z[0] = z[1] = 1
    coeffs, trace = reed_decode(z, 1, 3)
    same("errors at 0 and 1: check sums for x1", trace[0][1], "0000")
    same("errors at 0 and 1: check sums for x2", trace[1][1], "1100")
    require("errors at 0 and 1: the x2 vote ties and the decoder stops", coeffs is None and len(trace) == 2)
    same("errors at 0 and 1: check sums for x3, had it gone on", word([z[j] ^ z[j + 4] for j in range(4)]), "1100")
    same("errors at 0 and 1: codewords at distance 2", sum(1 for c in rm_code(1, 3) if bin(c ^ to_int(z)).count("1") == 2), 4)
    near = {sum(1 for c in rm_code(1, 3) if (c ^ e).bit_count() == 2) for e in patterns(8, 2)}
    same("RM(1,3): codewords at distance 2 from each of the 28 double errors", vec(sorted(near)), "4")
    pairs = {sum(1 for c in rm_code(1, 3) if c.bit_count() == 4 and c & e == e) for e in patterns(8, 2)}
    same("RM(1,3): weight-4 codewords holding each pair of positions", vec(sorted(pairs)), "3")
    same("errors at 0 and 1: correlations", vec(fht_stages(signs(z))[-1]), "4 0 -4 0 -4 0 -4 0")


def check_fht():
    section("The fast Hadamard transform")
    H2 = hadamard(1)
    H4 = hadamard(2)
    same("H4, row by row", " | ".join(vec(r) for r in H4), "1 1 1 1 | 1 -1 1 -1 | 1 1 -1 -1 | 1 -1 -1 1")
    require("H4 = H2 (x) H2", kron(H2, H2) == H4)
    require("H8 = H2 (x) H4", kron(H2, H4) == hadamard(3))
    require("H32 = H2 (x) H16", kron(H2, hadamard(4)) == hadamard(5))
    y = encode([1, 0, 1, 1], 1, 3)
    y[5] ^= 1
    Y = signs(y)
    same("received word as signs", vec(Y), "-1 -1 1 1 1 -1 -1 -1")
    st = fht_stages(Y)
    same("after stage 1, pairs (j, j + 4)", vec(st[0]), "0 -2 0 0 -2 0 2 2")
    same("after stage 2, pairs (j, j + 2)", vec(st[1]), "0 -2 0 -2 0 2 -4 -2")
    same("after stage 3, pairs (j, j + 1), the transform T", vec(st[2]), "-2 2 -2 2 2 -2 -6 -2")
    require("T equals the direct correlation", st[2] == direct_correlation(Y, 3))
    best, a0, T = fht_decode(y, 3)
    same("index of the largest |T|", best, 6)
    same("its value", T[best], -6)
    same("decoded a0 from the sign", a0, 1)
    same("agreements with 1 + x2 + x3, (8 + 6) / 2", (8 + abs(T[best])) // 2, 7)
    ok = all(fht_stages(signs(from_int(x, 8)))[-1] == direct_correlation(signs(from_int(x, 8)), 3) for x in range(256))
    require("FHT equals direct correlation for all 256 words of length 8", ok)
    ok = all(2 ** 3 + T == 2 * sum(1 for j in range(8) if from_int(x, 8)[j] == rm1_codeword(0, a, 3)[j])
             for x in range(256) for a, T in enumerate(fht_stages(signs(from_int(x, 8)))[-1]))
    require("agreements = (2^m + T) / 2 for every word and every a", ok)
    for m, sf, sd in [(2, 8, 12), (3, 24, 56), (4, 64, 240), (5, 160, 992), (10, 10240, 1047552)]:
        same(f"m = {m}: additions, fast and direct", f"{fht_additions(m)} {direct_additions(m)}", f"{sf} {sd}")
    same("m = 10: ratio, direct over fast", f"{direct_additions(10) / fht_additions(10):.1f}", "102.3")
    same("butterflies for m = 3", sum(1 for h in (4, 2, 1) for j in range(8) if not j & h), 12)
    # seven errors on a Mariner codeword
    c = rm1_codeword(1, 13, 5)
    same("codeword with a0 = 1, (a1..a5) = 1 0 1 1 0", word(c), "10100101010110101010010101011010")
    errs = [2, 7, 11, 16, 21, 26, 30]
    y = list(c)
    for p in errs:
        y[p] ^= 1
    same("received word with errors at 2 7 11 16 21 26 30", word(y), "10000100010010100010000101111000")
    best, a0, T = fht_decode(y, 5)
    same("T over all 32 linear functions", vec(T),
         "10 -2 -2 2 -2 2 2 6 6 2 2 -2 2 -18 -2 -6 2 -2 -10 2 6 2 -6 -10 -2 2 -6 -2 -6 -2 6 -6")
    same("index of the largest |T|", best, 13)
    same("T there", T[13], -18)
    same("largest |T| anywhere else", max(abs(T[a]) for a in range(32) if a != 13), 10)
    same("indices where |T| = 10", vec(a for a in range(32) if abs(T[a]) == 10), "0 18 23")
    same("correct |T| with 7 errors, 32 - 14", 32 - 2 * 7, 18)
    same("bound on every other |T|, 2 * 7", 2 * 7, 14)
    same("decoded a0", a0, 1)
    # three errors on RM(1,4)
    c = rm1_codeword(0, 11, 4)
    y = list(c)
    for p in (3, 8, 14):
        y[p] ^= 1
    same("RM(1,4) codeword x1 + x2 + x4", word(c), "0110011010011001")
    same("received with errors at 3, 8, 14", word(y), "0111011000011011")
    st = fht_stages(signs(y))
    same("RM(1,4) after stage 1", vec(st[0]), "2 0 0 -2 0 0 -2 0 0 -2 -2 0 2 -2 0 2")
    same("RM(1,4) after stage 2", vec(st[1]), "2 0 -2 -2 2 0 2 -2 2 -4 -2 2 -2 0 -2 -2")
    same("RM(1,4) after stage 3", vec(st[2]), "0 -2 4 2 4 -2 0 2 0 -2 4 -6 -4 -2 0 2")
    same("RM(1,4) after stage 4, T", vec(st[3]), "-2 2 6 2 2 6 2 -2 -2 2 -2 10 -6 -2 2 -2")
    best, a0, T = fht_decode(y, 4)
    same("RM(1,4): index of the largest |T|", best, 11)
    same("RM(1,4): T there", T[11], 10)
    same("RM(1,4): received signs", vec(signs(y)), "1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 -1")
    same("RM(1,4): largest |T| anywhere else", max(abs(T[a]) for a in range(16) if a != 11), 6)
    same("RM(1,4): correct |T| with 3 errors, 16 - 6", 16 - 2 * 3, 10)


def check_mariner():
    section("The Mariner code")
    same("(n, k, d) of RM(1,5)", vec([32, k_rm(1, 5), 16]), "32 6 16")
    same("rate 6/32", f"{6 / 32:.4f}", "0.1875")
    same("word error, RM(1,5) to 7 errors, p = 0.05", f"{word_error(32, 7, 0.05):.6f}", "0.000139")
    same("word error, 6 bits uncoded, p = 0.05", f"{1 - 0.95 ** 6:.4f}", "0.2649")


def check_hardware():
    section("Hardware")
    # parallel encoder: c_j is the XOR of the a_S whose monomial is 1 at point j
    for r, m, stated in [(1, 3, 12), (1, 5, 80), (2, 4, 56)]:
        gates = sum(max(0, sum(1 for S in monomials(r, m) if truth_table(S, m)[j]) - 1) for j in range(1 << m))
        same(f"two-input XOR gates in the parallel RM({r},{m}) encoder", gates, stated)
    # serial encoder: the netlist of the figure, run for every message and every clock
    net = serial_encoder_rm13()
    ok = True
    for msg in range(16):
        a = [(msg >> i) & 1 for i in range(4)]
        out = [run_netlist(net, {"a0": a[0], "a1": a[1], "a2": a[2], "a3": a[3],
                                 "x1": j & 1, "x2": (j >> 1) & 1, "x3": (j >> 2) & 1})["c"] for j in range(8)]
        ok &= out == encode(a, 1, 3)
    require("serial RM(1,3) encoder netlist gives every codeword, 16 messages times 8 clocks", ok)
    same("serial RM(1,3) encoder: AND gates, XOR gates",
         f"{sum(1 for g in net if g[0] == 'AND')} {sum(1 for g in net if g[0] == 'XOR')}", "3 3")
    for r, m, stated in [(1, 3, 12), (1, 5, 80), (2, 4, 72)]:
        top = sum(len(subcube(S, b, m)) - 1 for S in combinations(range(1, m + 1), r) for b in range(1 << (m - r)))
        same(f"two-input XOR gates in the top-degree check sums, RM({r},{m})", top, stated)
    for m, stated in [(3, "3 4 1 8"), (5, "5 16 1 32")]:
        plan = reed_plan(1, m)
        same(f"RM(1,{m}): majority gates and inputs at degree 1, then at degree 0",
             f"{len(plan[0])} {plan[0][0][3]} {len(plan[1])} {plan[1][0][3]}", stated)
    for m, stated in [(3, 12), (5, 80)]:
        same(f"m = {m}: two-input XOR gates that remove the degree-1 part", sum(j.bit_count() for j in range(1 << m)), stated)
    # the butterfly array of fixed-width two's-complement adders, checked against the exact transform
    s00, c00, _ = ripple(0b01, 0b11, 2, False)
    same("1 + (-1) in 2 bits: sum bits, carry, and the 3-bit word read as a number",
         f"{s00:02b} {c00} {from_bits(s00 | (c00 << 2), 3)}", "00 1 -4")
    tests = {3: [signs(from_int(x, 8)) for x in range(256)],
             5: [signs(from_int(c, 32)) for c in rm_code(1, 5)]}
    c = rm1_codeword(1, 13, 5)
    for p in [2, 7, 11, 16, 21, 26, 30]:
        c[p] ^= 1
    tests[5].append(signs(c))
    for m, stated in [(3, "2 3 4 5 24 96"), (5, "2 3 4 5 6 7 160 800")]:
        exact = all(butterfly_array(Y)[0] == direct_correlation(Y, m) for Y in tests[m])
        require(f"m = {m}: fixed-width butterflies equal the exact transform on all {len(tests[m])} test words", exact)
        _, widths, units, cells = butterfly_array(tests[m][-1])
        same(f"m = {m}: widths after the sign map and each stage, adders and subtracters, full-adder cells",
             vec(widths + [units, cells]), stated)
        wrong = sum(1 for Y in tests[m] if butterfly_array(Y, True)[0] != direct_correlation(Y, m))
        same(f"m = {m}: test words wrong when the carry is taken as the new bit", wrong, {3: 256, 5: 65}[m])
    same("m = 5: butterflies in all stages", sum(1 for h in (16, 8, 4, 2, 1) for j in range(32) if not j & h), 80)
    T = butterfly_array(tests[5][-1])[0]
    best, a0, comparators, levels = max_finder(T)
    same("m = 5: max-finder on the seven-error word: index, a0, comparators, levels", vec([best, a0, comparators, levels]), "13 1 31 5")
    require("the max-finder agrees with the decoder", (best, a0) == fht_decode(c, 5)[:2])
    same("m = 3: comparators and levels in the max-finder", vec(max_finder(tests[3][0])[2:]), "7 3")
    for m, stated in [(3, "8 4"), (5, "32 6")]:
        top = max(max(abs(t) for t in butterfly_array(Y)[0]) for Y in tests[m])
        same(f"largest |T| for m = {m}, and the bits in its magnitude", f"{top} {top.bit_length()}", stated)


def check_exercises():
    section("Exercises")
    same("Exercise: (n, k, d) of RM(2,5), d from all its words", vec([32, k_rm(2, 5), min_weight(rm_code(2, 5))]), "32 16 8")
    same("Exercise: (n, k, d) of RM(3,6), d from the formula", vec([64, k_rm(3, 6), 1 << (6 - 3)]), "64 42 8")
    f = add(add(truth_table((), 3), truth_table((2,), 3)), truth_table((1, 3), 3))
    same("Exercise: truth table of 1 + x2 + x1x3", word(f), "11001001")
    maj = [1 if j.bit_count() >= 2 else 0 for j in range(8)]
    same("Exercise: truth table of the majority of 3 bits", word(maj), "00010111")
    same("Exercise: its polynomial", anf_terms(mobius_stages(maj, 3)[-1], 3), "x1x2 + x1x3 + x2x3")
    plan = reed_plan(2, 5)
    for level, stated in [(0, "10 8 4"), (1, "5 16 2"), (2, "1 32 1")]:
        same(f"Exercise: RM(2,5) degree {2 - level}: coefficients, check sums each, points per sum",
             f"{len(plan[level])} {plan[level][0][3]} {1 << len(plan[level][0][0])}", stated)
    same("Exercise: RM(2,5) radius, half the top-stage votes less 1", plan[0][0][3] // 2 - 1, 3)
    same("Exercise: words of weight 8 in RM(1,4)", weight_distribution(rm_code(1, 4))[8], 30)
    B = macwilliams(weight_distribution(rm_code(1, 5)), 32)
    same("Exercise: (n, k, d) of the dual of RM(1,5), d by MacWilliams",
         vec([32, k_rm(3, 5), next(j for j in range(1, 33) if B[j])]), "32 26 4")
    same("Exercise: (n, k, d) of RM(1,10), d from the formula", vec([1 << 10, k_rm(1, 10), 1 << 9]), "1024 11 512")
    same("Exercise: RM(1,10) radius, from the formula", (1 << 8) - 1, 255)


def main():
    check_truth_tables()
    check_mobius()
    check_dimensions()
    check_recursion()
    check_distance()
    check_first_order()
    check_duality()
    check_reed()
    check_fht()
    check_mariner()
    check_hardware()
    check_exercises()
    print("\nAll chapter 26 numbers reproduced.")


if __name__ == "__main__":
    main()
