"""Reproduce every number in chapter 34 of Error Correction Coding, "Ideals and Quotient Rings".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch34_ideals_and_quotient_rings.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees. All
the arithmetic is on whole numbers, so every comparison is exact.

Every finite ring the chapter uses is held as an addition table and a
multiplication table on the elements 0, 1, ..., n - 1, with 0 the zero element.
Polynomials over GF(2) are whole numbers whose bit i is the coefficient of x^i,
so a 3-bit word and a polynomial of degree below 3 are the same number, addition
is XOR, and multiplication is the carry-less product followed by reduction. The
quotients of GF(3)[x] and GF(5)[x] by a quadratic x^2 + k hold the element
c + d*alpha as the number c + q*d, and the output writes alpha as the letter a.

The checks enumerate rather than sample. Each ring is tested against the axioms
of a commutative ring with 1 on every triple of its elements. Its ideals are
found by testing every nonempty subset of its elements, and its units, zero
divisors and nilpotent elements by testing every pair. Every matching between
two of the four-element rings is found by trying every map of the four elements
to the four elements, keeping the 24 that are one-to-one, and testing all 16
sums and all 16 products under each. The reduction networks are checked on all
32 products of two 3-bit words, and the whole multiplier on all 64 pairs.

The module also holds the unit counts that ch34_charts.py draws, so the chart
and this check share one source.
"""

import sys


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


def joined(values) -> str:
    return " ".join(str(v) for v in values)


def yes(flag: bool) -> str:
    return "yes" if flag else "no"


# ---------------------------------------------------------------------------
# Polynomials over GF(2) as whole numbers: bit i is the coefficient of x^i.

def deg2(a: int) -> int:
    """The degree, with -1 standing for the zero polynomial."""
    return a.bit_length() - 1


def clmul(a: int, b: int) -> int:
    """The carry-less product: shift and XOR, so no carry ever crosses a column."""
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def mod2(f: int, g: int) -> int:
    while f and deg2(f) >= deg2(g):
        f ^= g << (deg2(f) - deg2(g))
    return f


def div2(f: int, g: int) -> int:
    q = 0
    while f and deg2(f) >= deg2(g):
        s = deg2(f) - deg2(g)
        q |= 1 << s
        f ^= g << s
    return q


def gcd2(a: int, b: int) -> int:
    while b:
        a, b = b, mod2(a, b)
    return a


def mulmod(a: int, b: int, p: int) -> int:
    return mod2(clmul(a, b), p)


def bits(a: int, width: int) -> str:
    return format(a, f"0{width}b")


def words(values, width: int) -> str:
    out = [bits(v, width) for v in values]
    return " ".join(out) if out else "none"


def weight(a: int) -> int:
    return bin(a).count("1")


def is_irreducible2(p: int) -> bool:
    """No divisor of degree 1 up to half the degree of p."""
    m = deg2(p)
    for d in range(2, 1 << (m // 2 + 1)):
        if mod2(p, d) == 0:
            return False
    return True


def euclid_remainders(a: int, b: int) -> list[int]:
    out = []
    while b:
        a, b = b, mod2(a, b)
        if b:
            out.append(b)
    return out


# ---------------------------------------------------------------------------
# Finite rings as tables on 0, 1, ..., n - 1.

def ring_zn(n: int):
    add = [[(a + b) % n for b in range(n)] for a in range(n)]
    mul = [[a * b % n for b in range(n)] for a in range(n)]
    return add, mul, 1 % n


def ring_mod2(p: int):
    """GF(2)[x]/(p(x)), each coset held by its remainder of degree below deg p."""
    n = 1 << deg2(p)
    add = [[a ^ b for b in range(n)] for a in range(n)]
    mul = [[mulmod(a, b, p) for b in range(n)] for a in range(n)]
    return add, mul, 1


def quad_mul(e: int, f: int, q: int, k: int) -> int:
    """(c + d alpha)(g + h alpha) with alpha^2 = -k, in GF(q)[x]/(x^2 + k)."""
    c, d = e % q, e // q
    g, h = f % q, f // q
    return (c * g - k * d * h) % q + q * ((c * h + d * g) % q)


def ring_quad(q: int, k: int):
    n = q * q
    add = [[(a % q + b % q) % q + q * ((a // q + b // q) % q) for b in range(n)] for a in range(n)]
    mul = [[quad_mul(a, b, q, k) for b in range(n)] for a in range(n)]
    return add, mul, 1


def qtxt(e: int, q: int) -> str:
    c, d = e % q, e // q
    if d == 0:
        return str(c)
    dpart = "a" if d == 1 else f"{d}a"
    return dpart if c == 0 else f"{c}+{dpart}"


def bitwise_ring():
    """The 2-bit words of Rings and Fields: XOR for addition, AND for multiplication, one = 11."""
    add = [[a ^ b for b in range(4)] for a in range(4)]
    mul = [[a & b for b in range(4)] for a in range(4)]
    return add, mul, 3


def f4_ring():
    """F_4 exactly as Rings and Fields gives it, on the symbols 0, 1, a, b held as 0, 1, 2, 3."""
    add = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    mul = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
    return add, mul, 1


def is_comm_ring(add, mul, one) -> bool:
    n = len(add)
    for a in range(n):
        if add[a][0] != a or mul[a][one] != a:
            return False
        if not any(add[a][b] == 0 for b in range(n)):
            return False
        for b in range(n):
            if add[a][b] != add[b][a] or mul[a][b] != mul[b][a]:
                return False
            for c in range(n):
                if add[add[a][b]][c] != add[a][add[b][c]]:
                    return False
                if mul[mul[a][b]][c] != mul[a][mul[b][c]]:
                    return False
                if mul[a][add[b][c]] != add[mul[a][b]][mul[a][c]]:
                    return False
    return True


def unit_list(mul, one) -> list[int]:
    n = len(mul)
    return [a for a in range(n) if any(mul[a][b] == one for b in range(n))]


def zero_divisor_list(mul) -> list[int]:
    n = len(mul)
    return [a for a in range(1, n) if any(mul[a][b] == 0 for b in range(1, n))]


def nilpotent_list(mul) -> list[int]:
    n = len(mul)
    out = []
    for a in range(1, n):
        p = a
        for _ in range(n):
            if p == 0:
                break
            p = mul[p][a]
        if p == 0:
            out.append(a)
    return out


def inverse_in(mul, one, a: int) -> int:
    for b in range(len(mul)):
        if mul[a][b] == one:
            return b
    return -1


def is_field(mul, one) -> bool:
    return len(unit_list(mul, one)) == len(mul) - 1


def characteristic(add, one) -> int:
    k, s = 1, one
    while s != 0:
        s = add[s][one]
        k += 1
    return k


def ideal_masks(add, mul) -> list[int]:
    """Every nonempty subset closed under addition and negation that absorbs multiplication."""
    n = len(add)
    neg = [next(b for b in range(n) if add[a][b] == 0) for a in range(n)]
    found = []
    for mask in range(1, 1 << n):
        members = [a for a in range(n) if mask >> a & 1]
        if not all(mask >> add[a][b] & 1 for a in members for b in members):
            continue
        if not all(mask >> neg[a] & 1 for a in members):
            continue
        if all(mask >> mul[r][a] & 1 for r in range(n) for a in members):
            found.append(mask)
    return found


def subgroup_masks(add) -> list[int]:
    n = len(add)
    neg = [next(b for b in range(n) if add[a][b] == 0) for a in range(n)]
    found = []
    for mask in range(1, 1 << n):
        members = [a for a in range(n) if mask >> a & 1]
        if all(mask >> add[a][b] & 1 for a in members for b in members) and \
                all(mask >> neg[a] & 1 for a in members):
            found.append(mask)
    return found


def members_of(mask: int) -> list[int]:
    return [a for a in range(mask.bit_length()) if mask >> a & 1]


def maximal_masks(masks: list[int], n: int) -> list[int]:
    full = (1 << n) - 1
    return [m for m in masks if m != full and
            not any(j != m and j != full and j & m == m for j in masks)]


def prime_masks(masks: list[int], mul) -> list[int]:
    n = len(mul)
    full = (1 << n) - 1
    out = []
    for m in masks:
        if m == full:
            continue
        if all((m >> a & 1) or (m >> b & 1) for a in range(n) for b in range(n) if m >> mul[a][b] & 1):
            out.append(m)
    return out


def generator_zn(mask: int, n: int) -> int:
    """The smallest positive element of an ideal of Z_n, or 0 for the zero ideal."""
    positive = [a for a in members_of(mask) if a > 0]
    return min(positive) if positive else 0


def quotient(add, mul, mask: int):
    """The cosets a + I, each named by its smallest element, with the tables they inherit.

    Returns the two tables, the coset of each element, and whether every choice of
    representatives gave the same coset for each sum and each product."""
    n = len(add)
    ideal = members_of(mask)
    coset_of = [-1] * n
    reps = []
    for a in range(n):
        if coset_of[a] < 0:
            for i in ideal:
                coset_of[add[a][i]] = len(reps)
            reps.append(a)
    k = len(reps)
    qadd = [[-1] * k for _ in range(k)]
    qmul = [[-1] * k for _ in range(k)]
    ok = True
    for a in range(n):
        for b in range(n):
            s, p = coset_of[add[a][b]], coset_of[mul[a][b]]
            ca, cb = coset_of[a], coset_of[b]
            if qadd[ca][cb] < 0:
                qadd[ca][cb] = s
            elif qadd[ca][cb] != s:
                ok = False
            if qmul[ca][cb] < 0:
                qmul[ca][cb] = p
            elif qmul[ca][cb] != p:
                ok = False
    return qadd, qmul, coset_of, ok


def matchings(ring_a, ring_b) -> list[list[int]]:
    """Every one-to-one map of a four-element ring onto another that keeps all sums and products."""
    add_a, mul_a, _ = ring_a
    add_b, mul_b, _ = ring_b
    found = []
    for m0 in range(4):
        for m1 in range(4):
            for m2 in range(4):
                for m3 in range(4):
                    m = [m0, m1, m2, m3]
                    if len(set(m)) < 4:
                        continue
                    if all(m[add_a[a][b]] == add_b[m[a]][m[b]] and m[mul_a[a][b]] == mul_b[m[a]][m[b]]
                           for a in range(4) for b in range(4)):
                        found.append(m)
    return found


def bijection_count() -> int:
    return sum(1 for m0 in range(4) for m1 in range(4) for m2 in range(4) for m3 in range(4)
               if len({m0, m1, m2, m3}) == 4)


# ---------------------------------------------------------------------------
# Unit counts of every quotient GF(2)[x]/(p(x)) with p of degree 1 to 6.

def unit_count2(p: int) -> int:
    return sum(1 for f in range(1, 1 << deg2(p)) if gcd2(p, f) == 1)


def unit_counts(top: int = 6) -> list[list[int]]:
    """For each degree m from 1 to top, the unit count of the quotient by each p of degree m."""
    return [[unit_count2(p) for p in range(1 << m, 1 << (m + 1))] for m in range(1, top + 1)]


# ---------------------------------------------------------------------------

def check_opening() -> None:
    section("The opening: two rules for multiplying 3-bit words")
    p, q = 0b1011, 0b1001
    same("raw product of 010 and 101", bits(clmul(2, 5), 4), "1010")
    same("010 * 101 with x^3 = x + 1", bits(mulmod(2, 5, p), 3), "001")
    add_p, mul_p, _ = ring_mod2(p)
    same("inverses of 001 to 111 with x^3 = x + 1",
         words([inverse_in(mul_p, 1, a) for a in range(1, 8)], 3), "001 101 110 111 010 011 100")
    same("nonzero words with an inverse, x^3 = x + 1", len(unit_list(mul_p, 1)), 7)
    same("x^3 + x + 1 at x = 0 and x = 1", joined([p & 1, weight(p) % 2]), "1 1")
    same("raw product of 011 and 111", bits(clmul(3, 7), 4), "1001")
    same("011 * 111 with x^3 = 1", bits(mulmod(3, 7, q), 3), "000")
    same("(x + 1)(x^2 + x + 1)", bits(clmul(3, 7), 4), "1001")
    _, mul_q, _ = ring_mod2(q)
    same("units with x^3 = 1", words(unit_list(mul_q, 1), 3), "001 010 100")
    same("their inverses", words([inverse_in(mul_q, 1, a) for a in unit_list(mul_q, 1)], 3), "001 100 010")
    same("zero divisors with x^3 = 1", words(zero_divisor_list(mul_q), 3), "011 101 110 111")
    require("each zero divisor times its partner is 000",
            all(mul_q[a][b] == 0 for a, b in ((3, 7), (5, 7), (6, 7), (7, 3))))
    require("both rules satisfy the ring axioms on all 512 triples",
            is_comm_ring(*ring_mod2(p)) and is_comm_ring(*ring_mod2(q)))


def check_ideals() -> None:
    section("Ideals")
    z12 = ring_zn(12)
    masks = ideal_masks(z12[0], z12[1])
    same("ideals of Z_12 among its 4095 nonempty subsets", len(masks), 6)
    same("their generators", joined(sorted(generator_zn(m, 12) for m in masks)), "0 1 2 3 4 6")
    same("their sizes", joined(sorted(len(members_of(m)) for m in masks)), "1 2 3 4 6 12")
    same("divisors of 12", joined(d for d in range(1, 13) if 12 % d == 0), "1 2 3 4 6 12")
    require("each ideal of Z_12 is the multiples of its generator",
            all(set(members_of(m)) == {generator_zn(m, 12) * r % 12 for r in range(12)} for m in masks))
    require("the additive subgroups of Z_12 are exactly its ideals", subgroup_masks(z12[0]) == masks)
    units12 = unit_list(z12[1], 1)
    same("units of Z_12", joined(units12), "1 5 7 11")
    require("an ideal of Z_12 holding a unit is all of Z_12",
            all(m == 4095 for m in masks if any(m >> u & 1 for u in units12)))
    same("kernel of reduction mod 12 within -36..36",
         joined(k for k in range(-36, 37) if k % 12 == 0), "-36 -24 -12 0 12 24 36")
    two = []
    agree = True
    for n in range(2, 13):
        add, mul, one = ring_zn(n)
        exactly_two = len(ideal_masks(add, mul)) == 2
        agree = agree and exactly_two == is_field(mul, one)
        if exactly_two:
            two.append(n)
    same("n from 2 to 12 with exactly two ideals in Z_n", joined(two), "2 3 5 7 11")
    require("two ideals exactly when Z_n is a field, n = 2..12", agree)
    mult = sorted({clmul(0b101, u) for u in range(4)})
    same("multiples of x^2 + 1 of degree below 4", words(mult, 4), "0000 0101 1010 1111")
    require("they are closed under addition", all(a ^ b in mult for a in mult for b in mult))


def check_principal() -> None:
    section("Principal ideals")
    f, g = 0b1001, 0b101
    same("gcd(x^3 + 1, x^2 + 1)", bits(gcd2(f, g), 2), "11")
    same("Euclid's remainders", joined(bits(r, 2) for r in euclid_remainders(f, g)), "11")
    same("1 (x^3 + 1) + x (x^2 + 1)", bits(clmul(1, f) ^ clmul(0b10, g), 2), "11")
    combos = {clmul(u, f) ^ clmul(v, g) for u in range(8) for v in range(16)}
    same("distinct u f + v g, deg u < 3, deg v < 4", len(combos), 32)
    require("they are the multiples of x + 1 of degree below 6",
            combos == {f2 for f2 in range(64) if mod2(f2, 0b11) == 0})
    constants = [2 * a0 for a0 in range(-2, 3) for a1 in range(-2, 3) for b0 in range(-2, 3) for b1 in range(-2, 3)]
    same("combinations 2a + xb tried in Z[x]", len(constants), 625)
    require("each has an even constant term, so none is 1", all(c % 2 == 0 for c in constants))
    same("remainder of the coefficient 1 of x on division by 2", 1 % 2, 1)


def check_quotients() -> None:
    section("Cosets and quotient rings")
    p = 0b111
    for r in range(4):
        same(f"the coset {bits(r, 2)} + I among degree below 4",
             words([f2 for f2 in range(16) if mod2(f2, p) == r], 4),
             ["0000 0111 1001 1110", "0001 0110 1000 1111", "0010 0101 1011 1100",
              "0011 0100 1010 1101"][r])
    same("x^3 * x^2 reduced mod x^2 + x + 1", bits(mulmod(0b1000, 0b100, p), 2), "11")
    same("1 * (x + 1) reduced mod x^2 + x + 1", bits(mulmod(1, 0b11, p), 2), "11")
    z12 = ring_zn(12)
    four = 1 | 1 << 4 | 1 << 8
    qadd, qmul, coset_of, ok = quotient(z12[0], z12[1], four)
    same("cosets of (4) in Z_12", len(qadd), 4)
    same("size of each coset", joined(sorted({coset_of.count(c) for c in range(4)})), "3")
    require("sums and products of cosets of (4) are well defined", ok)
    z4 = ring_zn(4)
    require("a + (4) -> a mod 4 matches Z_12/(4) with Z_4",
            all(qadd[a][b] == z4[0][a][b] and qmul[a][b] == z4[1][a][b] for a in range(4) for b in range(4)))
    add_w, mul_w, one_w = bitwise_ring()
    s = [0, 3]
    require("{00, 11} is a subring of the 2-bit words",
            all(add_w[a][b] in s and mul_w[a][b] in s for a in s for b in s) and one_w in s)
    same("01 AND 11, which leaves {00, 11}", bits(mul_w[1][3], 2), "01")
    cosets = [[0, 3], [1, 2]]
    which = [0, 1, 1, 0]
    sum_ok = sum(1 for c in cosets for d in cosets if len({which[add_w[a][b]] for a in c for b in d}) == 1)
    mul_ok = sum(1 for c in cosets for d in cosets if len({which[mul_w[a][b]] for a in c for b in d}) == 1)
    same("coset pairs with a well-defined sum", sum_ok, 4)
    same("coset pairs with a well-defined product", mul_ok, 1)
    same("01 AND 01 and 01 AND 10", joined([bits(mul_w[1][1], 2), bits(mul_w[1][2], 2)]), "01 00")
    gf9 = ring_quad(3, 1)
    same("elements of GF(2)[x]/(x^2 + x + 1)", len(ring_mod2(0b111)[0]), 4)
    same("elements of GF(2)[x]/(x^3 + x + 1)", len(ring_mod2(0b1011)[0]), 8)
    same("elements of GF(3)[x]/(x^2 + 1)", len(gf9[0]), 9)


FOUR_MODULI = [0b111, 0b110, 0b101, 0b100]
FOUR_NAMES = ["x^2 + x + 1", "x^2 + x", "x^2 + 1", "x^2"]


def check_four() -> None:
    section("The four quotients with four elements")
    stated = {
        0b111: ("00 10 11 01", "00 11 01 10", "01 10 11", "none", "none", 2),
        0b110: ("00 10 10 00", "00 11 00 11", "01", "10 11", "none", 4),
        0b101: ("00 10 01 11", "00 11 11 00", "01 10", "11", "11", 3),
        0b100: ("00 10 00 10", "00 11 10 01", "01 11", "10", "10", 3),
    }
    rings = {}
    for p, name in zip(FOUR_MODULI, FOUR_NAMES):
        add, mul, one = ring_mod2(p)
        rings[p] = (add, mul, one)
        xr, x1r, u, z, nil, ideals = stated[p]
        require(f"mod {name}: the ring axioms on all 64 triples", is_comm_ring(add, mul, one))
        same(f"mod {name}: the row of x", words(mul[2], 2), xr)
        same(f"mod {name}: the row of x + 1", words(mul[3], 2), x1r)
        same(f"mod {name}: units", words(unit_list(mul, one), 2), u)
        same(f"mod {name}: zero divisors", words(zero_divisor_list(mul), 2), z)
        same(f"mod {name}: nonzero nilpotents", words(nilpotent_list(mul), 2), nil)
        same(f"mod {name}: ideals", len(ideal_masks(add, mul)), ideals)
    same("one-to-one maps of four elements onto four", bijection_count(), 24)
    sym = ["0", "1", "a", "b"]
    found = matchings(rings[0b111], f4_ring())
    same("matchings of mod x^2 + x + 1 with F_4", len(found), 2)
    same("images of x and x + 1, first matching", joined([sym[found[0][2]], sym[found[0][3]]]), "a b")
    same("images of x and x + 1, second matching", joined([sym[found[1][2]], sym[found[1][3]]]), "b a")
    found = matchings(rings[0b110], bitwise_ring())
    same("matchings of mod x^2 + x with the 2-bit words", len(found), 2)
    same("image of 1 in both", joined(sorted({bits(m[1], 2) for m in found})), "11")
    same("images of x and x + 1, first matching", words([found[0][2], found[0][3]], 2), "01 10")
    same("images of x and x + 1, second matching", words([found[1][2], found[1][3]], 2), "10 01")
    counts = [len(matchings(rings[a], rings[b])) for i, a in enumerate(FOUR_MODULI) for b in FOUR_MODULI[i + 1:]]
    same("matchings between each pair of the four", joined(counts), "0 0 0 0 0 1")
    found = matchings(rings[0b101], rings[0b100])
    same("mod x^2 + 1 to mod x^2: images of x, x + 1", words([found[0][2], found[0][3]], 2), "11 10")
    z4 = ring_zn(4)
    same("matchings of Z_4 with each of the four", joined(len(matchings(z4, rings[p])) for p in FOUR_MODULI),
         "0 0 0 0")
    same("characteristic of Z_4", characteristic(z4[0], 1), 4)
    same("characteristic of each of the four",
         joined(characteristic(rings[p][0], 1) for p in FOUR_MODULI), "2 2 2 2")
    same("1 + 1 in Z_4", z4[0][1][1], 2)


def check_first_iso() -> None:
    section("The first isomorphism theorem")

    def ev3(f: int, a: int) -> int:
        return (f % 3 + (f // 3 % 3) * a + (f // 9) * a * a) % 3

    def mul3_linear(u: int) -> int:
        # (x + 1)(u0 + u1 x) = u0 + (u0 + u1) x + u1 x^2, as base-3 digits
        u0, u1 = u % 3, u // 3
        return u0 + 3 * ((u0 + u1) % 3) + 9 * u1

    kernel = {f for f in range(27) if ev3(f, 2) == 0}
    same("kernel of evaluation at 2, degree below 3 over GF(3)", len(kernel), 9)
    require("it is the multiples of x + 1 = x - 2 of degree below 3", kernel == {mul3_linear(u) for u in range(9)})
    fibers = [{f for f in range(27) if ev3(f, 2) == v} for v in range(3)]
    same("fiber sizes over the values 0, 1, 2", joined(len(fb) for fb in fibers), "9 9 9")
    require("each fiber is the coset v + kernel",
            all(fibers[v] == {(k % 3 + v) % 3 + 3 * (k // 3) for k in kernel} for v in range(3)))

    def pair(f: int) -> int:
        return 2 * (weight(f) % 2) + (f & 1)

    kern = [f for f in range(16) if pair(f) == 0]
    same("kernel of f -> (f(1), f(0)), degree below 4", words(kern, 4), "0000 0110 1010 1100")
    require("it is the multiples of x^2 + x of degree below 4", set(kern) == {clmul(0b110, u) for u in range(4)})
    same("values taken by f -> (f(1), f(0))", words(sorted({pair(f) for f in range(16)}), 2), "00 01 10 11")
    induced = [pair(r) for r in range(4)]
    same("images of 1, x and x + 1", words(induced[1:], 2), "11 10 01")
    require("the induced map is the second matching found",
            matchings(ring_mod2(0b110), bitwise_ring())[1] == induced)


def check_maximal() -> None:
    section("Maximal and prime ideals")
    z12 = ring_zn(12)
    masks = ideal_masks(z12[0], z12[1])
    same("generators of the maximal ideals of Z_12",
         joined(sorted(generator_zn(m, 12) for m in maximal_masks(masks, 12))), "2 3")
    same("generators of the prime ideals of Z_12",
         joined(sorted(generator_zn(m, 12) for m in prime_masks(masks, z12[1]))), "2 3")
    verdicts = []
    for d in (2, 3, 4, 6):
        mask = sum(1 << e for e in {d * r % 12 for r in range(12)})
        qadd, qmul, _, ok = quotient(z12[0], z12[1], mask)
        verdicts.append(f"{len(qadd)}:{yes(ok and is_field(qmul, 1))}")
    same("size and field verdict of Z_12/(d), d = 2, 3, 4, 6", joined(verdicts), "2:yes 3:yes 4:no 6:no")
    above = [m for m in masks if m & 273 == 273]
    z4 = ring_zn(4)
    same("ideals of Z_12 containing (4)", len(above), 3)
    same("ideals of Z_12/(4), which is Z_4", len(ideal_masks(z4[0], z4[1])), 3)
    r8 = ring_mod2(0b1001)
    m8 = ideal_masks(r8[0], r8[1])
    same("ideals of GF(2)[x]/(x^3 + 1)", len(m8), 4)
    same("their sizes", joined(sorted(len(members_of(m)) for m in m8)), "1 2 4 8")
    same("sizes of the maximal ones", joined(sorted(len(members_of(m)) for m in maximal_masks(m8, 8))), "2 4")
    same("the ideal (x + 1)", words(sorted({mulmod(0b11, u, 0b1001) for u in range(8)}), 3), "000 011 101 110")
    same("monic divisors of x^3 + 1", " ".join(format(d, "b") for d in range(1, 16) if mod2(0b1001, d) == 0),
         "1 11 111 1001")
    agree = True
    irreducibles = []
    for m in range(1, 5):
        for p in range(1 << m, 1 << (m + 1)):
            _, mul, one = ring_mod2(p)
            agree = agree and is_field(mul, one) == is_irreducible2(p)
            if is_irreducible2(p):
                irreducibles.append(format(p, "b"))
    require("field exactly when p is irreducible, all 30 of degree 1 to 4", agree)
    same("the irreducible p of degree 1 to 4", " ".join(irreducibles), "10 11 111 1011 1101 10011 11001 11111")
    same("(x^3 + x + 1) + x (x^2 + 1)", bits(0b1011 ^ clmul(0b10, 0b101), 1), "1")
    _, mul_p, _ = ring_mod2(0b1011)
    same("inverse of 101 modulo x^3 + x + 1", bits(inverse_in(mul_p, 1, 0b101), 3), "010")
    same("inverse of 010 modulo x^3 + x + 1", bits(inverse_in(mul_p, 1, 0b010), 3), "101")
    _, mul_q, _ = ring_mod2(0b1001)
    same("units and zero divisors modulo x^3 + 1",
         joined([len(unit_list(mul_q, 1)), len(zero_divisor_list(mul_q))]), "3 4")
    same("words divisible by x + 1, by x^2 + x + 1, by both",
         joined([sum(1 for f in range(8) if mod2(f, 3) == 0), sum(1 for f in range(8) if mod2(f, 7) == 0),
                 sum(1 for f in range(8) if mod2(f, 3) == 0 and mod2(f, 7) == 0)]), "4 2 1")
    same("8 - 4 - 2 + 1", 8 - 4 - 2 + 1, 3)
    counts = unit_counts(6)
    same("moduli of each degree 1 to 6", joined(len(c) for c in counts), "2 4 8 16 32 64")
    same("moduli reaching 2^m - 1 units, each degree",
         joined(sum(1 for u in c if u == (1 << m) - 1) for m, c in enumerate(counts, start=1)), "2 1 2 3 6 9")
    require("2^m - 1 units exactly when p is irreducible, degree 1 to 6",
            all((counts[m - 1][p - (1 << m)] == (1 << m) - 1) == is_irreducible2(p)
                for m in range(1, 7) for p in range(1 << m, 1 << (m + 1))))
    six = counts[5]
    same("most units of a reducible p of degree 6",
         max(u for p, u in zip(range(64, 128), six) if not is_irreducible2(p)), 56)
    same("fewest units of any p of degree 6", min(six), 12)
    require("the gcd count agrees with the pair count, degree 1 to 4",
            all(counts[m - 1][p - (1 << m)] == len(unit_list(ring_mod2(p)[1], 1))
                for m in range(1, 5) for p in range(1 << m, 1 << (m + 1))))


def check_root() -> None:
    section("A root by construction")
    same("alpha^2 in GF(2)[x]/(x^2 + x + 1)", bits(mulmod(2, 2, 0b111), 2), "11")
    same("alpha^2 + alpha + 1 there", bits(mulmod(2, 2, 0b111) ^ 2 ^ 1, 2), "00")
    same("x^2 + 1 at 0, 1, 2 over GF(3)", joined((a * a + 1) % 3 for a in range(3)), "1 2 2")
    same("squares in GF(3)", joined(sorted({a * a % 3 for a in range(3)})), "0 1")
    gf9 = ring_quad(3, 1)
    require("GF(3)[x]/(x^2 + 1): the ring axioms on all 729 triples", is_comm_ring(*gf9))
    same("its units, by all 81 pairs", len(unit_list(gf9[1], 1)), 8)
    same("alpha^2 there", qtxt(gf9[1][3][3], 3), "2")
    same("(1 + 2 alpha)(2 + alpha)", qtxt(gf9[1][1 + 3 * 2][2 + 3], 3), "2a")
    same("inverse of 1 + alpha", qtxt(inverse_in(gf9[1], 1, 1 + 3), 3), "2+a")
    f = clmul(0b111, 0b1011)
    same("(x^2 + x + 1)(x^3 + x + 1)", format(f, "b"), "110001")
    same("x^5 + x^4 + 1 at x = 0 and x = 1", joined([f & 1, weight(f) % 2]), "1 1")
    same("x^5 + x^4 + 1 reduced mod x^2 + x + 1", mod2(f, 0b111), 0)
    same("x^5 + x^4 + 1 reduced mod x^3 + x + 1", mod2(f, 0b1011), 0)


def network(p: int):
    """The taps of the reduction network for a product of two words of degree below deg p."""
    m = deg2(p)
    images = [mod2(1 << k, p) for k in range(2 * m - 1)]
    taps = [[k for k in range(2 * m - 1) if images[k] >> j & 1] for j in range(m)]
    return images, taps


def run_network(c: int, taps) -> int:
    r = 0
    for j, t in enumerate(taps):
        bit = 0
        for k in t:
            bit ^= c >> k & 1
        r |= bit << j
    return r


def check_hardware() -> None:
    section("Hardware")
    same("XOR gates to add two 3-bit words", 3, 3)
    same("AND gates in the 3-bit product", 3 * 3, 9)
    columns = [sum(1 for i in range(3) for j in range(3) if i + j == k) for k in range(5)]
    same("partial products in columns 0 to 4", joined(columns), "1 2 3 2 1")
    same("XOR gates to sum the columns", sum(c - 1 for c in columns), 4)
    stated = {0b1011: ("011 110", "c0 c3|c1 c3 c4|c2 c4", 4), 0b1001: ("001 010", "c0 c3|c1 c4|c2", 2)}
    for p, name in ((0b1011, "x^3 + x + 1"), (0b1001, "x^3 + 1")):
        images, taps = network(p)
        img, tp, gates = stated[p]
        same(f"x^3 and x^4 modulo {name}", words(images[3:], 3), img)
        same(f"taps of r0|r1|r2 modulo {name}", "|".join(" ".join(f"c{k}" for k in t) for t in taps), tp)
        same(f"XOR gates in the reduction modulo {name}", sum(len(t) - 1 for t in taps), gates)
        require(f"the network matches reduction on all 32 products, {name}",
                all(run_network(c, taps) == mod2(c, p) for c in range(32)))
        require(f"AND array, XOR trees and network on all 64 pairs, {name}",
                all(run_network(clmul(a, b), taps) == mulmod(a, b, p) for a in range(8) for b in range(8)))
    same("XOR gates in the whole multiplier, x^3 + x + 1", 4 + 4, 8)
    same("XOR gates in the whole multiplier, x^3 + 1", 4 + 2, 6)


def check_worked() -> None:
    section("Worked examples")
    z8 = ring_zn(8)
    m8 = ideal_masks(z8[0], z8[1])
    same("generators of the ideals of Z_8", joined(sorted(generator_zn(m, 8) for m in m8)), "0 1 2 4")
    same("generators of its maximal ideals", joined(generator_zn(m, 8) for m in maximal_masks(m8, 8)), "2")
    qadd, qmul, _, ok = quotient(z8[0], z8[1], 1 | 1 << 2 | 1 << 4 | 1 << 6)
    same("Z_8/(2): size and field", joined([len(qadd), yes(ok and is_field(qmul, 1))]), "2 yes")
    qadd, qmul, _, ok = quotient(z8[0], z8[1], 1 | 1 << 4)
    same("Z_8/(4): size and field", joined([len(qadd), yes(ok and is_field(qmul, 1))]), "4 no")
    p = 0b10011
    same("Euclid's remainders for x^4 + x + 1 and x^2 + 1",
         joined(format(r, "b") for r in euclid_remainders(p, 0b101)), "10 1")
    same("quotients at the two steps", joined([format(div2(p, 0b101), "b"), format(div2(0b101, 0b10), "b")]),
         "101 10")
    same("(x^2 + 1)(x^3 + x + 1) + x (x^4 + x + 1)", clmul(0b101, 0b1011) ^ clmul(0b10, p), 1)
    same("inverse of x^2 + 1 modulo x^4 + x + 1", bits(inverse_in(ring_mod2(p)[1], 1, 0b101), 4), "1011")
    q = 0b10101
    same("(x^2 + x + 1)^2", format(clmul(0b111, 0b111), "b"), "10101")
    _, mul_q, _ = ring_mod2(q)
    same("units, zero divisors modulo x^4 + x^2 + 1",
         joined([len(unit_list(mul_q, 1)), len(zero_divisor_list(mul_q))]), "12 3")
    same("its zero divisors", words(zero_divisor_list(mul_q), 4), "0111 1001 1110")
    same("its nonzero nilpotents", words(nilpotent_list(mul_q), 4), "0111 1001 1110")
    same("squares in GF(5)", joined(sorted({a * a % 5 for a in range(5)})), "0 1 4")
    same("-2 in GF(5)", -2 % 5, 3)
    gf25 = ring_quad(5, 2)
    require("GF(5)[x]/(x^2 + 2): the ring axioms on all 15625 triples", is_comm_ring(*gf25))
    same("its units, by all 625 pairs", len(unit_list(gf25[1], 1)), 24)
    same("alpha^2 there", qtxt(gf25[1][5][5], 5), "3")
    same("inverse of alpha", qtxt(inverse_in(gf25[1], 1, 5), 5), "2a")
    same("(1 + alpha)^2", qtxt(gf25[1][6][6], 5), "4+2a")
    z12, z4 = ring_zn(12), ring_zn(4)
    require("a -> a mod 4 keeps every sum and product of Z_12",
            all((z12[0][a][b]) % 4 == z4[0][a % 4][b % 4] and (z12[1][a][b]) % 4 == z4[1][a % 4][b % 4]
                for a in range(12) for b in range(12)))
    same("its kernel", joined(a for a in range(12) if a % 4 == 0), "0 4 8")
    same("its fiber sizes", joined(sum(1 for a in range(12) if a % 4 == v) for v in range(4)), "3 3 3 3")


def check_exercises() -> None:
    section("Exercises")
    z10 = ring_zn(10)
    m10 = ideal_masks(z10[0], z10[1])
    same("generators of the ideals of Z_10", joined(sorted(generator_zn(m, 10) for m in m10)), "0 1 2 5")
    same("generators of its maximal ideals", joined(sorted(generator_zn(m, 10) for m in maximal_masks(m10, 10))),
         "2 5")
    require("f(1) = 0 exactly for the multiples of x + 1, degree below 6",
            all((weight(f) % 2 == 0) == (mod2(f, 0b11) == 0) for f in range(64)))
    same("inverse of x modulo x^4 + x + 1", bits(inverse_in(ring_mod2(0b10011)[1], 1, 0b10), 4), "1001")
    same("x^4 + x^3 + x^2 + x + 1 irreducible", yes(is_irreducible2(0b11111)), "yes")
    same("units modulo x^4 + x^3 + x^2 + x + 1", len(unit_list(ring_mod2(0b11111)[1], 1)), 15)
    same("(x + 1)^4", format(clmul(clmul(3, 3), clmul(3, 3)), "b"), "10001")
    _, mul_r, _ = ring_mod2(0b10001)
    same("units and nonzero nilpotents modulo x^4 + 1",
         joined([len(unit_list(mul_r, 1)), len(nilpotent_list(mul_r))]), "8 7")
    same("roots of x^2 + 1 in GF(5)", joined(a for a in range(5) if (a * a + 1) % 5 == 0), "2 3")
    gf5 = ring_quad(5, 1)
    same("(x + 2)(x + 3) modulo x^2 + 1 over GF(5)", qtxt(gf5[1][2 + 5][3 + 5], 5), "0")
    same("units of GF(5)[x]/(x^2 + 1)", len(unit_list(gf5[1], 1)), 16)
    z12 = ring_zn(12)

    def ideal(d: int) -> int:
        return sum(1 << e for e in {d * r % 12 for r in range(12)})

    same("(4) and (6) in Z_12 share", joined(members_of(ideal(4) & ideal(6))), "0")
    same("(2) and (3) in Z_12 share", joined(members_of(ideal(2) & ideal(3))), "0 6")
    require("both intersections are ideals",
            ideal(4) & ideal(6) in ideal_masks(z12[0], z12[1]) and ideal(2) & ideal(3) in ideal_masks(z12[0], z12[1]))


def main() -> None:
    check_opening()
    check_ideals()
    check_principal()
    check_quotients()
    check_four()
    check_first_iso()
    check_maximal()
    check_root()
    check_hardware()
    check_worked()
    check_exercises()
    print("\nAll chapter 34 numbers reproduced.")


if __name__ == "__main__":
    main()
