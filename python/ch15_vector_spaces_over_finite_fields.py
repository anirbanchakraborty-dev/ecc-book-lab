"""Reproduce every number in chapter 15 of Error Correction Coding, "Vector Spaces over Finite Fields".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch15_vector_spaces_over_finite_fields.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
A vector over GF(q), for q = 2, 3 or 5, is held as a tuple of residues, and
the small spaces are checked exhaustively: the vector space axioms on every
triple of vectors, the subspace test on every subset of GF(2)^3, every
subspace of GF(2)^4, GF(3)^3 and GF(5)^2 built and counted by dimension, and
every dual found by testing each vector of the space against the code. The
only numbers in the chapter that are not whole are the fractions of lists of
vectors that form a basis, which are products of ten factors and are printed
to three decimal places. A run that reaches the final line has reproduced
every number in the chapter.

The module also holds the counting functions that ch15_charts.py draws, so the
figures and this check share one source.
"""

import sys
from itertools import combinations, product
from math import comb, gcd, prod


# Helpers ------------------------------------------------------------------------

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


def vec(s: str) -> tuple:
    """A vector written as a string of digits, such as '0110' or '120'."""
    return tuple(int(ch) for ch in s)


def show(v: tuple, q: int) -> str:
    if q == 2:
        return "".join(str(a) for a in v)
    return "(" + ",".join(str(a) for a in v) + ")"


def as_set(vs, q: int) -> str:
    return "{" + ", ".join(show(v, q) for v in sorted(vs)) + "}"


def as_row(vs, q: int) -> str:
    return " ".join(show(v, q) for v in vs)


def as_ints(xs) -> str:
    return " ".join(str(x) for x in xs)


def weight(v: tuple) -> int:
    return sum(1 for a in v if a != 0)


# Vectors over GF(q) --------------------------------------------------------------

def space(q: int, n: int) -> list:
    """Every vector of GF(q)^n, in dictionary order."""
    return list(product(range(q), repeat=n))


def zero(n: int) -> tuple:
    return (0,) * n


def add(u: tuple, v: tuple, q: int) -> tuple:
    return tuple((a + b) % q for a, b in zip(u, v))


def scale(c: int, v: tuple, q: int) -> tuple:
    return tuple((c * a) % q for a in v)


def dot(u: tuple, v: tuple, q: int) -> int:
    return sum(a * b for a, b in zip(u, v)) % q


def combine(coeffs, vs, q: int, n: int) -> tuple:
    out = zero(n)
    for c, v in zip(coeffs, vs):
        out = add(out, scale(c, v, q), q)
    return out


def span(vs, q: int, n: int) -> list:
    """Every linear combination of the vectors vs, sorted."""
    out = {zero(n)}
    for v in vs:
        out = {add(s, scale(c, v, q), q) for s in out for c in range(q)}
    return sorted(out)


def independent(vs, q: int, n: int) -> bool:
    """True when only the all-zero choice of coefficients gives the zero vector."""
    k = len(vs)
    return all(combine(cs, vs, q, n) != zero(n) for cs in product(range(q), repeat=k) if any(cs))


def is_subspace(S, q: int) -> bool:
    members = set(S)
    if not members:
        return False
    n = len(next(iter(members)))
    if zero(n) not in members:
        return False
    return all(add(u, v, q) in members for u in members for v in members) and \
        all(scale(c, u, q) in members for c in range(q) for u in members)


def dimension(S, q: int) -> int:
    k, size = 0, 1
    while size < len(S):
        k, size = k + 1, size * q
    return k if size == len(S) else -1


def dual(S, q: int, n: int) -> list:
    """Every vector of GF(q)^n whose dot product with each member of S is 0."""
    return [y for y in space(q, n) if all(dot(y, s, q) == 0 for s in S)]


def subspaces(q: int, n: int) -> list:
    """All subspaces of GF(q)^n, grouped by dimension, each a sorted tuple of vectors."""
    levels = [[(zero(n),)]]
    for _ in range(n):
        found = set()
        for S in levels[-1]:
            members = set(S)
            for v in space(q, n):
                if v not in members:
                    found.add(tuple(span(list(S) + [v], q, n)))
        levels.append(sorted(found))
    return levels


def lists_independent(q: int, n: int, k: int) -> int:
    """Count ordered lists of k vectors of GF(q)^n that are linearly independent, by trying every list."""
    return sum(1 for vs in product(space(q, n), repeat=k) if independent(vs, q, n))


# Counting formulas ---------------------------------------------------------------

def ordered_independent(q: int, n: int, k: int) -> int:
    """(q^n - 1)(q^n - q)...(q^n - q^(k-1)), the ordered independent lists of k vectors in GF(q)^n."""
    out = 1
    for i in range(k):
        out *= q ** n - q ** i
    return out


def gaussian(n: int, k: int, q: int) -> int:
    """The number of k-dimensional subspaces of GF(q)^n."""
    return ordered_independent(q, n, k) // ordered_independent(q, k, k)


def basis_fraction(q: int, n: int) -> float:
    """The fraction of ordered lists of n vectors of GF(q)^n that are bases."""
    p = 1.0
    for i in range(n):
        p = p * (1 - 1 / q ** (n - i))
    return p


def fraction_text(num: int, den: int) -> str:
    g = gcd(num, den)
    return f"{num // g}/{den // g}"


# The codes of the chapter ----------------------------------------------------------

C5 = [vec(s) for s in ("00000", "01101", "10110", "11011")]
R3 = [vec("000"), vec("111")]
R5 = [vec("00000"), vec("11111")]
CHECKS5 = [vec("11100"), vec("10010"), vec("01001")]


def syndrome(x: tuple) -> tuple:
    """The syndrome of Groups and Cosets: (x1 + x2 + x3, x1 + x4, x2 + x5) over GF(2)."""
    return ((x[0] + x[1] + x[2]) % 2, (x[0] + x[3]) % 2, (x[1] + x[4]) % 2)


def even_words(n: int) -> list:
    return [v for v in space(2, n) if weight(v) % 2 == 0]


def weight_counts(S) -> str:
    counts = {}
    for v in S:
        counts[weight(v)] = counts.get(weight(v), 0) + 1
    return " ".join(f"{w}:{counts[w]}" for w in sorted(counts))


def axioms_hold(q: int, n: int) -> bool:
    V = space(q, n)
    z = zero(n)
    for u in V:
        if add(u, z, q) != u or scale(1, u, q) != u:
            return False
        if not any(add(u, w, q) == z for w in V):
            return False
        for v in V:
            if add(u, v, q) != add(v, u, q):
                return False
            for c in range(q):
                if scale(c, add(u, v, q), q) != add(scale(c, u, q), scale(c, v, q), q):
                    return False
            for w in V:
                if add(add(u, v, q), w, q) != add(u, add(v, w, q), q):
                    return False
        for c in range(q):
            for d in range(q):
                if scale((c + d) % q, u, q) != add(scale(c, u, q), scale(d, u, q), q):
                    return False
                if scale((c * d) % q, u, q) != scale(c, scale(d, u, q), q):
                    return False
    return True


def self_dual_count(n: int) -> int:
    """Binary self-dual codes of length n: ordered bases of n/2 pairwise orthogonal even words, over the bases per code."""
    k = n // 2
    candidates = [v for v in space(2, n) if any(v) and dot(v, v, 2) == 0]

    def extend(chosen: list) -> int:
        if len(chosen) == k:
            return 1
        members = set(span(chosen, 2, n))
        return sum(extend(chosen + [v]) for v in candidates
                   if v not in members and all(dot(v, u, 2) == 0 for u in chosen))

    return extend([]) // ordered_independent(2, k, k)


def has_information_set(S, n: int) -> bool:
    """True when some k = dim S positions take every pattern exactly once across S."""
    k = dimension(S, 2)
    for J in combinations(range(n), k):
        if len({tuple(v[j] for j in J) for v in S}) == len(S):
            return True
    return False


# The checks ------------------------------------------------------------------------

def check_scaling() -> None:
    section("Bit strings that can be scaled")
    same("vectors in GF(2)^3", len(space(2, 3)), 8)
    same("0 times 101", show(scale(0, vec("101"), 2), 2), "000")
    same("1 times 101", show(scale(1, vec("101"), 2), 2), "101")
    same("101 + 011 over GF(2)", show(add(vec("101"), vec("011"), 2), 2), "110")
    require("v + v = 0 for every v in GF(2)^3", all(add(v, v, 2) == zero(3) for v in space(2, 3)))
    same("vectors in GF(5)^2", len(space(5, 2)), 25)
    same("vectors in GF(3)^3", len(space(3, 3)), 27)
    same("(1,2) + (3,4) over GF(5)", show(add(vec("12"), vec("34"), 5), 5), "(4,1)")
    same("3 times (1,2) over GF(5)", show(scale(3, vec("12"), 5), 5), "(3,1)")
    same("(1,2) + (1,2) over GF(5)", show(add(vec("12"), vec("12"), 5), 5), "(2,4)")
    same("c times (1,2) for c = 0, 1, 2, 3, 4", as_row([scale(c, vec("12"), 5) for c in range(5)], 5),
         "(0,0) (1,2) (2,4) (3,1) (4,3)")
    same("3 times (1,2) over the integers", show(tuple(3 * a for a in vec("12")), 5), "(3,6)")
    same("(1,2) + (3,4) over the integers", show(tuple(a + b for a, b in zip(vec("12"), vec("34"))), 5), "(4,6)")
    ok = True
    for v in space(5, 2):
        if any(v):
            multiples = {scale(c, v, 5) for c in range(1, 5)}
            ok = ok and len(multiples) == 4 and zero(2) not in multiples
    require("each nonzero word of GF(5)^2 has four different nonzero multiples", ok)
    require("a nonzero binary word has one nonzero multiple, itself",
            all({scale(c, v, 2) for c in range(2)} - {zero(3)} == {v} for v in space(2, 3) if any(v)))
    E = even_words(3)
    same("the even-weight words of length 3", as_set(E, 2), "{000, 011, 101, 110}")
    require("they are closed under addition", all(add(u, v, 2) in E for u in E for v in E))
    same("the coset 001 + E", as_set([add(vec("001"), v, 2) for v in E], 2), "{001, 010, 100, 111}")
    require("every edge of the cube joins an even-weight word to an odd one",
            all(weight(u) % 2 != weight(v) % 2 for u in space(2, 3) for v in space(2, 3) if weight(add(u, v, 2)) == 1))


def check_fields() -> None:
    section("Fields as scalars")
    ok = True
    for n in range(2, 13):
        has_inverses = all(any(a * b % n == 1 for b in range(n)) for a in range(1, n))
        is_prime = all(n % d != 0 for d in range(2, n))
        ok = ok and has_inverses == is_prime
    require("Z_n is a field exactly when n is prime, n = 2 to 12", ok)
    same("an inverse of 2 modulo 4 exists", any(2 * b % 4 == 1 for b in range(4)), False)
    same("2 times 2 modulo 4", 2 * 2 % 4, 0)
    ok = True
    for q in (2, 3, 5):
        R = range(q)
        ok = ok and all(a * (b + c) % q == (a * b + a * c) % q and (a + b) * c % q == (a * c + b * c) % q
                        for a in R for b in R for c in R)
    require("both distributive laws in GF(2), GF(3) and GF(5)", ok)
    require("0 a = 0 for every a in GF(2), GF(3) and GF(5)", all(0 * a % q == 0 for q in (2, 3, 5) for a in range(q)))
    same("the negative of 1 in GF(2)", (-1) % 2, 1)
    same("the negative of 2 in GF(5)", (-2) % 5, 3)


def check_axioms() -> None:
    section("The vector space axioms")
    require("the eight axioms on every triple of GF(2)^3", axioms_hold(2, 3))
    require("the eight axioms on every triple of GF(3)^2", axioms_hold(3, 2))
    require("the eight axioms on every triple of GF(5)^2", axioms_hold(5, 2))
    V = space(5, 2)
    require("0 v = 0 for every v in GF(5)^2", all(scale(0, v, 5) == zero(2) for v in V))
    require("c 0 = 0 for every c in GF(5)", all(scale(c, zero(2), 5) == zero(2) for c in range(5)))
    require("(-1) v is the negative of v for every v in GF(5)^2", all(add(v, scale(4, v, 5), 5) == zero(2) for v in V))
    same("the negative of (1,2) in GF(5)^2", show(scale(4, vec("12"), 5), 5), "(4,3)")
    same("4 times (1,2) over the integers", show(tuple(4 * a for a in vec("12")), 5), "(4,8)")
    same("(1,2) + (4,3) over the integers", show(tuple(a + b for a, b in zip(vec("12"), vec("43"))), 5), "(5,5)")
    same("(1,2) + (4,3) over GF(5)", show(add(vec("12"), vec("43"), 5), 5), "(0,0)")


def check_subspaces() -> None:
    section("Subspaces")
    V = space(2, 3)
    subsets = [[v for i, v in enumerate(V) if mask >> i & 1] for mask in range(2 ** 8)]
    same("subsets of GF(2)^3", len(subsets), 256)
    same("subsets that pass the subspace test", sum(1 for S in subsets if is_subspace(S, 2)), 16)
    closed = [S for S in subsets if S and all(add(u, v, 2) in S for u in S for v in S)]
    same("nonempty subsets closed under addition alone", len(closed), 16)
    require("C5 contains 00000 and is closed under addition", is_subspace(C5, 2))
    odd3 = [v for v in V if weight(v) % 2 == 1]
    same("the odd-weight words of length 3 contain 000", zero(3) in odd3, False)
    same("100 + 010, two words of weight 1", show(add(vec("100"), vec("010"), 2), 2), "110")
    require("over GF(5), c v is v added to itself c times",
            all(scale(c, v, 5) == combine([1] * c, [v] * c, 5, 2) for c in range(5) for v in space(5, 2)))
    mul4 = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
    R4 = range(4)
    field4 = all(mul4[a][b] == mul4[b][a] for a in R4 for b in R4) and \
        all(mul4[a][mul4[b][c]] == mul4[mul4[a][b]][c] for a in R4 for b in R4 for c in R4) and \
        all(mul4[a][b ^ c] == (mul4[a][b] ^ mul4[a][c]) for a in R4 for b in R4 for c in R4) and \
        all(mul4[1][a] == a for a in R4) and all(any(mul4[a][b] == 1 for b in R4) for a in range(1, 4))
    require("GF(4) with XOR and the product table x^2 = x + 1 is a field", field4)
    same("{0, 1} closed under addition in GF(4)", all((a ^ b) in (0, 1) for a in (0, 1) for b in (0, 1)), True)
    same("the scalar 2 of GF(4) times the vector 1", mul4[2][1], 2)
    same("{0, 1} closed under scaling by GF(4)", all(mul4[c][a] in (0, 1) for c in R4 for a in (0, 1)), False)


def check_span() -> None:
    section("Span")
    gens = [vec("1100"), vec("0110"), vec("0011")]
    same("sizes of the spans of the first 0, 1, 2, 3 vectors", as_ints(len(span(gens[:j], 2, 4)) for j in range(4)), "1 2 4 8")
    same("span of 1100", as_set(span(gens[:1], 2, 4), 2), "{0000, 1100}")
    same("span of 1100, 0110", as_set(span(gens[:2], 2, 4), 2), "{0000, 0110, 1010, 1100}")
    same("span of 1100, 0110, 0011", as_set(span(gens, 2, 4), 2),
         "{0000, 0011, 0101, 0110, 1001, 1010, 1100, 1111}")
    require("that span is the set of even-weight words of length 4", span(gens, 2, 4) == even_words(4))
    same("span of (1,2) in GF(5)^2", as_set(span([vec("12")], 5, 2), 5), "{(0,0), (1,2), (2,4), (3,1), (4,3)}")
    same("size of the span of (1,2), (0,1)", len(span([vec("12"), vec("01")], 5, 2)), 25)
    same("size of the span of (1,2), (2,4)", len(span([vec("12"), vec("24")], 5, 2)), 5)
    require("every span in GF(2)^4 of up to 3 vectors is a subspace",
            all(is_subspace(span(vs, 2, 4), 2) for k in range(4) for vs in product(space(2, 4), repeat=k)))
    ok = True
    for q, n in ((2, 4), (3, 2), (5, 2)):
        for level in subspaces(q, n):
            for S in level:
                for w in space(q, n):
                    if w not in S:
                        cosets = {add(scale(c, w, q), s, q) for c in range(q) for s in S}
                        new = span(list(S) + [w], q, n)
                        ok = ok and set(new) == cosets and len(new) == q * len(S)
    require("a vector outside a subspace multiplies its size by q", ok)
    same("y - 2x mod 5 up the column x = 0", as_ints((y - 2 * 0) % 5 for y in range(5)), "0 1 2 3 4")
    same("y - 2x mod 5 up the column x = 1", as_ints((y - 2 * 1) % 5 for y in range(5)), "3 4 0 1 2")
    L = span([vec("12")], 5, 2)
    same("points of GF(5)^2 with digit 0 form the span of (1,2)",
         [v for v in space(5, 2) if (v[1] - 2 * v[0]) % 5 == 0] == L, True)
    same("the digit on the coset (0,1) + span of (1,2)",
         as_ints(sorted({(y - 2 * x) % 5 for (x, y) in [add(vec("01"), v, 5) for v in L]})), "1")


def check_independence() -> None:
    section("Linear independence")
    same("1100, 0110, 0011 independent", independent([vec("1100"), vec("0110"), vec("0011")], 2, 4), True)
    same("1100, 0110, 1010 independent", independent([vec("1100"), vec("0110"), vec("1010")], 2, 4), False)
    same("1100 + 0110", show(add(vec("1100"), vec("0110"), 2), 2), "1010")
    same("1100, 0110, 0011, 1111 independent",
         independent([vec("1100"), vec("0110"), vec("0011"), vec("1111")], 2, 4), False)
    same("1100 + 0011", show(add(vec("1100"), vec("0011"), 2), 2), "1111")
    same("(1,2), (2,4) independent in GF(5)^2", independent([vec("12"), vec("24")], 5, 2), False)
    same("(1,2), (0,1) independent in GF(5)^2", independent([vec("12"), vec("01")], 5, 2), True)
    ok = True
    for k in range(1, 4):
        for vs in product(space(2, 3), repeat=k):
            sums = [combine(cs, vs, 2, 3) for cs in product(range(2), repeat=k) if any(cs)]
            ok = ok and independent(vs, 2, 3) == all(s != zero(3) for s in sums)
    require("over GF(2), independent means no nonempty sublist sums to 0", ok)
    ok = True
    for q, n in ((2, 3), (3, 2)):
        for k in range(n):
            for vs in product(space(q, n), repeat=k):
                if independent(vs, q, n):
                    members = set(span(vs, q, n))
                    ok = ok and all(independent(list(vs) + [v], q, n) for v in space(q, n) if v not in members)
    require("adding a vector outside the span keeps a list independent", ok)
    same("1100 + 0110 + 1010 over GF(2)", show(add(add(vec("1100"), vec("0110"), 2), vec("1010"), 2), 2), "0000")
    same("2 (1,2) + 4 (2,4) over the integers", show(tuple(2 * a + 4 * b for a, b in zip(vec("12"), vec("24"))), 5),
         "(10,20)")
    same("2 (1,2) + 4 (2,4) over GF(5)", show(combine([2, 4], [vec("12"), vec("24")], 5, 2), 5), "(0,0)")


def check_basis() -> None:
    section("Basis and dimension")
    b1, b2 = vec("110"), vec("011")
    same("coordinates of 101 in the basis 110, 011",
         as_ints(cs for cs in product(range(2), repeat=2) if combine(cs, [b1, b2], 2, 3) == vec("101")), "(1, 1)")
    same("coordinates of 011 in the basis 110, 011",
         as_ints(cs for cs in product(range(2), repeat=2) if combine(cs, [b1, b2], 2, 3) == vec("011")), "(0, 1)")
    ok = True
    for q, n, top in ((2, 4, 3), (3, 3, 2), (5, 2, 2)):
        for k in range(top + 1):
            for vs in product(space(q, n), repeat=k):
                if independent(vs, q, n):
                    images = [combine(cs, vs, q, n) for cs in product(range(q), repeat=k)]
                    ok = ok and len(set(images)) == q ** k
    require("k independent vectors give q^k distinct combinations", ok)
    same("dimension of C5", dimension(C5, 2), 2)
    same("dimension of R5", dimension(R5, 2), 1)
    same("dimension of the even-weight words of length 4", dimension(even_words(4), 2), 3)
    S = span([vec("120"), vec("013")], 5, 3)
    same("size of the span of (1,2,0), (0,1,3) in GF(5)^3", len(S), 25)
    same("its dimension", dimension(S, 5), 2)
    sizes = set()
    for P in subspaces(2, 3)[2]:
        for k in range(1, 4):
            for vs in product(P, repeat=k):
                if independent(vs, 2, 3) and span(vs, 2, 3) == list(P):
                    sizes.add(k)
    same("sizes of the bases of the planes of GF(2)^3", as_ints(sorted(sizes)), "2")
    g = [vec("10110"), vec("01101")]
    same("01101 + 10110 over GF(2)", show(add(g[1], g[0], 2), 2), "11011")
    same("the message 11 with the basis 10110, 01101", show(combine([1, 1], g, 2, 5), 2), "11011")
    messages = [combine(cs, g, 2, 5) for cs in product(range(2), repeat=2)]
    same("the 4 messages encode to the words of C5, each once", sorted(messages) == sorted(C5), True)


def check_counting() -> None:
    section("Counting bases and subspaces")
    same("ordered bases of GF(2)^2, by trying every list", lists_independent(2, 2, 2), 6)
    same("ordered bases of GF(2)^3, by trying every list", lists_independent(2, 3, 3), 168)
    same("ordered bases of GF(2)^3 by the formula, 7 x 6 x 4", ordered_independent(2, 3, 3), 168)
    same("ordered independent pairs in GF(2)^3", lists_independent(2, 3, 2), 42)
    same("ordered bases of GF(2)^4, by trying every list", lists_independent(2, 4, 4), 20160)
    same("ordered bases of GF(3)^2, by trying every list", lists_independent(3, 2, 2), 48)
    same("ordered bases of GF(5)^2, by trying every list", lists_independent(5, 2, 2), 480)
    same("ordered bases of GF(3)^3, by trying every list", lists_independent(3, 3, 3), 11232)
    same("ordered bases of GF(3)^3 by the formula, 26 x 24 x 18", ordered_independent(3, 3, 3), 11232)
    same("ordered bases of each plane of GF(2)^3", ordered_independent(2, 2, 2), 6)
    stated = {(2, 2): "1 3 1", (2, 3): "1 7 7 1", (2, 4): "1 15 35 15 1",
              (3, 2): "1 4 1", (3, 3): "1 13 13 1", (5, 2): "1 6 1"}
    totals = {(2, 2): 5, (2, 3): 16, (2, 4): 67, (3, 2): 6, (3, 3): 28, (5, 2): 8}
    for (q, n), text in stated.items():
        levels = subspaces(q, n)
        same(f"subspaces of GF({q})^{n} by dimension, built one by one", as_ints(len(L) for L in levels), text)
        same(f"subspaces of GF({q})^{n} by dimension, by the formula",
             as_ints(gaussian(n, k, q) for k in range(n + 1)), text)
        same(f"subspaces of GF({q})^{n} in all", sum(len(L) for L in levels), totals[(q, n)])
    same("subspaces of GF(2)^7 by dimension", as_ints(gaussian(7, k, 2) for k in range(8)),
         "1 127 2667 11811 11811 2667 127 1")
    same("subspaces of GF(2)^7 in all", sum(gaussian(7, k, 2) for k in range(8)), 29212)
    same("sets of k of the 7 positions, k = 0 to 7", as_ints(comb(7, k) for k in range(8)), "1 7 21 35 35 21 7 1")
    require("the count for k equals the count for n - k, n up to 8, q = 2, 3, 5",
            all(gaussian(n, k, q) == gaussian(n, n - k, q) for q in (2, 3, 5) for n in range(9) for k in range(n + 1)))
    same("fraction of lists of 2 vectors of GF(2)^2 that are bases",
         fraction_text(ordered_independent(2, 2, 2), 2 ** 4), "3/8")
    same("fraction of lists of 3 vectors of GF(2)^3 that are bases",
         fraction_text(ordered_independent(2, 3, 3), 2 ** 9), "21/64")
    same("that fraction from the product, GF(2)^2", f"{basis_fraction(2, 2):.3f}", "0.375")
    same("fraction of lists of 10 vectors of GF(2)^10 that are bases", f"{basis_fraction(2, 10):.3f}", "0.289")
    same("fraction of lists of 10 vectors of GF(3)^10 that are bases", f"{basis_fraction(3, 10):.3f}", "0.560")
    same("fraction of lists of 10 vectors of GF(5)^10 that are bases", f"{basis_fraction(5, 10):.3f}", "0.760")
    same("the factors 2^3 - 2^i for i = 0, 1, 2", as_ints(2 ** 3 - 2 ** i for i in range(3)), "7 6 4")
    for (q, n), text in {(2, 2): "3 2", (3, 2): "8 6", (5, 2): "24 20", (3, 3): "26 24 18"}.items():
        same(f"the factors q^n - q^i for GF({q})^{n}", as_ints(q ** n - q ** i for i in range(n)), text)
    same("bases over lists, GF(2)^2, before reducing", f"{ordered_independent(2, 2, 2)}/{2 ** 4}", "6/16")
    same("bases over lists, GF(2)^3, before reducing", f"{ordered_independent(2, 3, 3)}/{2 ** 9}", "168/512")
    p, tail, ok_low, ok_down, ok_tail = 1.0, 1.0, True, True, True
    for n in range(1, 61):
        new = p * (1 - 2.0 ** (-n))
        ok_down = ok_down and new <= p
        p = new
        ok_low = ok_low and p >= 0.25
        if n >= 2:
            tail = tail * (1 - 2.0 ** (-n))
            ok_tail = ok_tail and tail >= 0.5
    require("over GF(2) the product never increases, n = 1 to 60", ok_down)
    require("the factors after the first keep it at least 1/2, n = 2 to 60", ok_tail)
    require("over GF(2) the product stays at or above 1/4, n = 1 to 60", ok_low)


def check_dot() -> None:
    section("Dot products that are not positive")
    same("11 . 11 over GF(2)", dot(vec("11"), vec("11"), 2), 0)
    same("11 . 11 over the integers", sum(a * b for a, b in zip(vec("11"), vec("11"))), 2)
    same("1111 . 1111 over GF(2)", dot(vec("1111"), vec("1111"), 2), 0)
    same("111 . 111 over GF(2)", dot(vec("111"), vec("111"), 2), 1)
    same("1100 . 1010 over GF(2)", dot(vec("1100"), vec("1010"), 2), 1)
    same("1011 . 1101 over GF(2)", dot(vec("1011"), vec("1101"), 2), 0)
    x = vec("1111")
    ands = tuple(a * b for a, b in zip(x, x))
    same("the four AND outputs for 1111 and 1111", show(ands, 2), "1111")
    same("the two first-level XOR outputs", as_ints([ands[0] ^ ands[1], ands[2] ^ ands[3]]), "0 0")
    same("the final XOR output", ands[0] ^ ands[1] ^ ands[2] ^ ands[3], 0)
    require("over GF(2), v . v = 0 exactly when v has even weight, length 1 to 8",
            all((dot(v, v, 2) == 0) == (weight(v) % 2 == 0) for n in range(1, 9) for v in space(2, n)))
    same("self-orthogonal vectors in GF(2)^4", sum(1 for v in space(2, 4) if dot(v, v, 2) == 0), 8)
    same("(1,1,1) . (1,1,1) over GF(3)", dot(vec("111"), vec("111"), 3), 0)
    same("nonzero self-orthogonal vectors in GF(3)^2",
         sum(1 for v in space(3, 2) if any(v) and dot(v, v, 3) == 0), 0)
    same("(1,2) . (1,2) over GF(5)", dot(vec("12"), vec("12"), 5), 0)
    so5 = [v for v in space(5, 2) if any(v) and dot(v, v, 5) == 0]
    same("nonzero self-orthogonal vectors in GF(5)^2", len(so5), 8)
    same("those vectors", as_set(so5, 5), "{(1,2), (1,3), (2,1), (2,4), (3,1), (3,4), (4,2), (4,3)}")
    ok = True
    for q, n in ((2, 3), (3, 2), (5, 2)):
        V = space(q, n)
        for u in V:
            for v in V:
                ok = ok and dot(u, v, q) == dot(v, u, q)
                for w in V:
                    ok = ok and dot(add(u, v, q), w, q) == (dot(u, w, q) + dot(v, w, q)) % q
                for c in range(q):
                    ok = ok and dot(scale(c, u, q), v, q) == (c * dot(u, v, q)) % q
    require("symmetric and linear in each slot on GF(2)^3, GF(3)^2, GF(5)^2", ok)
    require("the syndrome of C5 is x . 11100, x . 10010, x . 01001 for all 32 words",
            all(syndrome(x) == tuple(dot(x, h, 2) for h in CHECKS5) for x in space(2, 5)))
    same("positions where 1100 and 1010 both have a 1",
         as_ints(i + 1 for i in range(4) if vec("1100")[i] and vec("1010")[i]), "1")
    same("1^2 and 2^2 in GF(3)", as_ints([1 * 1 % 3, 2 * 2 % 3]), "1 1")
    require("over GF(3)^3, v . v is the nonzero-entry count mod 3",
            all(dot(v, v, 3) == weight(v) % 3 for v in space(3, 3)))
    same("(1,2) . (1,2) over the integers", sum(a * b for a, b in zip(vec("12"), vec("12"))), 5)
    same("they fill the lines through (1,2) and (1,3), less the origin",
         sorted((set(span([vec("12")], 5, 2)) | set(span([vec("13")], 5, 2))) - {zero(2)}) == so5, True)
    require("(x, y) . (3, 1) = y - 2x mod 5 at every point of GF(5)^2",
            all(dot(v, vec("31"), 5) == (v[1] - 2 * v[0]) % 5 for v in space(5, 2)))
    same("-2 modulo 5", (-2) % 5, 3)


def check_dual() -> None:
    section("Orthogonal complements and dual codes")
    D = dual(C5, 2, 5)
    same("the dual of C5", as_set(D, 2), "{00000, 00111, 01001, 01110, 10010, 10101, 11011, 11100}")
    same("its size", len(D), 8)
    same("its dimension", dimension(D, 2), 3)
    same("dim C5 + dim of its dual", dimension(C5, 2) + dimension(D, 2), 5)
    same("words in both C5 and its dual", as_set([v for v in C5 if v in D], 2), "{00000, 11011}")
    require("11100, 10010, 01001 lie in the dual of C5", all(h in D for h in CHECKS5))
    same("11100, 10010, 01001 independent", independent(CHECKS5, 2, 5), True)
    same("their span is the dual of C5", span(CHECKS5, 2, 5) == D, True)
    same("the dual of the dual of C5", as_set(dual(D, 2, 5), 2), "{00000, 01101, 10110, 11011}")
    same("C5 on positions 1 and 2", " ".join(show(v[:2], 2) for v in sorted(C5)), "00 01 10 11")
    same("free choices 100, 010, 001 for y3 y4 y5 give",
         as_row([((y[0] + y[1]) % 2, (y[0] + y[2]) % 2) + y for y in (vec("100"), vec("010"), vec("001"))], 2),
         "11100 10010 01001")
    same("the dual of R3", as_set(dual(R3, 2, 3), 2), "{000, 011, 101, 110}")
    same("size of the dual of R5", len(dual(R5, 2, 5)), 16)
    require("the dual of R5 is the set of even-weight words of length 5", dual(R5, 2, 5) == even_words(5))
    C = span([vec("1100"), vec("0111")], 2, 4)
    same("the code spanned by 1100, 0111", as_set(C, 2), "{0000, 0111, 1011, 1100}")
    same("its dual", as_set(dual(C, 2, 4), 2), "{0000, 0011, 1101, 1110}")
    same("words in both", as_set([v for v in C if v in dual(C, 2, 4)], 2), "{0000}")
    L = span([vec("12")], 5, 2)
    same("the dual of the span of (1,2) in GF(5)^2", as_set(dual(L, 5, 2), 5), "{(0,0), (1,2), (2,4), (3,1), (4,3)}")
    same("the coset (0,1) + span of (1,2)", as_set([add(vec("01"), v, 5) for v in L], 5),
         "{(0,1), (1,3), (2,0), (3,2), (4,4)}")
    require("y - 2x is the same on every point of each coset of the line",
            all(len({(y - 2 * x) % 5 for (x, y) in [add(vec("0" + str(c)), v, 5) for v in L]}) == 1 for c in range(5)))
    ok, count = True, 0
    for q, n in ((2, 4), (3, 3), (5, 2)):
        for level in subspaces(q, n):
            for S in level:
                D = dual(S, q, n)
                ok = ok and dimension(S, q) + dimension(D, q) == n and dual(D, q, n) == list(S)
                count += 1
    same("subspaces of GF(2)^4, GF(3)^3 and GF(5)^2 tested", count, 103)
    require("dim C + dim of its dual = n, and the dual of the dual is C", ok)
    require("every subspace of GF(2)^4 has an information set",
            all(has_information_set(S, 4) for level in subspaces(2, 4) for S in level))
    require("y . 10110 = y1 + y3 + y4 and y . 01101 = y2 + y3 + y5, all y",
            all(dot(y, vec("10110"), 2) == (y[0] + y[2] + y[3]) % 2 and dot(y, vec("01101"), 2) == (y[1] + y[2] + y[4]) % 2
                for y in space(2, 5)))
    same("11100, 10010, 01001 on positions 3, 4, 5", " ".join(show(h[2:], 2) for h in CHECKS5), "100 010 001")
    D5 = dual(C5, 2, 5)
    require("every y in the dual is y3 11100 + y4 10010 + y5 01001", all(combine(y[2:], CHECKS5, 2, 5) == y for y in D5))
    require("in C5, x3 = x1 + x2, x4 = x1 and x5 = x2",
            all(c[2] == (c[0] + c[1]) % 2 and c[3] == c[0] and c[4] == c[1] for c in C5))


def check_self_dual() -> None:
    section("Self-dual codes")
    two = [vec("00"), vec("11")]
    same("{00, 11} equals its dual", dual(two, 2, 2) == two, True)
    same("binary self-dual codes of length 2, 4, 6", as_ints(self_dual_count(n) for n in (2, 4, 6)), "1 3 15")
    found = [S for S in subspaces(2, 4)[2] if dual(S, 2, 4) == list(S)]
    same("the self-dual codes of length 4", " ".join(as_set(S, 2) for S in found),
         "{0000, 0011, 1100, 1111} {0000, 0101, 1010, 1111} {0000, 0110, 1001, 1111}")
    require("1111 lies in each of them", all(vec("1111") in S for S in found))
    G = [vec("11110000"), vec("00111100"), vec("00001111"), vec("01010101")]
    same("the ten dot products among the four rows", as_ints(dot(G[i], G[j], 2) for i in range(4) for j in range(i, 4)),
         "0 0 0 0 0 0 0 0 0 0")
    same("the four rows independent", independent(G, 2, 8), True)
    H = span(G, 2, 8)
    same("words in the length-8 code", len(H), 16)
    same("weights of its words, weight:count", weight_counts(H), "0:1 4:14 8:1")
    same("it equals its dual", dual(H, 2, 8) == H, True)
    same("it contains 11111111", vec("11111111") in H, True)
    T = [vec("1110"), vec("0121")]
    same("the three dot products among the two rows over GF(3)", as_ints([dot(T[0], T[0], 3), dot(T[0], T[1], 3), dot(T[1], T[1], 3)]), "0 0 0")
    tetra = span(T, 3, 4)
    same("words in the ternary code of length 4", len(tetra), 9)
    same("weights of its words, weight:count", weight_counts(tetra), "0:1 3:8")
    same("it equals its dual", dual(tetra, 3, 4) == tetra, True)
    lines5 = [S for S in subspaces(5, 2)[1] if dual(S, 5, 2) == list(S)]
    same("self-dual lines of GF(5)^2", len(lines5), 2)
    same("their nonzero points with first entry 1", as_row([v for S in lines5 for v in S if v[0] == 1], 5), "(1,2) (1,3)")
    same("self-dual lines of GF(3)^2", sum(1 for S in subspaces(3, 2)[1] if dual(S, 3, 2) == list(S)), 0)
    same("prod of 2^i + 1 for i = 1 to n/2 - 1, n = 2, 4, 6",
         as_ints(prod(2 ** i + 1 for i in range(1, n // 2)) for n in (2, 4, 6)), "1 3 15")
    E = even_words(3)
    require("every word of E has v . v = 0", all(dot(v, v, 2) == 0 for v in E))
    same("011 . 101 over GF(2)", dot(vec("011"), vec("101"), 2), 1)
    same("the ternary dot products over the integers",
         as_ints(sum(a * b for a, b in zip(T[i], T[j])) for i, j in ((0, 0), (0, 1), (1, 1))), "3 3 6")
    same("1110, 0121 independent over GF(3)", independent(T, 3, 4), True)


def check_worked_examples() -> None:
    section("Worked examples")
    T = [v for v in space(3, 3) if sum(v) % 3 == 0]
    same("vectors of GF(3)^3 whose entries sum to 0", len(T), 9)
    same("they pass the subspace test", is_subspace(T, 3), True)
    same("(1,2,0), (0,1,2) independent", independent([vec("120"), vec("012")], 3, 3), True)
    same("their span is that set", span([vec("120"), vec("012")], 3, 3) == T, True)
    same("that set is the dual of the span of (1,1,1)", dual(span([vec("111")], 3, 3), 3, 3) == T, True)
    same("(1,1,1) lies in it", vec("111") in T, True)
    v1, v2, v3 = vec("120"), vec("013"), vec("203")
    same("2 (1,2,0) + (0,1,3) over GF(5)", show(add(scale(2, v1, 5), v2, 5), 5), "(2,0,3)")
    same("(1,2,0), (0,1,3), (2,0,3) independent", independent([v1, v2, v3], 5, 3), False)
    same("size of their span", len(span([v1, v2, v3], 5, 3)), 25)
    same("lines of GF(3)^3", gaussian(3, 1, 3), 13)
    same("planes of GF(3)^3", gaussian(3, 2, 3), 13)
    same("lines inside each plane of GF(3)^3", gaussian(2, 1, 3), 4)
    planes = subspaces(3, 3)[2]
    line = span([vec("111")], 3, 3)
    same("planes of GF(3)^3 that contain the line through (1,1,1)",
         sum(1 for P in planes if set(line) <= set(P)), 4)
    same("those vectors", as_set(T, 3),
         "{(0,0,0), (0,1,2), (0,2,1), (1,0,2), (1,1,1), (1,2,0), (2,0,1), (2,1,0), (2,2,2)}")
    same("their dimension", dimension(T, 3), 2)
    same("2 (1,2,0) over GF(5)", show(scale(2, v1, 5), 5), "(2,4,0)")
    same("(2,4,0) + (0,1,3) over the integers", show(tuple(a + b for a, b in zip(scale(2, v1, 5), v2)), 5), "(2,5,3)")
    same("4 (2,0,3) over the integers", show(tuple(4 * a for a in v3), 5), "(8,0,12)")
    same("4 (2,0,3) over GF(5)", show(scale(4, v3, 5), 5), "(3,0,2)")
    same("2 v1 + v2 + 4 v3 over the integers",
         show(tuple(a + b + c for a, b, c in zip(scale(2, v1, 5), v2, scale(4, v3, 5))), 5), "(5,5,5)")
    same("2 v1 + v2 + 4 v3 over GF(5)", show(combine([2, 1, 4], [v1, v2, v3], 5, 3), 5), "(0,0,0)")
    same("(1,2,0), (0,1,3) independent", independent([v1, v2], 5, 3), True)
    same("coordinates of (2,0,3) in the basis (1,2,0), (0,1,3)",
         as_ints(cs for cs in product(range(5), repeat=2) if combine(cs, [v1, v2], 5, 3) == v3), "(2, 1)")
    same("dimension of the span of all three", dimension(span([v1, v2, v3], 5, 3), 5), 2)
    C4 = [vec("0000"), vec("1100"), vec("0111"), vec("1011")]
    same("0000, 1100, 0111, 1011 is the span of 1100, 0111", sorted(C4) == span([vec("1100"), vec("0111")], 2, 4), True)
    same("positions 1, 2 of 0000, 1100, 0111, 1011", " ".join(show(c[:2], 2) for c in C4), "00 11 01 10")
    words = [((a + b) % 2, (a + b) % 2, a, b) for a, b in (vec("00"), vec("10"), vec("01"), vec("11"))]
    same("free choices 00, 10, 01, 11 for y3 y4 give", as_row(words, 2), "0000 1110 1101 0011")
    same("those words are the dual of that code", sorted(words) == dual(C4, 2, 4), True)
    same("0111 . 0111 over GF(2)", dot(vec("0111"), vec("0111"), 2), 1)
    r = [vec("11110000"), vec("00111100"), vec("00001111"), vec("01010101")]
    pairs = ((0, 1), (0, 3), (1, 2), (1, 3), (2, 3), (0, 2))
    shared = ["".join(str(p + 1) for p in range(8) if r[i][p] and r[j][p]) or "-" for i, j in pairs]
    same("positions shared by r1r2 r1r4 r2r3 r2r4 r3r4 r1r3", " ".join(shared), "34 24 56 46 68 -")
    same("rows holding a 1 at positions 1, 7, 8 and 3",
         " ".join("".join(str(i + 1) for i in range(4) if r[i][p - 1]) for p in (1, 7, 8, 3)), "1 3 34 12")
    same("r1 + r3", show(add(r[0], r[2], 2), 2), "11111111")
    H = span(r, 2, 8)

    def common(x: tuple, y: tuple) -> int:
        return sum(1 for a, b in zip(x, y) if a and b)

    require("wt(x + y) = wt(x) + wt(y) - 2m on every pair of the code",
            all(weight(add(x, y, 2)) == weight(x) + weight(y) - 2 * common(x, y) for x in H for y in H))
    require("m is even on every pair of the code", all(common(x, y) % 2 == 0 for x in H for y in H))
    same("its words of weight 0 and 8", as_row([v for v in H if weight(v) in (0, 8)], 2), "00000000 11111111")


def check_exercises() -> None:
    section("Exercises")
    lines = subspaces(3, 2)[1]
    gens = sorted(next(v for v in S if any(v) and next(a for a in v if a) == 1) for S in lines)
    same("lines of GF(3)^2, each named by a point with leading entry 1", as_row(gens, 3), "(0,1) (1,0) (1,1) (1,2)")
    same("1101 + 0110", show(add(vec("1101"), vec("0110"), 2), 2), "1011")
    same("1101, 0110, 1011 independent", independent([vec("1101"), vec("0110"), vec("1011")], 2, 4), False)
    same("2-dimensional subspaces of GF(2)^4, 15 x 14 / (3 x 2)", gaussian(4, 2, 2), 35)
    B = [vec("11000"), vec("01100"), vec("00110"), vec("00011")]
    same("11000, 01100, 00110, 00011 independent", independent(B, 2, 5), True)
    same("their span is the dual of R5", span(B, 2, 5) == dual(R5, 2, 5), True)
    same("the squares in GF(3)", as_ints(sorted({a * a % 3 for a in range(3)})), "0 1")
    same("the solutions of b^2 = 4 in GF(5)", as_ints(b for b in range(5) if b * b % 5 == 4), "2 3")
    same("(1,3) . (1,3) over GF(5)", dot(vec("13"), vec("13"), 5), 0)
    same("points on each line of GF(3)^2", as_ints(len(S) for S in lines), "3 3 3 3")
    same("nonzero points of GF(3)^2", sum(1 for v in space(3, 2) if any(v)), 8)
    require("each nonzero point of GF(3)^2 lies on exactly one line",
            all(sum(1 for S in lines if v in S) == 1 for v in space(3, 2) if any(v)))
    same("the factors q^n - q^i for GF(2)^4", as_ints(2 ** 4 - 2 ** i for i in range(4)), "15 14 12 8")
    same("weights of 11000, 01100, 00110, 00011", as_ints(weight(v) for v in B), "2 2 2 2")
    same("0^2, 1^2, 2^2 in GF(3)", as_ints(a * a % 3 for a in range(3)), "0 1 1")
    w2 = [v for v in space(2, 4) if weight(v) == 2]
    same("words of weight 2 in GF(2)^4", len(w2), 6)
    flips = [as_set([v, add(v, vec("1111"), 2)], 2) for v in w2 if v < add(v, vec("1111"), 2)]
    same("each paired with its flip", " ".join(flips), "{0011, 1100} {0101, 1010} {0110, 1001}")
    same("(0,1) . (0,1) over GF(5)", dot(vec("01"), vec("01"), 5), 1)
    same("(1,3) . (1,3) over the integers", sum(a * b for a, b in zip(vec("13"), vec("13"))), 10)


def main() -> None:
    check_scaling()
    check_fields()
    check_axioms()
    check_subspaces()
    check_span()
    check_independence()
    check_basis()
    check_counting()
    check_dot()
    check_dual()
    check_self_dual()
    check_worked_examples()
    check_exercises()
    print("\nAll chapter 15 numbers reproduced.")


if __name__ == "__main__":
    main()
