"""Reproduce every number in chapter 32 of Error Correction Coding, "Rings and Fields".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch32_rings_and_fields.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees. All
the arithmetic is on whole numbers, so every comparison is exact.

Every finite ring the chapter uses is checked against the ring axioms by
enumerating every pair and every triple of its elements: the integers modulo
4, 12, 13 and 256, the words of 2 and 8 bits under XOR and AND, the 2 x 2
matrices over GF(2) and the four-element field whose tables the chapter
prints. Beside the worked numbers, each theorem the chapter proves is also
tested on small cases: the first consequences of the axioms, that reduction
modulo n preserves sums and products, that the units form a group, that every
nonzero element of these finite rings is a unit or a zero divisor, the gcd test
for zero divisors, the degree of a product, that the integers modulo n form a
field exactly when n is prime, the characteristic, the divisibility of the
binomial coefficients, the identity (a + b)^p = a^p + b^p, squaring in GF(2)[x]
as a spreading of bits, and which sizes up to 64 are prime powers. A run that
reaches the final line has reproduced every number in the chapter.

The module also holds the arithmetic that ch32_charts.py draws, so the charts
and this check share one source.
"""

import math
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


def yes(flag: bool) -> str:
    return "yes" if flag else "no"


# ---------------------------------------------------------------------------
# Whole numbers.

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def is_prime_power(n: int) -> bool:
    if n < 2:
        return False
    p = 2
    while n % p != 0:
        p += 1
    while n % p == 0:
        n //= p
    return n == 1


def phi(n: int) -> int:
    return sum(1 for a in range(1, n + 1) if math.gcd(a, n) == 1)


def pascal_row(n: int) -> list[int]:
    row = [1]
    for _ in range(n):
        row = [1] + [row[i] + row[i + 1] for i in range(len(row) - 1)] + [1]
    return row


# ---------------------------------------------------------------------------
# The finite rings of the chapter. Each is a dictionary holding its addition
# and multiplication tables on the elements 0 .. n-1, with 0 the zero element.

F4_MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
F4_NAMES = "01ab"


def ring_zn(n: int) -> dict:
    return {"name": f"Z_{n}", "kind": "zn", "n": n, "k": 0, "one": 1,
            "add": [[(a + b) % n for b in range(n)] for a in range(n)],
            "mul": [[a * b % n for b in range(n)] for a in range(n)]}


def ring_words(k: int) -> dict:
    n = 1 << k
    return {"name": f"{k}-bit words", "kind": "words", "n": n, "k": k, "one": n - 1,
            "add": [[a ^ b for b in range(n)] for a in range(n)],
            "mul": [[a & b for b in range(n)] for a in range(n)]}


def mat_mul(x: int, y: int) -> int:
    """Multiply two 2 x 2 matrices over GF(2), each packed as the bits a b c d of [[a, b], [c, d]]."""
    a, b, c, d = (x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1
    e, f, g, h = (y >> 3) & 1, (y >> 2) & 1, (y >> 1) & 1, y & 1
    return (((a * e + b * g) % 2) << 3 | ((a * f + b * h) % 2) << 2
            | ((c * e + d * g) % 2) << 1 | ((c * f + d * h) % 2))


def ring_matrices() -> dict:
    return {"name": "M2(GF(2))", "kind": "mat", "n": 16, "k": 0, "one": 0b1001,
            "add": [[a ^ b for b in range(16)] for a in range(16)],
            "mul": [[mat_mul(a, b) for b in range(16)] for a in range(16)]}


def ring_f4() -> dict:
    return {"name": "F_4", "kind": "f4", "n": 4, "k": 0, "one": 1,
            "add": [[a ^ b for b in range(4)] for a in range(4)],
            "mul": [row[:] for row in F4_MUL]}


def show(ring: dict, a: int) -> str:
    kind = ring["kind"]
    if kind == "words":
        return format(a, f"0{ring['k']}b")
    if kind == "mat":
        return f"{(a >> 3) & 1}{(a >> 2) & 1}|{(a >> 1) & 1}{a & 1}"
    if kind == "f4":
        return F4_NAMES[a]
    return str(a)


def show_all(ring: dict, values: list[int]) -> str:
    return " ".join(show(ring, a) for a in values) if values else "none"


def chapter_rings() -> list[dict]:
    return [ring_zn(4), ring_zn(12), ring_zn(13), ring_zn(256),
            ring_words(2), ring_words(8), ring_matrices(), ring_f4()]


def _tbl(row: list[int]) -> bytes:
    """A row as a 256-entry translation table, so bytes.translate looks it up for a whole row at once."""
    return bytes(row) + bytes(256 - len(row))


def associative(t: list[list[int]]) -> bool:
    n = len(t)
    rows = [bytes(r) for r in t]
    tables = [_tbl(r) for r in t]
    for a in range(n):
        ta = t[a]
        for b in range(n):
            if rows[ta[b]] != rows[b].translate(tables[a]):   # (ab)c against a(bc), every c
                return False
    return True


def distributive(add: list[list[int]], mul: list[list[int]]) -> bool:
    n = len(add)
    add_rows = [bytes(r) for r in add]
    add_tables = [_tbl(r) for r in add]
    for a in range(n):
        ma = mul[a]
        col = [mul[x][a] for x in range(n)]
        ma_tbl, col_tbl = _tbl(ma), _tbl(col)
        ma_row, col_row = bytes(ma), bytes(col)
        for b in range(n):
            if add_rows[b].translate(ma_tbl) != ma_row.translate(add_tables[ma[b]]):     # a(b+c) = ab + ac
                return False
            if add_rows[b].translate(col_tbl) != col_row.translate(add_tables[col[b]]):  # (b+c)a = ba + ca
                return False
    return True


def commutative(t: list[list[int]]) -> bool:
    n = len(t)
    return all(t[a][b] == t[b][a] for a in range(n) for b in range(n))


def units(ring: dict) -> list[int]:
    n, mul, one = ring["n"], ring["mul"], ring["one"]
    return [u for u in range(n) if any(mul[u][v] == one and mul[v][u] == one for v in range(n))]


def inverse_in(ring: dict, u: int) -> int:
    mul, one = ring["mul"], ring["one"]
    return next(v for v in range(ring["n"]) if mul[u][v] == one and mul[v][u] == one)


def zero_divisors(ring: dict) -> list[int]:
    n, mul = ring["n"], ring["mul"]
    return [a for a in range(1, n) if any(mul[a][b] == 0 or mul[b][a] == 0 for b in range(1, n))]


def characteristic(ring: dict) -> int:
    add, one = ring["add"], ring["one"]
    s, k = one, 1
    while s != 0:
        s = add[s][one]
        k += 1
    return k


def multiples_of_one(ring: dict) -> list[int]:
    add, one = ring["add"], ring["one"]
    seen, s = [0], one
    while s != 0:
        seen.append(s)
        s = add[s][one]
    return sorted(seen)


def zn_units(n: int) -> list[int]:
    return [a for a in range(1, n) if any(a * b % n == 1 for b in range(n))]


def zn_zero_divisors(n: int) -> list[int]:
    return [a for a in range(1, n) if any(a * b % n == 0 for b in range(1, n))]


def unit_and_zero_divisor_counts(n: int) -> tuple[int, int]:
    """Number of units and of zero divisors of Z_n, the data of the chart."""
    return len(zn_units(n)), len(zn_zero_divisors(n))


def solutions(a: int, c: int, n: int) -> list[int]:
    return [x for x in range(n) if a * x % n == c % n]


def times_map(a: int, n: int) -> list[int]:
    return [a * x % n for x in range(n)]


# ---------------------------------------------------------------------------
# Polynomials. Integer coefficient lists are written constant term first.
# A polynomial over GF(2) is a whole number whose bit i is the coefficient of x^i.

def polymul(f: list[int], g: list[int], modulus: int = 0) -> list[int]:
    out = [0] * (len(f) + len(g) - 1)
    for i, fi in enumerate(f):
        for j, gj in enumerate(g):
            out[i + j] += fi * gj
    return [c % modulus for c in out] if modulus else out


def clmul(a: int, b: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def spread(a: int) -> int:
    r, i = 0, 0
    while a:
        if a & 1:
            r |= 1 << (2 * i)
        a >>= 1
        i += 1
    return r


def degree(a: int) -> int:
    return a.bit_length() - 1


def bits(a: int) -> str:
    return format(a, "b")


# ---------------------------------------------------------------------------
# The checks, in chapter order.

def check_opening() -> None:
    section("The opening: two rings inside an 8-bit processor")
    same("3 * 171", 3 * 171, 513)
    same("2 * 256 + 1", 2 * 256 + 1, 513)
    same("3 * 171 mod 256", 3 * 171 % 256, 1)
    require("multiplying by 171 undoes multiplying by 3, all 256 bytes",
            all(3 * x * 171 % 256 == x for x in range(256)))
    same("2 * 128 mod 256", 2 * 128 % 256, 0)
    same("2 * 0 mod 256", 2 * 0 % 256, 0)
    same("different values of 3x mod 256", len({3 * x % 256 for x in range(256)}), 256)
    same("different values of 2x mod 256", len({2 * x % 256 for x in range(256)}), 128)
    same("0x0F AND 0xF0", 0x0F & 0xF0, 0)
    same("0x0F XOR 0xF0", 0x0F ^ 0xF0, 255)


def check_axioms(rings: list[dict]) -> None:
    section("The ring axioms, checked on every pair and every triple")
    for r in rings:
        n, add, mul, one = r["n"], r["add"], r["mul"], r["one"]
        name = r["name"]
        require(f"{name}: + associative on all {n ** 3} triples", associative(add))
        require(f"{name}: + commutative, 0 an identity, negatives exist",
                commutative(add) and all(add[0][a] == a for a in range(n))
                and all(0 in add[a] for a in range(n)))
        require(f"{name}: x associative on all {n ** 3} triples", associative(mul))
        require(f"{name}: both distributive laws on all {n ** 3} triples", distributive(add, mul))
        require(f"{name}: {show(r, one)} is a two-sided identity for x",
                all(mul[one][a] == a and mul[a][one] == a for a in range(n)))
        same(f"{name}: number of elements, commutative", f"{n} {yes(commutative(mul))}",
             {"Z_4": "4 yes", "Z_12": "12 yes", "Z_13": "13 yes", "Z_256": "256 yes",
              "2-bit words": "4 yes", "8-bit words": "256 yes", "M2(GF(2))": "16 no",
              "F_4": "4 yes"}[name])


def check_gallery() -> None:
    section("A gallery of rings")
    same("(x^2 + x + 1)(x + 1) in Z[x], constant first", polymul([1, 1, 1], [1, 1]), [1, 2, 2, 1])
    same("(x^2 + x + 1)(x + 1) in GF(2)[x], bits", bits(clmul(0b111, 0b11)), "1001")
    words8 = ring_words(8)
    require("a AND a = a for every 8-bit word a", all(words8["mul"][a][a] == a for a in range(256)))
    require("a XOR a = 0 for every 8-bit word a", all(words8["add"][a][a] == 0 for a in range(256)))
    A, B = 0b1101, 0b1011
    m = ring_matrices()
    same("A", show(m, A), "11|01")
    same("B", show(m, B), "10|11")
    same("AB over GF(2)", show(m, mat_mul(A, B)), "01|11")
    same("BA over GF(2)", show(m, mat_mul(B, A)), "11|10")
    same("number of 2 x 2 matrices over GF(2)", 2 ** 4, 16)


def check_consequences(rings: list[dict]) -> None:
    section("First consequences of the axioms")
    ok_zero = ok_neg = ok_ids = ok_unique_neg = True
    for r in rings:
        n, add, mul, one = r["n"], r["add"], r["mul"], r["one"]
        neg = [add[a].index(0) for a in range(n)]
        ok_zero &= all(mul[0][a] == 0 and mul[a][0] == 0 for a in range(n))
        ok_neg &= all(mul[neg[a]][b] == neg[mul[a][b]] and mul[a][neg[b]] == neg[mul[a][b]]
                      for a in range(n) for b in range(n))
        zeros = [z for z in range(n) if all(add[z][a] == a for a in range(n))]
        ones = [e for e in range(n) if all(mul[e][a] == a and mul[a][e] == a for a in range(n))]
        ok_ids &= zeros == [0] and ones == [one]
        ok_unique_neg &= all(add[a].count(0) == 1 for a in range(n))
    require("0a = a0 = 0 for every a, all eight rings", ok_zero)
    require("(-a)b = a(-b) = -(ab) for every a and b, all eight rings", ok_neg)
    require("exactly one zero and one multiplicative identity, all eight", ok_ids)
    require("every element has exactly one negative, all eight rings", ok_unique_neg)


def check_homomorphism() -> None:
    section("Reduction modulo n is a ring homomorphism")
    require("(a + b) mod n = (a mod n + b mod n) mod n, n = 2..12, |a|, |b| <= 30",
            all((a + b) % n == (a % n + b % n) % n
                for n in range(2, 13) for a in range(-30, 31) for b in range(-30, 31)))
    require("(ab) mod n = (a mod n)(b mod n) mod n, n = 2..12, |a|, |b| <= 30",
            all(a * b % n == (a % n) * (b % n) % n
                for n in range(2, 13) for a in range(-30, 31) for b in range(-30, 31)))
    same("47 * 86", 47 * 86, 4042)
    same("47 mod 9, 86 mod 9", [47 % 9, 86 % 9], [2, 5])
    same("(2 * 5) mod 9", 2 * 5 % 9, 1)
    same("4042 mod 9", 4042 % 9, 1)
    same("4052 mod 9, a wrong product", 4052 % 9, 2)


def check_units(rings: list[dict]) -> None:
    section("Units and zero divisors")
    same("units of Z among -20..20", [a for a in range(-20, 21) if a != 0 and 1 % a == 0], [-1, 1])
    group_ok = disjoint_ok = cover_ok = True
    for r in rings:
        u, z = units(r), zero_divisors(r)
        mul = r["mul"]
        group_ok &= all(mul[a][b] in u for a in u for b in u) and all(inverse_in(r, a) in u for a in u)
        disjoint_ok &= not set(u) & set(z)
        cover_ok &= sorted(u + z) == list(range(1, r["n"]))
    require("units closed under x and under inverses, all eight rings", group_ok)
    require("no unit is a zero divisor, all eight rings", disjoint_ok)
    require("every nonzero element is a unit or a zero divisor, all eight", cover_ok)

    z12 = ring_zn(12)
    same("units of Z_12", units(z12), [1, 5, 7, 11])
    same("zero divisors of Z_12", zero_divisors(z12), [2, 3, 4, 6, 8, 9, 10])
    same("squares of the units of Z_12", [u * u % 12 for u in units(z12)], [1, 1, 1, 1])
    same("4 * 3 mod 12, 6 * 2 mod 12", [4 * 3 % 12, 6 * 2 % 12], [0, 0])

    z256 = ring_zn(256)
    u256, zd256 = units(z256), zero_divisors(z256)
    same("units of Z_256, all odd", f"{len(u256)} {yes(all(a % 2 == 1 for a in u256))}", "128 yes")
    same("zero divisors of Z_256, all even", f"{len(zd256)} {yes(all(a % 2 == 0 for a in zd256))}", "127 yes")
    same("inverse of 3 in Z_256", inverse_in(z256, 3), 171)

    w2, w8 = ring_words(2), ring_words(8)
    same("units of the 2-bit words", show_all(w2, units(w2)), "11")
    same("zero divisors of the 2-bit words", show_all(w2, zero_divisors(w2)), "01 10")
    same("units of the 8-bit words", show_all(w8, units(w8)), "11111111")
    same("zero divisors of the 8-bit words", len(zero_divisors(w8)), 254)
    require("a AND (NOT a) = 0 for every 8-bit word a", all(a & (255 ^ a) == 0 for a in range(256)))

    m = ring_matrices()
    um = units(m)
    same("invertible 2 x 2 matrices over GF(2)", len(um), 6)
    same("(2^2 - 1)(2^2 - 2), the count from the matrix chapter", (4 - 1) * (4 - 2), 6)
    same("nonzero matrices that are not invertible", len(zero_divisors(m)), 9)
    require("each of the nine has a nonzero partner on each side",
            all(any(mat_mul(a, b) == 0 for b in range(1, 16)) and any(mat_mul(b, a) == 0 for b in range(1, 16))
                for a in zero_divisors(m)))

    gcd_ok = phi_ok = True
    for n in range(2, 41):
        zu, zz = zn_units(n), zn_zero_divisors(n)
        gcd_ok &= zz == [a for a in range(1, n) if math.gcd(a, n) > 1]
        phi_ok &= zu == [a for a in range(1, n) if math.gcd(a, n) == 1] and len(zu) == phi(n)
        gcd_ok &= len(zu) + len(zz) == n - 1
    require("units of Z_n are the residues coprime to n, phi(n), n = 2..40", phi_ok)
    require("zero divisors of Z_n: nonzero a with gcd(a, n) > 1, n = 2..40", gcd_ok)
    same("x -> 3x on Z_7", times_map(3, 7), [0, 3, 6, 2, 5, 1, 4])
    same("x -> 2x on Z_6", times_map(2, 6), [0, 2, 4, 0, 2, 4])


def check_solving() -> None:
    section("Solving ax = c")
    same("4x = 8 in Z_12", solutions(4, 8, 12), [2, 5, 8, 11])
    same("4x = 6 in Z_12", solutions(4, 6, 12), [])
    same("c in Z_12 for which 4x = c has a solution", [c for c in range(12) if solutions(4, c, 12)], [0, 4, 8])
    same("gcd(4, 12)", math.gcd(4, 12), 4)
    counts = [len(solutions(a, c, 13)) for a in range(1, 13) for c in range(13)]
    same("equations ax = c in Z_13 with a nonzero", len(counts), 156)
    require("each of them has exactly one solution", all(k == 1 for k in counts))
    same("4x = 6 in Z_13", solutions(4, 6, 13), [8])
    same("inverse of 4 in Z_13", zn_units(13) and next(b for b in range(13) if 4 * b % 13 == 1), 10)


def check_domains(rings: list[dict]) -> None:
    section("Integral domains")
    pairs = [(f, g) for f in range(1, 64) for g in range(1, 64)]
    same("pairs of nonzero GF(2)[x] polynomials of degree below 6", len(pairs), 3969)
    require("deg(fg) = deg f + deg g for each of them",
            all(degree(clmul(f, g)) == degree(f) + degree(g) for f, g in pairs))
    same("pairs among them with fg = 1", [[f, g] for f, g in pairs if clmul(f, g) == 1], [[1, 1]])
    same("(2x)(2x) in Z_4[x], constant first", polymul([0, 2], [0, 2], 4), [0, 0, 0])
    same("(2x + 1)(2x + 1) in Z_4[x], constant first", polymul([1, 2], [1, 2], 4), [1, 0, 0])
    ok = True
    for r in rings:
        n, mul = r["n"], r["mul"]
        cancel = all(len(set(mul[a])) == n for a in range(1, n))
        ok &= cancel == (not zero_divisors(r))
    require("cancellation holds exactly when there are no zero divisors", ok)
    same("2 * 3 mod 6", 2 * 3 % 6, 0)


def check_fields(rings: list[dict]) -> None:
    section("Fields")
    require("Z_n is an integral domain exactly when n is prime, n = 2..64",
            all((not zn_zero_divisors(n)) == is_prime(n) for n in range(2, 65)))
    require("Z_n is a field exactly when n is prime, n = 2..64",
            all((len(zn_units(n)) == n - 1) == is_prime(n) for n in range(2, 65)))
    stated = {"Z_4": "no no", "Z_12": "no no", "Z_13": "yes yes", "Z_256": "no no",
              "2-bit words": "no no", "8-bit words": "no no", "M2(GF(2))": "no no", "F_4": "yes yes"}
    finite_domain_ok = True
    for r in rings:
        comm = commutative(r["mul"])
        u, z = units(r), zero_divisors(r)
        domain = comm and not z
        field = comm and len(u) == r["n"] - 1
        finite_domain_ok &= (not domain) or field
        same(f"{r['name']}: integral domain, field", f"{yes(domain)} {yes(field)}", stated[r["name"]])
    require("every finite integral domain among them is a field", finite_domain_ok)
    z4, w2, f4 = ring_zn(4), ring_words(2), ring_f4()
    same("Z_4: 2 * 2", z4["mul"][2][2], 0)
    same("Z_4: units, zero divisors", f"{show_all(z4, units(z4))}; {show_all(z4, zero_divisors(z4))}", "1 3; 2")
    same("2-bit words: 01 AND 10", show(w2, w2["mul"][1][2]), "00")
    same("F_4: a*a, a*b, b*b", " ".join(F4_NAMES[F4_MUL[x][y]] for x, y in [(2, 2), (2, 3), (3, 3)]), "b 1 a")
    same("F_4: 1+a, a+b, 1+b", " ".join(F4_NAMES[x ^ y] for x, y in [(1, 2), (2, 3), (1, 3)]), "b 1 a")
    same("F_4: inverses of 1, a, b", show_all(f4, [inverse_in(f4, u) for u in (1, 2, 3)]), "1 b a")
    same("F_4: units", show_all(f4, units(f4)), "1 a b")


def check_characteristic(rings: list[dict]) -> None:
    section("Characteristic")
    stated = {"Z_4": 4, "Z_12": 12, "Z_13": 13, "Z_256": 256, "2-bit words": 2,
              "8-bit words": 2, "M2(GF(2))": 2, "F_4": 2}
    all_ok = True
    for r in rings:
        c = characteristic(r)
        same(f"characteristic of {r['name']}", c, stated[r["name"]])
        n, add = r["n"], r["add"]
        s = [0] * n
        for _ in range(c):
            s = [add[s[a]][a] for a in range(n)]
        all_ok &= all(v == 0 for v in s)
    require("char * a = 0 for every element a, all eight rings", all_ok)
    require("f + f = 0, so -f = f, every f in GF(2)[x] of degree below 6",
            all(f ^ f == 0 for f in range(64)))
    char_zn = {n: characteristic(ring_zn(n)) for n in range(2, 65)}
    require("characteristic of Z_n is n, n = 2..64", all(c == n for n, c in char_zn.items()))
    require("characteristic of Z_n prime whenever Z_n is a domain, n = 2..64",
            all(is_prime(c) for n, c in char_zn.items() if not zn_zero_divisors(n)))
    same("7 + 1 as integers, bits", format(7 + 1, "04b"), "1000")
    same("0111 XOR 0001", format(0b0111 ^ 0b0001, "04b"), "0110")
    carries, c = [], 0
    for i in range(4):
        x, y = (0b0111 >> i) & 1, (0b0001 >> i) & 1
        c = (x & y) | (x & c) | (y & c)
        carries.append(c)
    same("carries out of bits 0..3 for 0111 + 0001", carries, [1, 1, 1, 0])


def check_freshman() -> None:
    section("The binomial theorem in characteristic p")
    primes = [p for p in range(2, 50) if is_prime(p)]
    same("primes below 50", primes, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47])
    require("p divides C(p, i) for 0 < i < p, every prime p below 50",
            all(all(pascal_row(p)[i] % p == 0 for i in range(1, p)) for p in primes))
    require("n fails to divide some C(n, i), every composite n from 4 to 49",
            all(any(pascal_row(n)[i] % n != 0 for i in range(1, n)) for n in range(4, 50) if not is_prime(n)))
    same("C(5, i) for i = 0..5", pascal_row(5), [1, 5, 10, 10, 5, 1])
    same("C(4, 2)", pascal_row(4)[2], 6)
    require("(a + b)^p = a^p + b^p in Z_p for all a, b, every prime p below 50",
            all(pow(a + b, p, p) == (pow(a, p, p) + pow(b, p, p)) % p
                for p in primes for a in range(p) for b in range(p)))
    same("(1 + 1)^2 and 1^2 + 1^2 in Z_4", [(1 + 1) ** 2 % 4, (1 + 1) % 4], [0, 2])
    require("(f + g)^2 = f^2 + g^2 in GF(2)[x], all f, g of degree below 6",
            all(clmul(f ^ g, f ^ g) == clmul(f, f) ^ clmul(g, g) for f in range(64) for g in range(64)))
    require("squaring spreads bit i to bit 2i, every f of degree below 8",
            all(clmul(f, f) == spread(f) for f in range(256)))
    same("(x^3 + x + 1)^2 in GF(2)[x], bits", f"{bits(0b1011)} -> {bits(clmul(0b1011, 0b1011))}", "1011 -> 1000101")
    same("(x^3 + x + 1)^2 in Z[x], constant first", polymul([1, 1, 0, 1], [1, 1, 0, 1]), [1, 2, 1, 2, 2, 0, 1])
    f4 = ring_f4()
    require("(u + v)^2 = u^2 + v^2 in F_4 for all 16 pairs",
            all(f4["mul"][u ^ v][u ^ v] == f4["mul"][u][u] ^ f4["mul"][v][v] for u in range(4) for v in range(4)))


def check_sizes() -> None:
    section("How big can a finite field be")
    f4, z13 = ring_f4(), ring_zn(13)
    same("multiples of 1 in F_4", show_all(f4, multiples_of_one(f4)), "0 1")
    same("c1 * a + c0 * 1 in F_4 for c1 c0 = 00, 01, 10, 11",
         show_all(f4, [F4_MUL[c1][2] ^ F4_MUL[c0][1] for c1, c0 in [(0, 0), (0, 1), (1, 0), (1, 1)]]), "0 1 a b")
    same("multiples of 1 in Z_13", len(multiples_of_one(z13)), 13)
    same("sizes from 2 to 32 that are not prime powers", [n for n in range(2, 33) if not is_prime_power(n)],
         [6, 10, 12, 14, 15, 18, 20, 21, 22, 24, 26, 28, 30])
    same("sizes from 2 to 64: primes, higher prime powers, excluded",
         [sum(1 for n in range(2, 65) if is_prime(n)),
          sum(1 for n in range(2, 65) if is_prime_power(n) and not is_prime(n)),
          sum(1 for n in range(2, 65) if not is_prime_power(n))], [18, 9, 36])
    same("higher prime powers up to 64", [n for n in range(2, 65) if is_prime_power(n) and not is_prime(n)],
         [4, 8, 9, 16, 25, 27, 32, 49, 64])


def check_examples() -> None:
    section("Worked examples")
    z15 = ring_zn(15)
    u15 = units(z15)
    same("units of Z_15", u15, [1, 2, 4, 7, 8, 11, 13, 14])
    same("their inverses", [inverse_in(z15, u) for u in u15], [1, 8, 4, 13, 2, 11, 7, 14])
    same("zero divisors of Z_15", zero_divisors(z15), [3, 5, 6, 9, 10, 12])
    same("3 * 5, 6 * 5, 9 * 5, 10 * 3, 12 * 5 mod 15",
         [3 * 5 % 15, 6 * 5 % 15, 9 * 5 % 15, 10 * 3 % 15, 12 * 5 % 15], [0, 0, 0, 0, 0])
    same("6x = 10 in Z_256", solutions(6, 10, 256), [87, 215])
    same("6x = 7 in Z_256", solutions(6, 7, 256), [])
    same("inverse of 3 modulo 128", next(b for b in range(128) if 3 * b % 128 == 1), 43)
    same("5 * 43 mod 128", 5 * 43 % 128, 87)
    same("6 * 87 and 6 * 215 mod 256", [6 * 87 % 256, 6 * 215 % 256], [10, 10])
    m = ring_matrices()
    N, I = 0b1111, 0b1001
    same("N * N for N = 11|11", show(m, mat_mul(N, N)), "00|00")
    same("I + N", show(m, I ^ N), "01|10")
    same("(I + N)^2", show(m, mat_mul(I ^ N, I ^ N)), "10|01")
    f = 0b11010011
    same("a byte squared by wiring, bits", f"{bits(f)} -> {bits(spread(f))}", "11010011 -> 101000100000101")
    same("the same square by multiplying", bits(clmul(f, f)), "101000100000101")
    same("degree of the square", degree(spread(f)), 14)


def check_exercises() -> None:
    section("Exercises")
    z20 = ring_zn(20)
    same("units of Z_20", units(z20), [1, 3, 7, 9, 11, 13, 17, 19])
    same("zero divisors of Z_20", zero_divisors(z20), [2, 4, 5, 6, 8, 10, 12, 14, 15, 16, 18])
    same("12x = 6 in Z_18", solutions(12, 6, 18), [2, 5, 8, 11, 14, 17])
    same("12x = 4 in Z_18", solutions(12, 4, 18), [])
    same("gcd(12, 18)", math.gcd(12, 18), 6)
    w3 = ring_words(3)
    same("3-bit words: units, zero divisors, characteristic",
         f"{show_all(w3, units(w3))}; {len(zero_divisors(w3))}; {characteristic(w3)}", "111; 6; 2")
    same("C(7, i) for i = 0..7", pascal_row(7), [1, 7, 21, 35, 35, 21, 7, 1])
    same("C(6, i) for i = 0..6", pascal_row(6), [1, 6, 15, 20, 15, 6, 1])
    same("C(6, 2) mod 6, C(6, 3) mod 6", [15 % 6, 20 % 6], [3, 2])
    same("49, 50, 64, 100, 121 prime powers", " ".join(yes(is_prime_power(n)) for n in (49, 50, 64, 100, 121)),
         "yes no yes no yes")
    same("(3 + 2x)(3 + 2x) in Z_4[x], constant first", polymul([3, 2], [3, 2], 4), [1, 0, 0])
    same("F_4: solution of a x = 1 + a", F4_NAMES[next(x for x in range(4) if F4_MUL[2][x] == 1 ^ 2)], "a")


def main() -> None:
    rings = chapter_rings()
    check_opening()
    check_axioms(rings)
    check_gallery()
    check_consequences(rings)
    check_homomorphism()
    check_units(rings)
    check_solving()
    check_domains(rings)
    check_fields(rings)
    check_characteristic(rings)
    check_freshman()
    check_sizes()
    check_examples()
    check_exercises()
    print("\nAll chapter 32 numbers reproduced.")


if __name__ == "__main__":
    main()
