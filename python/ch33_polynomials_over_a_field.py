"""Reproduce every number in chapter 33 of Error Correction Coding, "Polynomials over a Field".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch33_polynomials_over_a_field.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees. All
the arithmetic is on whole numbers, so every comparison is exact.

The finite rings the chapter computes in, the integers modulo 4 and 8 and the
fields GF(2), GF(3) and GF(5), are first checked against the ring axioms by
enumerating every pair and every triple of their elements, and the three
fields against the field axioms. Polynomials over GF(2) are held as whole
numbers whose bit i is the coefficient of x^i, so that addition is XOR and
multiplication is the carry-less product. Polynomials over GF(3) and GF(5) are
held as lists of coefficients, constant term first.

Beside the worked numbers, each theorem the chapter proves is tested on every
small case: the degree of a product, that a polynomial is not its function,
the existence and uniqueness of quotient and remainder, that remainders
respect sums and products, that Euclid's algorithm returns the greatest common
divisor and that its cofactors satisfy Bezout's identity, the remainder and
factor theorems, the bound on the number of roots, the two root tests over
GF(2), the root test for irreducibility in degrees 2 and 3, the lemma that an
irreducible dividing a product divides a factor, and unique factorization. A
run that reaches the final line has reproduced every number in the chapter.

The module also holds the counts that ch33_charts.py draws, so the chart and
this check share one source.
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


# ---------------------------------------------------------------------------
# The finite rings the coefficients live in. Each is the integers modulo n.

def ring_axioms(n: int) -> bool:
    add = [[(a + b) % n for b in range(n)] for a in range(n)]
    mul = [[a * b % n for b in range(n)] for a in range(n)]
    for a in range(n):
        if add[a][0] != a or mul[a][1] != a:
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


def units(n: int) -> list[int]:
    return [a for a in range(1, n) if any(a * b % n == 1 for b in range(n))]


def zero_divisors(n: int) -> list[int]:
    return [a for a in range(1, n) if any(a * b % n == 0 for b in range(1, n))]


def inverse(a: int, p: int) -> int:
    for b in range(1, p):
        if a * b % p == 1:
            return b
    raise ValueError(f"{a} has no inverse modulo {p}")


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


def divmod2(f: int, g: int) -> tuple[int, int]:
    q = 0
    while f and deg2(f) >= deg2(g):
        s = deg2(f) - deg2(g)
        q |= 1 << s
        f ^= g << s
    return q, f


def division_steps2(f: int, g: int) -> tuple[list[int], list[str]]:
    """The shift of each step of long division and the remainder left after it."""
    shifts, rems = [], []
    while f and deg2(f) >= deg2(g):
        s = deg2(f) - deg2(g)
        f ^= g << s
        shifts.append(s)
        rems.append(bits(f))
    return shifts, rems


def mod2(f: int, g: int) -> int:
    return divmod2(f, g)[1]


def divides2(d: int, f: int) -> bool:
    return mod2(f, d) == 0


def bits(a: int) -> str:
    return format(a, "b")


def weight(a: int) -> int:
    return bin(a).count("1")


def ptxt2(a: int) -> str:
    if a == 0:
        return "0"
    terms = []
    for k in range(deg2(a), -1, -1):
        if a >> k & 1:
            terms.append("1" if k == 0 else "x" if k == 1 else f"x^{k}")
    return "+".join(terms)


def eval2(f: int, a: int) -> int:
    return weight(f) % 2 if a == 1 else f & 1


def xgcd2(a: int, b: int) -> list[list[int]]:
    """Rows [r, q, s, t] of the extended Euclidean table. Row i carries the
    quotient q_i of r_(i-1) by r_i, so r_(i+1) = r_(i-1) - q_i r_i, and the
    last row holds the remainder 0 and q = -1."""
    rows = [[a, -1, 1, 0], [b, -1, 0, 1]]
    while rows[-1][0] != 0:
        r0, _, s0, t0 = rows[-2]
        r1, _, s1, t1 = rows[-1]
        q, r = divmod2(r0, r1)
        rows[-1][1] = q
        rows.append([r, -1, s0 ^ clmul(q, s1), t0 ^ clmul(q, t1)])
    return rows


def gcd2(a: int, b: int) -> int:
    return xgcd2(a, b)[-2][0]


def factor_list2(f: int) -> list[int]:
    """The irreducible factors of f by trial division, smallest first, with repeats."""
    out = []
    d = 2
    while deg2(f) > 0:
        if deg2(d) * 2 > deg2(f):
            out.append(f)
            break
        q, r = divmod2(f, d)
        if r == 0:
            out.append(d)
            f = q
        else:
            d += 1
    return out


def factor_text2(f: int) -> str:
    fs = factor_list2(f)
    out = ""
    i = 0
    while i < len(fs):
        j = i
        while j < len(fs) and fs[j] == fs[i]:
            j += 1
        out += f"({ptxt2(fs[i])})" + (f"^{j - i}" if j - i > 1 else "")
        i = j
    return out


def irreducible_by_trial2(f: int) -> bool:
    if deg2(f) < 1:
        return False
    for d in range(2, 1 << (deg2(f) // 2 + 1)):
        if mod2(f, d) == 0:
            return False
    return True


def sieve2(top: int) -> list[int]:
    """Cross out every product of a known irreducible with a polynomial of degree
    at least 1. What survives up to degree top is irreducible."""
    limit = 1 << (top + 1)
    crossed = [False] * limit
    for a in range(2, limit):
        if crossed[a] or 2 * deg2(a) > top:
            continue
        for b in range(2, 1 << (top - deg2(a) + 1)):
            crossed[clmul(a, b)] = True
    return [a for a in range(2, limit) if not crossed[a]]


def irreducible_counts(top: int = 10) -> list[int]:
    irr = sieve2(top)
    return [sum(1 for a in irr if deg2(a) == m) for m in range(1, top + 1)]


# ---------------------------------------------------------------------------
# Polynomials over GF(p) or Z_n as coefficient lists, constant term first.

def trim(f: list[int]) -> list[int]:
    f = list(f)
    while f and f[-1] == 0:
        f.pop()
    return f


def pdeg(f: list[int]) -> int:
    return len(trim(f)) - 1


def padd(f: list[int], g: list[int], n: int) -> list[int]:
    m = max(len(f), len(g))
    return trim([((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % n for i in range(m)])


def pneg(f: list[int], n: int) -> list[int]:
    return trim([(-c) % n for c in f])


def conv(f: list[int], g: list[int]) -> list[int]:
    """The coefficients of the product before any reduction: c_k = sum of f_i g_(k-i)."""
    if not f or not g:
        return []
    c = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            c[i + j] += a * b
    return c


def pmul(f: list[int], g: list[int], n: int) -> list[int]:
    return trim([c % n for c in conv(f, g)])


def pdivmod(f: list[int], g: list[int], p: int) -> tuple[list[int], list[int]]:
    f, g = trim(f), trim(g)
    inv = inverse(g[-1], p)
    q = [0] * max(len(f) - len(g) + 1, 0)
    while f and len(f) >= len(g):
        s = len(f) - len(g)
        c = f[-1] * inv % p
        q[s] = c
        f = padd(f, pneg([0] * s + [c * x % p for x in g], p), p)
    return trim(q), f


def division_steps(f: list[int], g: list[int], p: int) -> tuple[list[str], list[str]]:
    f, g = trim(f), trim(g)
    inv = inverse(g[-1], p)
    terms, rems = [], []
    while f and len(f) >= len(g):
        s = len(f) - len(g)
        c = f[-1] * inv % p
        f = padd(f, pneg([0] * s + [c * x % p for x in g], p), p)
        terms.append(ptxt([0] * s + [c]))
        rems.append(ptxt(f))
    return terms, rems


def peval(f: list[int], a: int, n: int) -> int:
    v = 0
    for c in reversed(f):
        v = (v * a + c) % n
    return v


def roots(f: list[int], n: int) -> list[int]:
    return [a for a in range(n) if peval(f, a, n) == 0]


def ptxt(f: list[int]) -> str:
    f = trim(f)
    if not f:
        return "0"
    terms = []
    for k in range(len(f) - 1, -1, -1):
        c = f[k]
        if c == 0:
            continue
        if k == 0:
            terms.append(str(c))
        else:
            terms.append(("" if c == 1 else str(c)) + ("x" if k == 1 else f"x^{k}"))
    return "+".join(terms)


def pxgcd(a: list[int], b: list[int], p: int) -> list[list[list[int]]]:
    rows = [[trim(a), [], [1], []], [trim(b), [], [], [1]]]
    while rows[-1][0]:
        r0, _, s0, t0 = rows[-2]
        r1, _, s1, t1 = rows[-1]
        q, r = pdivmod(r0, r1, p)
        rows[-1][1] = q
        rows.append([r, [], padd(s0, pneg(pmul(q, s1, p), p), p), padd(t0, pneg(pmul(q, t1, p), p), p)])
    return rows


def scale(f: list[int], c: int, p: int) -> list[int]:
    return trim([c * x % p for x in f])


def from_index(i: int, p: int, length: int) -> list[int]:
    out = []
    for _ in range(length):
        out.append(i % p)
        i //= p
    return trim(out)


def monics(p: int, d: int) -> list[list[int]]:
    return [[i // p ** k % p for k in range(d)] + [1] for i in range(p ** d)]


# ---------------------------------------------------------------------------
# The checks, in the order the chapter states them.

def check_rings() -> None:
    section("The coefficient rings")
    for n, name in [(4, "Z_4"), (8, "Z_8"), (2, "GF(2)"), (3, "GF(3)"), (5, "GF(5)")]:
        require(f"ring axioms on all {n ** 3} triples of {name}", ring_axioms(n))
    for p in (2, 3, 5):
        same(f"nonzero elements of GF({p}) with an inverse", len(units(p)), p - 1)
    same("zero divisors of Z_4", joined(zero_divisors(4)), "2")
    same("zero divisors of Z_8", joined(zero_divisors(8)), "2 4 6")
    same("inverses of 1 2 3 4 in GF(5)", joined(inverse(a, 5) for a in range(1, 5)), "1 3 2 4")


def check_opening() -> None:
    section("The opening product")
    a, b = 0b1011, 0b0111
    same("1011 and 0111 as whole numbers", f"{a} {b}", "11 7")
    same("their integer product", a * b, 77)
    same("77 in binary", bits(a * b), "1001101")
    rows = [a << i for i in range(4) if b >> i & 1]
    same("the shifted copies of 1011", joined(bits(r) for r in rows), "1011 10110 101100")
    sums = [sum(r >> k & 1 for r in rows) for k in range(5, -1, -1)]
    same("column sums, from x^5 down to the constant", joined(sums), "1 1 2 2 2 1")
    carry, carries = 0, []
    for k in range(6):
        total = sum(r >> k & 1 for r in rows) + carry
        carry = total // 2
        carries.append(carry)
    same("carries out of the six columns, constant column first", joined(carries), "0 1 1 1 1 1")
    p = clmul(a, b)
    same("carry-less product of 1011 and 0111", bits(p), "110001")
    same("the carry-less product as a whole number", p, 49)
    same("(x^3+x+1)(x^2+x+1) over GF(2)", ptxt2(p), "x^5+x^4+1")
    same("column sums reduced modulo 2", joined(s % 2 for s in sums), "1 1 0 0 0 1")


def check_arithmetic() -> None:
    section("Addition and multiplication")
    same("(x^3+x+1) + (x^2+x+1) over GF(2)", ptxt2(0b1011 ^ 0b0111), "x^3+x^2")
    same("1011 XOR 0111", format(0b1011 ^ 0b0111, "04b"), "1100")
    f3, g3 = [1, 1, 2], [2, 1]
    same("(2x^2+x+1) + (x+2) over GF(3)", ptxt(padd(f3, g3, 3)), "2x^2+2x")
    same("convolution sums of 2x^2+x+1 and x+2, top first", joined(reversed(conv(f3, g3))), "2 5 3 2")
    same("(2x^2+x+1)(x+2) over GF(3)", ptxt(pmul(f3, g3, 3)), "2x^3+2x^2+2")
    f5, g5 = [2, 4, 3], [3, 2]
    same("(3x^2+4x+2) + (2x+3) over GF(5)", ptxt(padd(f5, g5, 5)), "3x^2+x")
    same("convolution sums of 3x^2+4x+2 and 2x+3, top first", joined(reversed(conv(f5, g5))), "6 17 16 6")
    same("(3x^2+4x+2)(2x+3) over GF(5)", ptxt(pmul(f5, g5, 5)), "x^3+2x^2+x+1")
    pairs, ok = 0, True
    for i in range(1, 27):
        for j in range(1, 27):
            f, g = from_index(i, 3, 3), from_index(j, 3, 3)
            pairs += 1
            ok = ok and pdeg(pmul(f, g, 3)) == pdeg(f) + pdeg(g)
    same("pairs of nonzero GF(3) polynomials of degree below 3", pairs, 676)
    require("deg(fg) = deg f + deg g on every one of those pairs", ok)
    same("(2x+1)^2 over Z_4", ptxt(pmul([1, 2], [1, 2], 4)), "1")
    sq = [ptxt([a0, a1]) for a1 in range(1, 4) for a0 in range(4) if pmul([a0, a1], [a0, a1], 4) == [1]]
    same("degree-1 polynomials over Z_4 whose square is 1", joined(sq), "2x+1 2x+3")
    same("(2x)(2x) over Z_4", ptxt(pmul([0, 2], [0, 2], 4)), "0")


def check_functions() -> None:
    section("A polynomial is not its function")
    f = 0b110
    same("x^2+x over GF(2) at 0 and at 1", f"{eval2(f, 0)} {eval2(f, 1)}", "0 0")
    vanish = [ptxt2(g) for g in range(16) if eval2(g, 0) == 0 and eval2(g, 1) == 0]
    same("polynomials of degree below 4 that vanish at 0 and 1", joined(vanish), "0 x^2+x x^3+x x^3+x^2")
    same("how many of the 16 that is", len(vanish), 4)
    same("functions from GF(2) to GF(2)", 2 ** 2, 4)
    require("each of those four is divisible by x^2+x", all(divides2(0b110, int_of(v)) for v in vanish))
    same("x^3-x over GF(3) at 0, 1 and 2", joined(peval([0, 2, 0, 1], a, 3) for a in range(3)), "0 0 0")
    same("x^5-x over GF(5) at 0 through 4", joined(peval([0, 4, 0, 0, 0, 1], a, 5) for a in range(5)), "0 0 0 0 0")


def int_of(text: str) -> int:
    """Read back a polynomial over GF(2) written by ptxt2."""
    if text == "0":
        return 0
    v = 0
    for term in text.split("+"):
        v |= 1 << (0 if term == "1" else 1 if term == "x" else int(term[2:]))
    return v


def check_division() -> None:
    section("The division algorithm")
    pairs, ok, unique = 0, True, True
    for f in range(64):
        for g in range(1, 16):
            q, r = divmod2(f, g)
            pairs += 1
            ok = ok and clmul(q, g) ^ r == f and deg2(r) < deg2(g)
            count = sum(1 for rr in range(1 << deg2(g)) if divides2(g, f ^ rr))
            unique = unique and count == 1
    same("pairs f of degree below 6, g nonzero of degree below 4", pairs, 960)
    require("f = qg + r with deg r < deg g on every pair", ok)
    require("exactly one remainder of degree below deg g on every pair", unique)
    require("2c = 1 has no integer solution, since 2c is even", all((2 * c) % 2 == 0 for c in range(-50, 51)))
    f, g = 0b10010011, 0b1011
    same("x^7+x^4+x+1 as a bit string", bits(f), "10010011")
    shifts, rems = division_steps2(f, g)
    same("shifts of 1011 at each step", joined(shifts), "4 2 0")
    same("what is left after each step", joined(rems), "100011 1111 100")
    q, r = divmod2(f, g)
    same("quotient of x^7+x^4+x+1 by x^3+x+1", ptxt2(q), "x^4+x^2+1")
    same("remainder", ptxt2(r), "x^2")
    same("(x^4+x^2+1)(x^3+x+1) over GF(2)", ptxt2(clmul(q, g)), "x^7+x^4+x^2+x+1")
    f5, g5 = [2, 4, 1, 3], [1, 3, 2]
    same("inverse of the leading coefficient 2 in GF(5)", inverse(2, 5), 3)
    terms, rems5 = division_steps(f5, g5, 5)
    same("quotient terms of 3x^3+x^2+4x+2 by 2x^2+3x+1 over GF(5)", joined(terms), "4x 2")
    same("what is left after each step", joined(rems5), "4x^2+2 4x")
    q5, r5 = pdivmod(f5, g5, 5)
    same("quotient", ptxt(q5), "4x+2")
    same("remainder", ptxt(r5), "4x")
    same("(4x+2)(2x^2+3x+1) + 4x over GF(5)", ptxt(padd(pmul(q5, g5, 5), r5, 5)), "3x^3+x^2+4x+2")


def check_remainders() -> None:
    section("Remainders respect sums and products")
    g = 0b1011
    a, b = 0b10000, 0b1100
    ra, rb = mod2(a, g), mod2(b, g)
    same("x^4 mod (x^3+x+1)", ptxt2(ra), "x^2+x")
    same("(x^3+x^2) mod (x^3+x+1)", ptxt2(rb), "x^2+x+1")
    same("(x^4)(x^3+x^2) = x^7+x^6, reduced", ptxt2(mod2(clmul(a, b), g)), "x^2")
    same("(x^2+x)(x^2+x+1), before reduction", ptxt2(clmul(ra, rb)), "x^4+x")
    same("(x^2+x)(x^2+x+1), reduced", ptxt2(mod2(clmul(ra, rb), g)), "x^2")
    same("x^4 + x^3 + x^2, reduced", ptxt2(mod2(a ^ b, g)), "1")
    same("(x^2+x) + (x^2+x+1)", ptxt2(ra ^ rb), "1")
    pairs, ok = 0, True
    for a in range(32):
        for b in range(32):
            pairs += 1
            ok = ok and mod2(a ^ b, g) == mod2(a, g) ^ mod2(b, g)
            ok = ok and mod2(clmul(a, b), g) == mod2(clmul(mod2(a, g), mod2(b, g)), g)
    same("pairs of polynomials of degree below 5", pairs, 1024)
    require("both rules hold modulo x^3+x+1 on every pair", ok)


def show_table2(name: str, rows: list[list[int]]) -> None:
    same(f"{name}: remainders r_i", joined(ptxt2(r[0]) for r in rows), stated_tables[name][0])
    same(f"{name}: quotients q_i", joined(ptxt2(r[1]) for r in rows if r[1] >= 0), stated_tables[name][1])
    same(f"{name}: cofactors s_i", joined(ptxt2(r[2]) for r in rows), stated_tables[name][2])
    same(f"{name}: cofactors t_i", joined(ptxt2(r[3]) for r in rows), stated_tables[name][3])


stated_tables = {
    "GF(2) table": ["x^5+x^4+1 x^4+x^2+1 x^3+x^2+x x^2+x+1 0",
                    "x+1 x+1 x",
                    "1 0 1 x+1 x^2+x+1",
                    "0 1 x+1 x^2 x^3+x+1"],
    "inverse of x^2+1": ["x^3+x+1 x^2+1 1 0", "x x^2+1", "1 0 1 x^2+1", "0 1 x x^3+x+1"],
    "inverse of x^2": ["x^4+x+1 x^2 x+1 1 0", "x^2 x+1 x+1", "1 0 1 x+1 x^2",
                       "0 1 x^2 x^3+x^2+1 x^4+x+1"],
    "exercise gcd": ["x^5+x+1 x^4+x^2+x+1 x^3+x^2+1 0", "x x+1", "1 0 1 x+1", "0 1 x x^2+x+1"],
}


def check_gcd() -> None:
    section("Greatest common divisors and Euclid's algorithm")
    a, b = 0b110001, 0b10101
    same("x^5+x^4+1, factored", factor_text2(a), "(x^2+x+1)(x^3+x+1)")
    same("x^4+x^2+1, factored", factor_text2(b), "(x^2+x+1)^2")
    rows = xgcd2(a, b)
    show_table2("GF(2) table", rows)
    g, s, t = rows[-2][0], rows[-2][2], rows[-2][3]
    same("gcd(x^5+x^4+1, x^4+x^2+1)", ptxt2(g), "x^2+x+1")
    same("(x+1)(x^5+x^4+1) + x^2(x^4+x^2+1)", ptxt2(clmul(s, a) ^ clmul(t, b)), "x^2+x+1")
    same("degrees of the remainders", joined(deg2(r[0]) for r in rows[:-1]), "5 4 3 2")
    divides = [[False] * 32] + [[divides2(d, f) for f in range(32)] for d in range(1, 32)]
    pairs, ok = 0, True
    for x in range(1, 32):
        for y in range(1, 32):
            pairs += 1
            common = [d for d in range(1, 32) if divides[d][x] and divides[d][y]]
            top = max(deg2(d) for d in common)
            best = [d for d in common if deg2(d) == top]
            rows = xgcd2(x, y)
            gg, ss, tt = rows[-2][0], rows[-2][2], rows[-2][3]
            ok = ok and best == [gg] and clmul(ss, x) ^ clmul(tt, y) == gg
            ok = ok and all(divides2(d, gg) for d in common)
    same("pairs of nonzero GF(2) polynomials of degree below 5", pairs, 961)
    require("Euclid's answer is the one common divisor of top degree", ok)
    require("and every common divisor divides it, with s a + t b equal to it", ok)
    a3, b3 = [1, 1, 1, 1], [1, 0, 2]
    same("x^3+x^2+x+1 over GF(3), factored by hand", ptxt(pmul([1, 1], [1, 0, 1], 3)), "x^3+x^2+x+1")
    same("(x+1)(2x+1) over GF(3)", ptxt(pmul([1, 1], [1, 2], 3)), "2x^2+1")
    rows3 = pxgcd(a3, b3, 3)
    same("GF(3) table: remainders r_i", joined(ptxt(r[0]) for r in rows3), "x^3+x^2+x+1 2x^2+1 2x+2 0")
    same("GF(3) table: quotients q_i", joined(ptxt(r[1]) for r in rows3[1:-1]), "2x+2 x+2")
    same("GF(3) table: cofactors s_i", joined(ptxt(r[2]) for r in rows3), "1 0 1 2x+1")
    same("GF(3) table: cofactors t_i", joined(ptxt(r[3]) for r in rows3), "0 1 x+1 2x^2+2")
    g3, s3, t3 = rows3[-2][0], rows3[-2][2], rows3[-2][3]
    same("1(x^3+x^2+x+1) + (x+1)(2x^2+1) over GF(3)", ptxt(padd(pmul(s3, a3, 3), pmul(t3, b3, 3), 3)), "2x+2")
    c = inverse(g3[-1], 3)
    same("inverse of the leading coefficient 2 in GF(3)", c, 2)
    same("the monic gcd", ptxt(scale(g3, c, 3)), "x+1")
    same("the monic cofactors s and t", f"{ptxt(scale(s3, c, 3))} {ptxt(scale(t3, c, 3))}", "2 2x+2")
    same("2(x^3+x^2+x+1) + (2x+2)(2x^2+1) over GF(3)",
         ptxt(padd(pmul(scale(s3, c, 3), a3, 3), pmul(scale(t3, c, 3), b3, 3), 3)), "x+1")
    rows = xgcd2(0b1011, 0b101)
    show_table2("inverse of x^2+1", rows)
    same("inverse of x^2+1 modulo x^3+x+1", ptxt2(rows[-2][3]), "x")
    same("x(x^2+1) = x^3+x, reduced modulo x^3+x+1", ptxt2(mod2(clmul(0b10, 0b101), 0b1011)), "1")
    m3, u3 = [1, 0, 1], [1, 1]
    rows = pxgcd(m3, u3, 3)
    same("x^2+1 and x+1 over GF(3): remainders r_i", joined(ptxt(r[0]) for r in rows), "x^2+1 x+1 2 0")
    same("x^2+1 and x+1 over GF(3): q_1, s_2 and t_2",
         f"{ptxt(rows[1][1])} {ptxt(rows[2][2])} {ptxt(rows[2][3])}", "x+2 1 2x+1")
    inv3 = scale(rows[2][3], inverse(rows[2][0][-1], 3), 3)
    same("inverse of x+1 modulo x^2+1 over GF(3), 2 times t_2", ptxt(inv3), "x+2")
    same("(x+1)(x+2) modulo x^2+1 over GF(3)", ptxt(pdivmod(pmul(u3, inv3, 3), m3, 3)[1]), "1")


def check_roots() -> None:
    section("Roots and factors")
    f = [1, 2, 0, 1]
    same("x^3+2x+1 at 2 over GF(5)", peval(f, 2, 5), 3)
    q, r = pdivmod(f, [3, 1], 5)
    same("x - 2 written over GF(5)", ptxt([3, 1]), "x+3")
    same("quotient of x^3+2x+1 by x+3 over GF(5)", ptxt(q), "x^2+2x+1")
    same("remainder", ptxt(r), "3")
    polys, ok, bound = 0, True, True
    for i in range(625):
        g = from_index(i, 5, 4)
        polys += 1
        for a in range(5):
            _, rr = pdivmod(g, [(-a) % 5, 1], 5)
            ok = ok and (rr[0] if rr else 0) == peval(g, a, 5)
        if g:
            bound = bound and len(roots(g, 5)) <= pdeg(g)
    same("polynomials over GF(5) of degree below 4", polys, 625)
    require("f mod (x - a) = f(a) for every one of them and every a", ok)
    require("no nonzero one has more roots than its degree", bound)
    same("x^4-1 over GF(5): its roots", joined(roots([4, 0, 0, 0, 1], 5)), "1 2 3 4")
    same("x^2-1 over Z_8: its roots", joined(roots([7, 0, 1], 8)), "1 3 5 7")
    same("(x-3)(x+3) over Z_8", ptxt(pmul([5, 1], [3, 1], 8)), "x^2+7")
    same("(x-1)(x+1) over Z_8", ptxt(pmul([7, 1], [1, 1], 8)), "x^2+7")
    root1, root0 = True, True
    for g in range(256):
        root1 = root1 and (eval2(g, 1) == 0) == (weight(g) % 2 == 0)
        root0 = root0 and (eval2(g, 0) == 0) == (g & 1 == 0)
        if g and eval2(g, 1) == 0:
            root1 = root1 and divides2(0b11, g)
    require("over GF(2): 1 is a root exactly when the weight is even", root1)
    require("over GF(2): 0 is a root exactly when the constant is 0", root0)
    same("x^3+x^2+x+1 over GF(2): weight and root 1", f"{weight(0b1111)} {yes(eval2(0b1111, 1) == 0)}", "4 yes")


def yes(flag: bool) -> str:
    return "yes" if flag else "no"


def check_irreducibles() -> None:
    section("Irreducible polynomials")
    for p in (2, 3, 5):
        ok = True
        for d in (2, 3):
            products = set()
            for g in monics(p, 1):
                for h in monics(p, d - 1):
                    products.add(tuple(pmul(g, h, p)))
            for f in monics(p, d):
                ok = ok and ((tuple(f) not in products) == (len(roots(f, p)) == 0))
        require(f"GF({p}), degrees 2 and 3: irreducible exactly when rootless", ok)
    same("x^4+x^2+1 over GF(2): roots", len([a for a in (0, 1) if eval2(0b10101, a) == 0]), 0)
    same("(x^2+x+1)^2 over GF(2)", ptxt2(clmul(0b111, 0b111)), "x^4+x^2+1")
    irr = sieve2(10)
    for d, stated in [(1, "10 11"), (2, "111"), (3, "1011 1101"), (4, "10011 11001 11111"),
                      (5, "100101 101001 101111 110111 111011 111101")]:
        same(f"irreducibles of degree {d} over GF(2)", joined(bits(a) for a in irr if deg2(a) == d), stated)
    counts = irreducible_counts(10)
    same("counts for degrees 1 to 5", joined(counts[:5]), "2 1 2 3 6")
    same("counts for degrees 6 to 10", joined(counts[5:]), "9 18 30 56 99")
    trial = [sum(1 for f in range(1 << m, 1 << (m + 1)) if irreducible_by_trial2(f)) for m in range(1, 9)]
    same("counts for degrees 1 to 8 by trial division", joined(trial), "2 1 2 3 6 9 18 30")
    by_trial = [f for f in range(2, 1 << 9) if irreducible_by_trial2(f)]
    require("trial division and the sieve list the same irreducibles through degree 8",
            by_trial == [a for a in irr if deg2(a) <= 8])
    same("polynomials of degree 10 over GF(2)", 2 ** 10, 1024)
    same("99 out of 1024, times 10, to three places", f"{10 * 99 / 1024:.3f}", "0.967")
    for d, stated in [(2, "2 1 0 1"), (3, "4 2 0 2"), (4, "8 4 1 3"), (5, "16 8 2 6")]:
        top = [f for f in range(1 << d, 1 << (d + 1))]
        by_x = sum(1 for f in top if f & 1 == 0)
        by_x1 = sum(1 for f in top if f & 1 and weight(f) % 2 == 0)
        by_q = sum(1 for f in top if f & 1 and weight(f) % 2 and f != 0b111 and divides2(0b111, f))
        kept = sum(1 for f in top if f in irr)
        same(f"degree {d}: crossed by x, by x+1, by x^2+x+1, kept", joined([by_x, by_x1, by_q, kept]), stated)
    same("x^2+x+1 times x^3+x^2+1", ptxt2(clmul(0b111, 0b1101)), "x^5+x+1")
    irr3 = [f for f in monics(3, 2) if not roots(f, 3)]
    same("monic irreducible quadratics over GF(3)", joined(ptxt(f) for f in irr3), "x^2+1 x^2+x+2 x^2+2x+2")
    same("monic quadratics over GF(3)", len(monics(3, 2)), 9)
    split3 = [(a, b) for a in range(3) for b in range(a, 3)]
    same("products (x-a)(x-b) with a <= b over GF(3)", len(split3), 6)
    same("x^2+x+2 over GF(3) at 0, 1 and 2", joined(peval([2, 1, 1], a, 3) for a in range(3)), "2 1 2")


def check_factorization() -> None:
    section("Unique factorization")
    irr = [a for a in sieve2(3)]
    checks, ok = 0, True
    for p in irr:
        for a in range(1, 16):
            for b in range(1, 16):
                checks += 1
                if divides2(p, clmul(a, b)):
                    ok = ok and (divides2(p, a) or divides2(p, b))
    same("irreducibles p of degree at most 3", len(irr), 5)
    same("triples p, a, b with a, b nonzero of degree below 4", checks, 1125)
    require("p | ab forces p | a or p | b on every triple", ok)
    same("x^2+1 = (x+1)^2 divides (x+1)(x+1) but not x+1", yes(divides2(0b101, clmul(0b11, 0b11)) and not divides2(0b101, 0b11)), "yes")
    irr6 = sieve2(6)
    hits = [0] * 128

    def build(start: int, product: int) -> None:
        hits[product] += 1
        for i in range(start, len(irr6)):
            if deg2(product) + deg2(irr6[i]) <= 6:
                build(i, clmul(product, irr6[i]))

    build(0, 1)
    same("products of irreducibles of total degree at most 6", sum(hits) - 1, 126)
    require("every polynomial of degree 1 to 6 arises exactly once", all(h == 1 for h in hits[2:]))
    same("x^4+1 over GF(2), factored", factor_text2(0b10001), "(x+1)^4")
    same("(x+1)^4 over GF(2)", ptxt2(clmul(clmul(0b11, 0b11), clmul(0b11, 0b11))), "x^4+1")
    same("x^5+1 over GF(2), factored", factor_text2(0b100001), "(x+1)(x^4+x^3+x^2+x+1)")
    require("x^4+x^3+x^2+x+1 is irreducible", irreducible_by_trial2(0b11111))


def check_hardware() -> None:
    section("Polynomial arithmetic in hardware")
    for m, stated in [(3, "3 9 4"), (8, "8 64 49")]:
        products = m * m
        xors = sum(len([(i, k - i) for i in range(m) if 0 <= k - i < m]) - 1 for k in range(2 * m - 1))
        same(f"{m}-bit operands: adder XORs, multiplier ANDs and XORs", f"{m} {products} {xors}", stated)
    same("x^3 mod (x^3+x+1)", ptxt2(mod2(0b1000, 0b1011)), "x+1")
    same("x^4 mod (x^3+x+1)", ptxt2(mod2(0b10000, 0b1011)), "x^2+x")
    ok = True
    for c in range(32):
        bit = [c >> k & 1 for k in range(5)]
        wired = (bit[0] ^ bit[3]) | (bit[1] ^ bit[3] ^ bit[4]) << 1 | (bit[2] ^ bit[4]) << 2
        ok = ok and wired == mod2(c, 0b1011)
    require("r0 = c0+c3, r1 = c1+c3+c4, r2 = c2+c4 on all 32 inputs", ok)
    same("XOR gates in that reduction, one per tap bit",
         weight(mod2(0b1000, 0b1011)) + weight(mod2(0b10000, 0b1011)), 4)
    # A shift-register divider for 1011: s2' = s1, s1' = s0 + s2, s0' = in + s2.
    state, quotient = 0, ""
    for k in range(7, -1, -1):
        fb = state >> 2 & 1
        state = (state << 1 & 0b111) | (0b10010011 >> k & 1)
        if fb:
            state ^= 0b011
        quotient += str(fb)
    same("divider for 1011 on 10010011: quotient and remainder", f"{quotient[3:]} {bits(state)}", "10101 100")
    same("XOR gates in that divider, ones of 1011 below the top", weight(0b1011 ^ 0b1000), 2)
    same("x^3 mod (x^3+x^2+1) and x^4 mod (x^3+x^2+1)",
         f"{ptxt2(mod2(0b1000, 0b1101))} {ptxt2(mod2(0b10000, 0b1101))}", "x^2+1 x^2+x+1")


def check_examples() -> None:
    section("Worked examples")
    f, g = [2, 1, 0, 2, 1], [1, 1, 2]
    terms, rems = division_steps(f, g, 3)
    same("x^4+2x^3+x+2 by 2x^2+x+1 over GF(3): quotient terms", joined(terms), "2x^2 2")
    same("what is left after each step", joined(rems), "x^2+x+2 2x")
    q, r = pdivmod(f, g, 3)
    same("quotient and remainder", f"{ptxt(q)} {ptxt(r)}", "2x^2+2 2x")
    rows = xgcd2(0b10011, 0b100)
    show_table2("inverse of x^2", rows)
    inv = rows[-2][3]
    same("inverse of x^2 modulo x^4+x+1", ptxt2(inv), "x^3+x^2+1")
    same("x^2(x^3+x^2+1) before reduction", ptxt2(clmul(0b100, inv)), "x^5+x^4+x^2")
    same("reduced modulo x^4+x+1", ptxt2(mod2(clmul(0b100, inv), 0b10011)), "1")
    f = 0b1010011
    same("x^6+x^4+x+1: weight and constant term", f"{weight(f)} {f & 1}", "4 1")
    q1, r1 = divmod2(f, 0b11)
    same("divided by x+1", f"{ptxt2(q1)} {ptxt2(r1)}", "x^5+x^4+1 0")
    same("x^5+x^4+1 at 1", eval2(q1, 1), 1)
    q2, r2 = divmod2(q1, 0b111)
    same("x^5+x^4+1 divided by x^2+x+1", f"{ptxt2(q2)} {ptxt2(r2)}", "x^3+x+1 0")
    same("x^6+x^4+x+1, factored", factor_text2(f), "(x+1)(x^2+x+1)(x^3+x+1)")
    f = [1, 1, 0, 1]
    same("x^3+x+1 over GF(5) at 0 to 4", joined(peval(f, a, 5) for a in range(5)), "1 3 1 1 4")
    same("x^3+x+1 over GF(3) at 0 to 2", joined(peval(f, a, 3) for a in range(3)), "1 0 2")
    q, r = pdivmod(f, [2, 1], 3)
    same("x^3+x+1 divided by x-1 = x+2 over GF(3)", f"{ptxt(q)} {ptxt(r)}", "x^2+x+2 0")
    same("x^3+x+1 over GF(3), the quadratic factor at 0 to 2", joined(peval(q, a, 3) for a in range(3)), "2 1 2")


def check_exercises() -> None:
    section("Exercises")
    same("101 times 111, carry-less", bits(clmul(0b101, 0b111)), "11011")
    same("as a polynomial", ptxt2(clmul(0b101, 0b111)), "x^4+x^3+x+1")
    same("5 times 7 as integers, in binary", bits(5 * 7), "100011")
    q, r = divmod2(0b1000011, 0b111)
    same("x^6+x+1 by x^2+x+1: quotient and remainder", f"{ptxt2(q)} {ptxt2(r)}", "x^4+x^3+x+1 x")
    show_table2("exercise gcd", xgcd2(0b100011, 0b10111))
    same("x^5+x+1, factored", factor_text2(0b100011), "(x^2+x+1)(x^3+x^2+1)")
    same("x^4+x^2+x+1, factored", factor_text2(0b10111), "(x+1)(x^3+x^2+1)")
    same("x^3+2x+1 over GF(3) at 0, 1 and 2", joined(peval([1, 2, 0, 1], a, 3) for a in range(3)), "1 1 1")
    same("roots of x^2-1 in Z_15", joined(roots([14, 0, 1], 15)), "1 4 11 14")
    rows = xgcd2(0b1011, 0b10)
    same("inverse of x modulo x^3+x+1", ptxt2(rows[-2][3]), "x^2+1")
    same("x^4+x^3+x^2+1 over GF(2), factored", factor_text2(0b11101), "(x+1)(x^3+x+1)")
    same("(2x+1)(3x+1) over Z_6", ptxt(pmul([1, 2], [1, 3], 6)), "5x+1")
    irr5 = [f for f in monics(5, 2) if not roots(f, 5)]
    same("monic irreducible quadratics over GF(5)", len(irr5), 10)
    split5 = sum(1 for f in monics(5, 2) if roots(f, 5))
    same("monic quadratics over GF(5), and the split ones", f"{len(monics(5, 2))} {split5}", "25 15")


def main() -> None:
    check_rings()
    check_opening()
    check_arithmetic()
    check_functions()
    check_division()
    check_remainders()
    check_gcd()
    check_roots()
    check_irreducibles()
    check_factorization()
    check_hardware()
    check_examples()
    check_exercises()
    print("\nAll chapter 33 numbers reproduced.")


if __name__ == "__main__":
    main()
