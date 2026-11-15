"""Reproduce every number in chapter 13 of Error Correction Coding, "Integers and Modular Arithmetic".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch13_integers_and_modular_arithmetic.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees. All
the arithmetic is on whole numbers, so every comparison is exact. Beside the
worked numbers, each theorem the chapter proves is also tested by brute force
on small cases: the division theorem, the step bound and Lame's theorem for
Euclid's algorithm, Bezout's identity, the rules of congruence arithmetic, the
existence and uniqueness of inverses, the theorems of Fermat and Euler, the
Chinese remainder theorem and the multiplicativity of the totient, and the fact
that the integers modulo n form a field exactly when n is prime. A run that
reaches the final line has reproduced every number in the chapter.

The module also holds the arithmetic that ch13_charts.py draws, so the charts
and this check share one source.
"""

import math
import sys


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


# The arithmetic ---------------------------------------------------------------

def floor_divmod(a: int, b: int) -> list:
    """Quotient rounded down, as Python's // and % give it."""
    q = a // b
    return [q, a - q * b]


def trunc_divmod(a: int, b: int) -> list:
    """Quotient rounded toward zero, as C99 and Octave's fix and rem give it."""
    q = abs(a) // abs(b)
    if (a < 0) != (b < 0):
        q = -q
    return [q, a - q * b]


def euclid_divmod(a: int, b: int) -> list:
    """The remainder that is never negative, whatever the signs."""
    r = a % abs(b)
    return [(a - r) // b, r]


def divisors(n: int) -> list:
    return [d for d in range(1, n + 1) if n % d == 0]


def euclid_steps(a: int, b: int) -> list:
    """Each division of Euclid's algorithm as [dividend, quotient, divisor, remainder]."""
    steps = []
    while b != 0:
        q, r = a // b, a % b
        steps.append([a, q, b, r])
        a, b = b, r
    return steps


def step_bound(b: int) -> int:
    """2(floor(log2 b) + 1), the most divisions the halving lemma allows."""
    k = 0
    while 2 ** (k + 1) <= b:
        k += 1
    return 2 * (k + 1)


def extended_rows(a: int, b: int) -> list:
    """Rows [r, q, s, t] of the extended algorithm, with r = s a + t b on every row.

    The q on a row is the quotient that row's r divides into the r above it, so
    the first and last rows carry none, which is written as 0 here.
    """
    rows = [[a, 0, 1, 0], [b, 0, 0, 1]]
    while rows[-1][0] != 0:
        r0, _, s0, t0 = rows[-2]
        r1, _, s1, t1 = rows[-1]
        q = r0 // r1
        rows[-1][1] = q
        rows.append([r0 - q * r1, 0, s0 - q * s1, t0 - q * t1])
    return rows


def bezout(a: int, b: int) -> list:
    """[gcd, s, t] with gcd = s a + t b, read from the last nonzero row."""
    rows = extended_rows(a, b)
    r, _, s, t = rows[-2]
    return [r, s, t]


def inverse(a: int, n: int):
    """The inverse of a modulo n in 0..n-1, or the word none when gcd(a, n) > 1."""
    g, s, _ = bezout(a % n, n) if a % n else [n, 0, 0]
    if g != 1:
        return "none"
    return s % n


def power_mod(a: int, e: int, n: int) -> int:
    """a^e mod n by repeated squaring, one squaring per bit of e."""
    result, square = 1 % n, a % n
    while e:
        if e & 1:
            result = result * square % n
        square = square * square % n
        e >>= 1
    return result


def squares(a: int, count: int, n: int) -> list:
    """a, a^2, a^4, ..., a^(2^(count-1)) modulo n."""
    out, x = [], a % n
    for _ in range(count):
        out.append(x)
        x = x * x % n
    return out


def binary(e: int) -> str:
    return format(e, "b")


def units(n: int) -> list:
    return [x for x in range(1, n + 1) if math.gcd(x, n) == 1] if n > 1 else [1]


def phi(n: int) -> int:
    return len(units(n))


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))


def crt(a: int, m: int, b: int, n: int) -> int:
    """The x in 0..mn-1 with x = a mod m and x = b mod n, for coprime m and n."""
    _, s, t = bezout(m, n)
    return (a * t * n + b * s * m) % (m * n)


def solutions(c: int, d: int, n: int) -> list:
    """Every x in 0..n-1 with c x = d mod n."""
    return [x for x in range(n) if (c * x - d) % n == 0]


def order(a: int, n: int) -> int:
    k, x = 1, a % n
    while x != 1:
        x = x * a % n
        k += 1
    return k


def is_field(n: int) -> bool:
    """Every nonzero element of Z_n has an inverse, and 1 differs from 0."""
    return n >= 2 and all(any(a * b % n == 1 for b in range(1, n)) for a in range(1, n))


def mult_table(n: int) -> list:
    return [[a * b % n for b in range(n)] for a in range(n)]


def isbn10_terms(digits: str) -> list:
    return [(10 - i) * int(x) for i, x in enumerate(digits)]


def luhn_double(x: int) -> int:
    return 2 * x if x <= 4 else 2 * x - 9


def luhn_terms(number: str) -> list:
    """The Luhn terms from the right: digit 1 as it is, digit 2 doubled, and so on."""
    digits = [int(c) for c in reversed(number)]
    return [luhn_double(y) if j % 2 == 1 else y for j, y in enumerate(digits)]


def fibonacci(k: int) -> int:
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


# The checks -------------------------------------------------------------------

def check_division() -> None:
    section("Division with remainder")
    same("17 divided by 5, quotient and remainder", floor_divmod(17, 5), [3, 2])
    table = {(17, 5): ([3, 2], [3, 2], [3, 2]),
             (-17, 5): ([-4, 3], [-3, -2], [-4, 3]),
             (17, -5): ([-4, -3], [-3, 2], [-3, 2]),
             (-17, -5): ([3, -2], [3, -2], [4, 3])}
    for (a, b), (fl, tr, eu) in table.items():
        same(f"{a} by {b}, floored q and r", floor_divmod(a, b), fl)
        same(f"{a} by {b}, truncated q and r", trunc_divmod(a, b), tr)
        same(f"{a} by {b}, Euclidean q and r", euclid_divmod(a, b), eu)
    same("-17 by 5 with the built-in floored operators", [-17 // 5, -17 % 5], [-4, 3])
    same("C-style remainder of -17 by 5 (built-in)", int(math.fmod(-17, 5)), -2)
    unique = all(sum(1 for q in range(-60, 61) for r in range(b) if a == q * b + r) == 1
                 for a in range(-50, 51) for b in range(1, 13))
    require("one q and r with 0 <= r < b, a in -50..50, b in 1..12", unique)
    lemma = True
    for a in range(-30, 31):
        for b in [v for v in range(-30, 31) if v != 0]:
            rf, rt = floor_divmod(a, b)[1], trunc_divmod(a, b)[1]
            lemma &= rf == (rt if rt == 0 or (rt > 0) == (b > 0) else rt + b)
    require("floored r = truncated r, plus b when the signs differ", lemma)


def check_gcd() -> None:
    section("The greatest common divisor")
    same("divisors of 1147", divisors(1147), [1, 31, 37, 1147])
    same("divisors of 899", divisors(899), [1, 29, 31, 899])
    same("common divisors of 1147 and 899",
         [d for d in divisors(1147) if 899 % d == 0], [1, 31])
    steps = euclid_steps(1147, 899)
    same("Euclid on 1147 and 899, quotients", [s[1] for s in steps], [1, 3, 1, 1, 1, 2])
    same("Euclid on 1147 and 899, remainders", [s[3] for s in steps], [248, 155, 93, 62, 31, 0])
    same("gcd(1147, 899)", math.gcd(1147, 899), 31)
    same("divisions taken", len(steps), 6)
    same("1147 / 31 and 899 / 31", [1147 // 31, 899 // 31], [37, 29])
    same("step bound 2(floor(log2 899) + 1)", step_bound(899), 20)
    halving = True
    for a, b in [(1147, 899), (4181, 2584)]:
        rs = [a, b] + [s[3] for s in euclid_steps(a, b)]
        halving &= all(2 * rs[i + 2] < rs[i] for i in range(len(rs) - 2))
    require("each remainder under half the one two places before", halving)
    bound_ok = all(len(euclid_steps(a, b)) <= step_bound(b)
                   for a in range(1, 300) for b in range(1, 300))
    require("divisions <= 2(floor(log2 b) + 1) for a, b below 300", bound_ok)
    smallest = {}
    for a in range(2, 234):
        for b in range(1, a):
            n = len(euclid_steps(a, b))
            smallest.setdefault(n, []).append((a, b))
    lame = all(min(p[0] for p in pairs) == fibonacci(n + 2)
               and [p for p in pairs if p[0] == fibonacci(n + 2)] == [(fibonacci(n + 2), fibonacci(n + 1))]
               for n, pairs in smallest.items())
    require("smallest pair needing n divisions is F(n+2), F(n+1)", lame)
    same("Fibonacci numbers F(1) to F(12)", [fibonacci(k) for k in range(1, 13)],
         [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144])


def check_bezout() -> None:
    section("Bezout's identity and the extended algorithm")
    same("11 * 1147 and 14 * 899", [11 * 1147, 14 * 899], [12617, 12586])
    same("11 * 1147 - 14 * 899", 11 * 1147 - 14 * 899, 31)
    rows = extended_rows(1147, 899)
    same("extended table, r", [r[0] for r in rows], [1147, 899, 248, 155, 93, 62, 31, 0])
    same("extended table, q (0 where none)", [r[1] for r in rows], [0, 1, 3, 1, 1, 1, 2, 0])
    same("extended table, s", [r[2] for r in rows], [1, 0, 1, -3, 4, -7, 11, -29])
    same("extended table, t", [r[3] for r in rows], [0, 1, -1, 4, -5, 9, -14, 37])
    require("r = 1147 s + 899 t on every row", all(r == 1147 * s + 899 * t for r, _, s, t in rows))
    same("-29 * 1147 + 37 * 899", -29 * 1147 + 37 * 899, 0)
    smallest = min(v for s in range(-60, 61) for t in range(-60, 61) if (v := 1147 * s + 899 * t) > 0)
    same("smallest positive 1147 s + 899 t, |s|, |t| <= 60", smallest, 31)
    bez = all(bezout(a, b)[0] == math.gcd(a, b) == bezout(a, b)[1] * a + bezout(a, b)[2] * b
              for a in range(1, 80) for b in range(1, 80))
    require("bezout gives gcd = s a + t b for a, b below 80", bez)
    same("36 mod 6, 4 mod 6, 9 mod 6", [36 % 6, 4 % 6, 9 % 6], [0, 4, 3])
    euclid_lemma = all((a * b) % p != 0 or a % p == 0 or b % p == 0
                       for p in range(2, 40) if is_prime(p) for a in range(60) for b in range(60))
    require("p | ab forces p | a or p | b, primes below 40", euclid_lemma)


def check_congruences() -> None:
    section("Congruences")
    same("9 + 5 and its remainder modulo 12", [9 + 5, (9 + 5) % 12], [14, 2])
    same("-17 - 3", -17 - 3, -20)
    same("remainders of -17 and 3 modulo 5", [-17 % 5, 3 % 5], [3, 3])
    same("37 mod 7 and 52 mod 7", [37 % 7, 52 % 7], [2, 3])
    same("37 * 52 and its remainder modulo 7", [37 * 52, 37 * 52 % 7], [1924, 6])
    same("1924 - 6", 1924 - 6, 274 * 7)
    respects = all(((a + b) - (c + d)) % n == 0 and (a * b - c * d) % n == 0
                   for n in range(1, 13) for a in range(-12, 13) for c in range(-12, 13)
                   if (a - c) % n == 0 for b in range(-6, 7) for d in range(-6, 7) if (b - d) % n == 0)
    require("congruence respects + and x for n up to 12", respects)
    same("13 in binary", binary(13), "1101")
    same("3, 3^2, 3^4, 3^8 modulo 11", squares(3, 4, 11), [3, 9, 4, 5])
    same("3^8 * 3^4 * 3 modulo 11", 5 * 4 * 3 % 11, 5)
    same("5 * 4 * 3", 5 * 4 * 3, 60)
    same("3^13", 3 ** 13, 1594323)
    same("3^13 mod 11", power_mod(3, 13, 11), 5)
    terms = isbn10_terms("0130426725")
    same("ISBN-10 terms of 0-13-042672-5", terms, [0, 9, 24, 0, 24, 10, 24, 21, 4, 5])
    same("ISBN-10 weighted sum and its remainder modulo 11", [sum(terms), sum(terms) % 11], [121, 0])
    same("Luhn sum of 3141592653 and its remainder modulo 10",
         [sum(luhn_terms("3141592653")), sum(luhn_terms("3141592653")) % 10], [40, 0])
    same("2 * 3 mod 6 and 2 * 0 mod 6", [2 * 3 % 6, 2 * 0 % 6], [0, 0])
    cancel = all((c * x - c * y) % n != 0 or (x - y) % n == 0
                 for n in range(2, 31) for c in range(1, n) if math.gcd(c, n) == 1
                 for x in range(n) for y in range(n))
    require("c x = c y forces x = y when gcd(c, n) = 1, n up to 30", cancel)


def check_inverses() -> None:
    section("Inverses modulo n")
    same("3 * 5 mod 7", 3 * 5 % 7, 1)
    same("inverses of 1 to 6 modulo 7", [inverse(a, 7) for a in range(1, 7)], [1, 4, 5, 2, 3, 6])
    steps = euclid_steps(37, 29)
    same("Euclid on 37 and 29, quotients", [s[1] for s in steps], [1, 3, 1, 1, 1, 2])
    same("Euclid on 37 and 29, remainders", [s[3] for s in steps], [8, 5, 3, 2, 1, 0])
    same("gcd, s, t for 37 and 29", bezout(37, 29), [1, 11, -14])
    same("11 * 37 - 14 * 29", 11 * 37 - 14 * 29, 1)
    same("-14 mod 37", -14 % 37, 23)
    same("29 * 23 and 18 * 37", [29 * 23, 18 * 37], [667, 666])
    same("inverse of 29 modulo 37", inverse(29, 37), 23)
    same("inverse of 29 modulo 37 (built-in)", pow(29, -1, 37), 23)
    terms = isbn10_terms("0130426725")
    same("ISBN-10 sum without position 5", sum(terms) - terms[4], 97)
    same("-97 mod 11", -97 % 11, 2)
    same("inverse of 6 modulo 11", inverse(6, 11), 2)
    same("recovered digit 2 * 2 mod 11", 2 * 2 % 11, 4)
    same("Z_6 elements with inverses", [a for a in range(1, 6) if inverse(a, 6) != "none"], [1, 5])
    same("5 * 5 mod 6", 5 * 5 % 6, 1)
    same("2 * 3 mod 6 and 3 * 4 mod 6", [2 * 3 % 6, 3 * 4 % 6], [0, 0])
    same("gcd(899, 1147)", math.gcd(899, 1147), 31)
    same("inverse of 899 modulo 1147", inverse(899, 1147), "none")
    exists = all((sum(1 for b in range(n) if a * b % n == 1) == 1) == (math.gcd(a, n) == 1)
                 and sum(1 for b in range(n) if a * b % n == 1) <= 1
                 for n in range(2, 41) for a in range(n))
    require("one inverse when gcd(a, n) = 1, none otherwise, n up to 40", exists)


def check_fermat_euler() -> None:
    section("Fermat's and Euler's theorems")
    rows = [[power_mod(a, k, 7) for k in range(1, 7)] for a in range(1, 7)]
    same("powers of 1 modulo 7", rows[0], [1, 1, 1, 1, 1, 1])
    same("powers of 2 modulo 7", rows[1], [2, 4, 1, 2, 4, 1])
    same("powers of 3 modulo 7", rows[2], [3, 2, 6, 4, 5, 1])
    same("powers of 4 modulo 7", rows[3], [4, 2, 1, 4, 2, 1])
    same("powers of 5 modulo 7", rows[4], [5, 4, 6, 2, 3, 1])
    same("powers of 6 modulo 7", rows[5], [6, 1, 6, 1, 6, 1])
    same("3k mod 7 for k = 1 to 6", [3 * k % 7 for k in range(1, 7)], [3, 6, 2, 5, 1, 4])
    same("6! and its remainder modulo 7", [math.factorial(6), math.factorial(6) % 7], [720, 6])
    fermat = all(power_mod(a, p - 1, p) == 1 for p in range(2, 60) if is_prime(p) for a in range(1, p))
    require("a^(p-1) = 1 mod p for primes below 60", fermat)
    little = all(power_mod(a, p, p) == a % p for p in range(2, 60) if is_prime(p) for a in range(-20, 40))
    require("a^p = a mod p for every a, primes below 60", little)
    same("units modulo 12", units(12), [1, 5, 7, 11])
    same("phi(12), phi(7), phi(9), phi(8)", [phi(12), phi(7), phi(9), phi(8)], [4, 6, 6, 4])
    same("units modulo 9", units(9), [1, 2, 4, 5, 7, 8])
    prime_power = all(phi(p ** k) == p ** k - p ** (k - 1)
                      for p in range(2, 1001) if is_prime(p) for k in range(1, 11) if p ** k <= 1000)
    require("phi(p^k) = p^k - p^(k-1) for p^k up to 1000", prime_power)
    same("2^6 and its remainder modulo 9", [2 ** 6, 2 ** 6 % 9], [64, 1])
    same("squares of 1, 5, 7, 11 modulo 12", [a * a % 12 for a in units(12)], [1, 1, 1, 1])
    euler = all(power_mod(a, phi(n), n) == 1 for n in range(2, 61) for a in units(n))
    require("a^phi(n) = 1 mod n when gcd(a, n) = 1, n up to 60", euler)
    same("100 = 6q + r", list(divmod(100, 6)), [16, 4])
    same("3^4 and its remainder modulo 7", [3 ** 4, 3 ** 4 % 7], [81, 4])
    same("3^100 mod 7", power_mod(3, 100, 7), 4)
    same("6, 6^2, 6^4, 6^8 modulo 11", squares(6, 4, 11), [6, 3, 9, 4])
    same("6^9 and its remainder modulo 11", [6 ** 9, 6 ** 9 % 11], [10077696, 2])
    same("11 * 31", 11 * 31, 341)
    same("2^10 and its remainder modulo 341", [2 ** 10, 2 ** 10 % 341], [1024, 1])
    same("3 * 341", 3 * 341, 1023)
    same("2^340 mod 341", power_mod(2, 340, 341), 1)
    same("3^340 mod 341", power_mod(3, 340, 341), 56)


def check_crt() -> None:
    section("The Chinese remainder theorem")
    same("x mod 3 for x = 3, 8, 13", [x % 3 for x in (3, 8, 13)], [0, 2, 1])
    same("x = 2 mod 3 and x = 3 mod 5", crt(2, 3, 3, 5), 8)
    grid = sorted((x % 3, x % 5) for x in range(15))
    require("each pair (x mod 3, x mod 5) once for x = 0 to 14",
            grid == [(i, j) for i in range(3) for j in range(5)])
    same("gcd, s, t for 3 and 5", bezout(3, 5), [1, 2, -1])
    same("2 * (-1) * 5 + 3 * 2 * 3", 2 * (-1) * 5 + 3 * 2 * 3, 8)
    same("gcd, s, t for 4 and 9", bezout(4, 9), [1, -2, 1])
    same("3 * 1 * 9 + 5 * (-2) * 4", 3 * 1 * 9 + 5 * (-2) * 4, -13)
    same("-13 mod 36", -13 % 36, 23)
    same("23 mod 4 and 23 mod 9", [23 % 4, 23 % 9], [3, 5])
    same("x = 3 mod 4 and x = 5 mod 9", crt(3, 4, 5, 9), 23)
    bijective = all(sorted((x % m, x % n) for x in range(m * n)) == [(i, j) for i in range(m) for j in range(n)]
                    and all(crt(a, m, b, n) % m == a and crt(a, m, b, n) % n == b
                            for a in range(m) for b in range(n))
                    for m in range(1, 16) for n in range(1, 16) if math.gcd(m, n) == 1)
    require("x -> (x mod m, x mod n) one to one when gcd(m, n) = 1", bijective)
    lemma = all(x % (m * n) == 0 for m in range(1, 13) for n in range(1, 13) if math.gcd(m, n) == 1
                for x in range(0, 400) if x % m == 0 and x % n == 0)
    require("m | x and n | x give mn | x when gcd(m, n) = 1", lemma)
    same("phi(4), phi(3), phi(12)", [phi(4), phi(3), phi(12)], [2, 2, 4])
    same("phi(4), phi(9), phi(36)", [phi(4), phi(9), phi(36)], [2, 6, 12])
    mult = all(phi(m * n) == phi(m) * phi(n) for m in range(1, 31) for n in range(1, 31) if math.gcd(m, n) == 1)
    require("phi(mn) = phi(m) phi(n) when gcd(m, n) = 1, up to 30", mult)
    same("x = 1 mod 4 and x = 2 mod 6, x in 0..11", [x for x in range(12) if x % 4 == 1 and x % 6 == 2], [])
    same("x = 1 mod 4 and x = 3 mod 6, x in 0..11", [x for x in range(12) if x % 4 == 1 and x % 6 == 3], [9])


def check_fields() -> None:
    section("The integers modulo a prime")
    same("elements of Z_7 and Z_6 with inverses", [phi(7), phi(6)], [6, 2])
    t7 = mult_table(7)
    require("each nonzero row of Z_7 holds 1 exactly once", all(row.count(1) == 1 for row in t7[1:]))
    t6 = mult_table(6)
    zero_pairs = [[a, b] for a in range(1, 6) for b in range(1, 6) if t6[a][b] == 0]
    same("pairs of nonzero elements of Z_6 with product 0", zero_pairs, [[2, 3], [3, 2], [3, 4], [4, 3]])
    same("rows of Z_6 holding no 1", [a for a in range(1, 6) if 1 not in t6[a]], [2, 3, 4])
    require("Z_n is a field exactly when n is prime, n = 2 to 60",
            all(is_field(n) == is_prime(n) for n in range(2, 61)))
    require("a + b mod 2 is a XOR b for all four pairs",
            all((a + b) % 2 == a ^ b for a in (0, 1) for b in (0, 1)))
    require("a b mod 2 is a AND b for all four pairs",
            all((a * b) % 2 == a & b for a in (0, 1) for b in (0, 1)))
    same("1 + 0 + 1 + 1 and its remainder modulo 2", [1 + 0 + 1 + 1, (1 + 0 + 1 + 1) % 2], [3, 1])
    same("2 * 128 mod 256", 2 * 128 % 256, 0)
    same("phi(256)", phi(256), 128)
    require("257 is prime", is_prime(257))


def check_examples() -> None:
    section("Worked examples")
    steps = euclid_steps(4181, 2584)
    same("F(19) and F(18)", [fibonacci(19), fibonacci(18)], [4181, 2584])
    same("Euclid on 4181 and 2584, divisions", len(steps), 17)
    same("Euclid on 4181 and 2584, quotients", [s[1] for s in steps], [1] * 16 + [2])
    same("Euclid on 4181 and 2584, remainders",
         [s[3] for s in steps], [fibonacci(k) for k in range(17, 1, -1)] + [0])
    same("gcd(4181, 2584)", math.gcd(4181, 2584), 1)
    same("step bound 2(floor(log2 2584) + 1)", step_bound(2584), 24)
    same("Luhn doubling of 0 to 9", [luhn_double(x) for x in range(10)], [0, 2, 4, 6, 8, 1, 3, 5, 7, 9])
    terms = luhn_terms("3141592653")
    same("Luhn terms of 3141592653 from the right", terms, [3, 1, 6, 4, 9, 1, 1, 8, 1, 6])
    same("Luhn sum without digit 2 from the right", sum(terms) - terms[1], 39)
    same("-39 mod 10", -39 % 10, 1)
    same("digits x with D(x) = 1", [x for x in range(10) if luhn_double(x) == 1], [5])
    same("solutions of 2x = 1 mod 10", solutions(2, 1, 10), [])
    require("D(x) = 2x mod 9 for x = 0 to 8", all(luhn_double(x) % 9 == 2 * x % 9 for x in range(9)))
    same("inverse of 2 modulo 9", inverse(2, 9), 5)
    same("gcd(6, 10)", math.gcd(6, 10), 2)
    same("solutions of 6x = 4 mod 10", solutions(6, 4, 10), [4, 9])
    same("inverse of 3 modulo 5", inverse(3, 5), 2)
    same("2 * 2 mod 5", 2 * 2 % 5, 4)
    same("6 * 4 and 6 * 9", [6 * 4, 6 * 9], [24, 54])
    same("solutions of 6x = 5 mod 10", solutions(6, 5, 10), [])
    same("x = 1 mod 4 and x = 2 mod 5", crt(1, 4, 2, 5), 17)
    same("17 mod 7 and 20 mod 7", [17 % 7, 20 % 7], [3, 6])
    same("inverse of 6 modulo 7", inverse(6, 7), 6)
    same("k from 6k = 2 mod 7", solutions(6, 2, 7), [5])
    same("17 + 20 * 5", 17 + 20 * 5, 117)
    same("117 mod 4, 5, 7", [117 % 4, 117 % 5, 117 % 7], [1, 2, 5])
    same("x = 1 mod 4, 2 mod 5, 5 mod 7", crt(17, 20, 5, 7), 117)
    same("117 in binary", binary(117), "1110101")
    same("5^(2^j) mod 19 for j = 0 to 6", squares(5, 7, 19), [5, 6, 17, 4, 16, 9, 5])
    same("5 * 9, 7 * 16, 17 * 17, 4 * 5 modulo 19",
         [5 * 9 % 19, 7 * 16 % 19, 17 * 17 % 19, 4 * 5 % 19], [7, 17, 4, 1])
    same("5^117 mod 19", power_mod(5, 117, 19), 1)
    same("117 = 18q + r", list(divmod(117, 18)), [6, 9])
    same("5^9 mod 19", power_mod(5, 9, 19), 1)
    same("smallest k with 5^k = 1 mod 19", order(5, 19), 9)


def check_exercises() -> None:
    section("Exercises")
    same("-23 by 6, floored q and r", floor_divmod(-23, 6), [-4, 1])
    same("-23 by 6, truncated q and r", trunc_divmod(-23, 6), [-3, -5])
    steps = euclid_steps(252, 198)
    same("Euclid on 252 and 198, quotients", [s[1] for s in steps], [1, 3, 1, 2])
    same("Euclid on 252 and 198, remainders", [s[3] for s in steps], [54, 36, 18, 0])
    same("gcd, s, t for 252 and 198", bezout(252, 198), [18, 4, -5])
    same("4 * 252 - 5 * 198", 4 * 252 - 5 * 198, 18)
    same("inverse of 7 modulo 26", inverse(7, 26), 15)
    same("7 * 15 and 4 * 26", [7 * 15, 4 * 26], [105, 104])
    same("gcd(13, 26) and the inverse of 13 modulo 26", [math.gcd(13, 26), inverse(13, 26)], [13, "none"])
    terms = isbn10_terms("0130426725")
    same("ISBN-10 sum without position 8", sum(terms) - terms[7], 100)
    same("-100 mod 11", -100 % 11, 10)
    same("inverse of 3 modulo 11", inverse(3, 11), 4)
    same("4 * 10 mod 11", 4 * 10 % 11, 7)
    same("100 = 12q + r", list(divmod(100, 12)), [8, 4])
    same("2^100 mod 13", power_mod(2, 100, 13), 3)
    same("phi(4), phi(25), phi(100)", [phi(4), phi(25), phi(100)], [2, 20, 40])
    same("2026 mod 40", 2026 % 40, 26)
    same("3, 3^2, 3^4, 3^8, 3^16 modulo 100", squares(3, 5, 100), [3, 9, 81, 61, 21])
    same("21 * 61 mod 100 and 81 * 9 mod 100", [21 * 61 % 100, 81 * 9 % 100], [81, 29])
    same("3^26 mod 100 and 3^2026 mod 100", [power_mod(3, 26, 100), power_mod(3, 2026, 100)], [29, 29])
    same("inverse of 7 modulo 11", inverse(7, 11), 8)
    same("k from 7k = 2 mod 11", solutions(7, 2, 11), [5])
    same("4 + 7 * 5", 4 + 7 * 5, 39)
    same("x = 4 mod 7 and x = 6 mod 11", crt(4, 7, 6, 11), 39)
    same("units modulo 9 and their inverses",
         [units(9), [inverse(a, 9) for a in units(9)]], [[1, 2, 4, 5, 7, 8], [1, 5, 7, 2, 4, 8]])
    same("3 * 3 mod 9", 3 * 3 % 9, 0)
    same("3 * 11 * 17", 3 * 11 * 17, 561)
    same("560 mod 2, 10, 16", [560 % 2, 560 % 10, 560 % 16], [0, 0, 0])
    same("phi(561)", phi(561), 320)
    require("a^560 = 1 mod 561 for all 320 a coprime to 561",
            all(power_mod(a, 560, 561) == 1 for a in units(561)))


def main() -> None:
    check_division()
    check_gcd()
    check_bezout()
    check_congruences()
    check_inverses()
    check_fermat_euler()
    check_crt()
    check_fields()
    check_examples()
    check_exercises()
    print("\nAll chapter 13 numbers reproduced.")


if __name__ == "__main__":
    main()
