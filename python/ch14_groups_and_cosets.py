"""Reproduce every number in chapter 14 of Error Correction Coding, "Groups and Cosets".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch14_groups_and_cosets.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in this chapter is a whole number, a bit string or a list of
them, so every comparison is exact. The small groups are checked
exhaustively: the group axioms on every triple of elements, the subgroup
tests on every subset of Z_n for n up to 12 and of the 3-bit strings, the
coset rules on every pair of elements, and the syndrome map of the code C_5
on all 32 words of five bits. A run that reaches the final line has
reproduced every number in the chapter.

The module also holds the group, coset and syndrome functions that
ch14_charts.py draws, so the figures and this check share one source.
"""

import sys
from itertools import combinations, permutations
from math import gcd


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


def show(x: int, n: int) -> str:
    return format(x, f"0{n}b")


def word(s: str) -> int:
    return int(s, 2)


def weight(x: int) -> int:
    return bin(x).count("1")


def as_list(xs, n: int = 0) -> str:
    return ", ".join(show(x, n) if n else str(x) for x in xs)


def as_set(xs, n: int = 0) -> str:
    return "{" + as_list(xs, n) + "}"


def as_row(xs, n: int = 0) -> str:
    return " ".join(show(x, n) if n else str(x) for x in xs)


# Groups -------------------------------------------------------------------------
# A finite group is held as the list of its elements and a function that
# combines two of them.

def z_add(n: int):
    """Z_n, the residues 0, 1, ..., n - 1 under addition mod n."""
    return list(range(n)), lambda a, b: (a + b) % n


def z_units(n: int):
    """The residues below n that share no factor with n, under multiplication mod n."""
    return [a for a in range(1, n) if gcd(a, n) == 1], lambda a, b: (a * b) % n


def bit_strings(n: int):
    """The 2^n strings of n bits under XOR."""
    return list(range(2 ** n)), lambda a, b: a ^ b


def identities(elems, op) -> list:
    return [e for e in elems if all(op(e, a) == a and op(a, e) == a for a in elems)]


def identity(elems, op) -> int:
    return identities(elems, op)[0]


def is_closed(elems, op) -> bool:
    members = set(elems)
    return all(op(a, b) in members for a in elems for b in elems)


def is_group(elems, op) -> bool:
    """Closure, associativity, an identity and an inverse for every element."""
    if not elems or not is_closed(elems, op):
        return False
    if any(op(op(a, b), c) != op(a, op(b, c)) for a in elems for b in elems for c in elems):
        return False
    found = identities(elems, op)
    if not found:
        return False
    e = found[0]
    return all(any(op(a, b) == e and op(b, a) == e for b in elems) for a in elems)


def inverses_of(elems, op, a) -> list:
    e = identity(elems, op)
    return [b for b in elems if op(a, b) == e and op(b, a) == e]


def inverse(elems, op, a) -> int:
    return inverses_of(elems, op, a)[0]


def table_rows(elems, op, n: int = 0) -> list:
    return [as_row([op(a, b) for b in elems], n) for a in elems]


def is_latin(elems, op) -> bool:
    target = sorted(elems)
    rows = all(sorted(op(a, b) for b in elems) == target for a in elems)
    cols = all(sorted(op(a, b) for a in elems) == target for b in elems)
    return rows and cols


def powers(elems, op, a) -> list:
    """a, a^2, a^3, ... up to and including the first power equal to the identity."""
    e = identity(elems, op)
    out = [a]
    while out[-1] != e:
        out.append(op(out[-1], a))
    return out


def order(elems, op, a) -> int:
    return len(powers(elems, op, a))


def power(elems, op, a, k: int) -> int:
    x = identity(elems, op)
    for _ in range(k):
        x = op(x, a)
    return x


def subsets(elems):
    for r in range(1, len(elems) + 1):
        yield from (list(c) for c in combinations(elems, r))


def subgroups(elems, op) -> list:
    return [h for h in subsets(elems) if is_closed(h, op) and is_group(h, op)]


def subgroup_test(elems, op, h) -> bool:
    members = set(h)
    return bool(h) and all(op(a, inverse(elems, op, b)) in members for a in h for b in h)


def coset(op, a, h) -> list:
    return sorted({op(a, x) for x in h})


def cosets(elems, op, h) -> list:
    out = []
    for a in elems:
        c = coset(op, a, h)
        if c not in out:
            out.append(c)
    return out


def is_hom(dom, op_dom, op_cod, f) -> bool:
    return all(f(op_dom(a, b)) == op_cod(f(a), f(b)) for a in dom for b in dom)


def kernel(dom, f, e_cod) -> list:
    return [a for a in dom if f(a) == e_cod]


def image(dom, f) -> list:
    return sorted({f(a) for a in dom})


# The code C_5 and its syndromes ---------------------------------------------------

C5 = [word("00000"), word("01101"), word("10110"), word("11011")]


def bit(x: int, i: int, n: int = 5) -> int:
    """Bit i of an n-bit word, counting from 1 at the left."""
    return (x >> (n - i)) & 1


def syndrome(x: int) -> int:
    """The three parity checks of C_5: x1+x2+x3, x1+x4 and x2+x5, each mod 2."""
    s1 = bit(x, 1) ^ bit(x, 2) ^ bit(x, 3)
    s2 = bit(x, 1) ^ bit(x, 4)
    s3 = bit(x, 2) ^ bit(x, 5)
    return (s1 << 2) | (s2 << 1) | s3


def leader(row) -> int:
    """Least weight, and among equals the smaller word read as a binary number."""
    return min(row, key=lambda x: (weight(x), x))


def standard_array(code: list, n: int) -> list:
    """Each coset as [leader ^ c for c in code], rows by leader weight, then leftmost bit first."""
    seen, rows = set(), []
    for x in range(2 ** n):
        if x in seen:
            continue
        row = sorted({x ^ c for c in code})
        seen.update(row)
        rows.append(leader(row))
    rows.sort(key=lambda x: (weight(x), -x))
    return [[ld ^ c for c in code] for ld in rows]


def primes_below(n: int) -> list:
    return [p for p in range(2, n) if all(p % q for q in range(2, int(p ** 0.5) + 1))]


def order_of_two(p: int) -> int:
    return order(*z_units(p), 2)


# The test groups ------------------------------------------------------------------

def test_groups() -> list:
    groups = [z_add(n) for n in range(1, 13)]
    groups += [bit_strings(n) for n in range(1, 4)]
    groups += [z_units(n) for n in (5, 7, 8, 9, 12)]
    return groups


# The chapter's numbers ------------------------------------------------------------

def check_tables() -> None:
    section("Three tables that behave alike")
    z4, add4 = z_add(4)
    b2, xor = bit_strings(2)
    u5, mul5 = z_units(5)
    for row, stated in zip(table_rows(z4, add4), ["0 1 2 3", "1 2 3 0", "2 3 0 1", "3 0 1 2"]):
        same(f"Z_4 under addition, row {row[0]}", row, stated)
    for row, stated in zip(table_rows(b2, xor, 2), ["00 01 10 11", "01 00 11 10", "10 11 00 01", "11 10 01 00"]):
        same(f"two-bit strings under XOR, row {row[:2]}", row, stated)
    for row, stated in zip(table_rows(u5, mul5), ["1 2 3 4", "2 4 1 3", "3 1 4 2", "4 3 2 1"]):
        same(f"nonzero residues mod 5 under multiplication, row {row[0]}", row, stated)
    rename = [pow(2, k, 5) for k in z4]
    same("2^k mod 5 for k = 0, 1, 2, 3", as_list(rename), "1, 2, 4, 3")
    require("2^(a + b mod 4) = 2^a 2^b mod 5 for all a, b", all(rename[add4(a, b)] == mul5(rename[a], rename[b]) for a in z4 for b in z4))
    same("1 + 1 in Z_4", add4(1, 1), 2)
    same("elements a of Z_4 with a + a = 0", as_list([a for a in z4 if add4(a, a) == 0]), "0, 2")
    same("two-bit strings x with x XOR x = 00", len([x for x in b2 if xor(x, x) == 0]), 4)
    require("every row and column of all three tables holds each element once", all(is_latin(*g) for g in [(z4, add4), (b2, xor), (u5, mul5)]))


def check_axioms() -> None:
    section("The group axioms")
    groups = test_groups()
    require("Z_n, n = 1..12, bit strings of 1..3 bits, Z_5*, Z_7*, Z_8*, Z_9*, Z_12* are groups", all(is_group(*g) for g in groups))
    require("each of them has exactly one identity", all(len(identities(*g)) == 1 for g in groups))
    require("each element of each of them has exactly one inverse", all(len(inverses_of(e, op, a)) == 1 for e, op in groups for a in e))
    require("the inverse of the inverse of a is a", all(inverse(e, op, inverse(e, op, a)) == a for e, op in groups for a in e))
    require("cancellation: a*b = a*c forces b = c", all(op(a, b) != op(a, c) for e, op in groups for a in e for b in e for c in e if b != c))
    require("every Cayley table is a Latin square", all(is_latin(*g) for g in groups))
    require("(a*b)^-1 = b^-1 * a^-1 in every test group", all(inverse(e, op, op(a, b)) == op(inverse(e, op, b), inverse(e, op, a)) for e, op in groups for a in e for b in e))
    same("2 x 3 mod 6", (2 * 3) % 6, 0)
    require("the nonzero residues mod 6 are not closed under multiplication", not is_closed(list(range(1, 6)), lambda a, b: (a * b) % 6))
    same("(5 - 3) - 1", (5 - 3) - 1, 1)
    same("5 - (3 - 1)", 5 - (3 - 1), 3)
    u5, mul5 = z_units(5)
    same("inverses of 1, 2, 3, 4 mod 5", as_list(inverse(u5, mul5, a) for a in u5), "1, 3, 2, 4")
    b3, xor = bit_strings(3)
    require("every 3-bit string is its own inverse under XOR", all(xor(x, x) == 0 for x in b3))


def check_subgroups() -> None:
    section("Subgroups")
    z12, add12 = z_add(12)
    evens = [0, 2, 4, 6, 8, 10]
    require("the even residues {0, 2, ..., 10} form a subgroup of Z_12", is_group(evens, add12))
    same("-4 in Z_12", inverse(z12, add12, 4), 8)
    require("{0, 4, 8} is a subgroup of Z_12", is_group([0, 4, 8], add12))
    require("{1, 5, 9} is not a subgroup of Z_12", not is_group([1, 5, 9], add12))
    same("1 + 5 in Z_12", add12(1, 5), 6)
    require("{0, 1, 2, 3, 4, 5} is not a subgroup of Z_12", not is_group([0, 1, 2, 3, 4, 5], add12))
    same("3 + 4 in Z_12", add12(3, 4), 7)
    ok_test, ok_finite = True, True
    for n in range(1, 13):
        elems, op = z_add(n)
        for h in subsets(elems):
            truth = is_group(h, op)
            ok_test &= subgroup_test(elems, op, h) == truth
            ok_finite &= is_closed(h, op) == truth
    b3, xor = bit_strings(3)
    for h in subsets(b3):
        truth = is_group(h, xor)
        ok_test &= subgroup_test(b3, xor, h) == truth
        ok_finite &= is_closed(h, xor) == truth
    require("subgroup test agrees with the axioms on every subset tried", ok_test)
    require("closure alone decides every finite subset tried", ok_finite)
    b5, xor5 = bit_strings(5)
    same("01101 XOR 10110", show(xor5(word("01101"), word("10110")), 5), "11011")
    same("01101 XOR 11011", show(xor5(word("01101"), word("11011")), 5), "10110")
    same("10110 XOR 11011", show(xor5(word("10110"), word("11011")), 5), "01101")
    require("C_5 is a subgroup of the 5-bit strings", is_group(C5, xor5))
    require("the even-weight words {000, 011, 101, 110} form a subgroup", is_group([0b000, 0b011, 0b101, 0b110], xor))


def check_cyclic() -> None:
    section("Cyclic groups and the order of an element")
    u7, mul7 = z_units(7)
    same("powers of 3 mod 7", as_list(powers(u7, mul7, 3)), "3, 2, 6, 4, 5, 1")
    same("powers of 2 mod 7", as_list(powers(u7, mul7, 2)), "2, 4, 1")
    same("order of 6 mod 7", order(u7, mul7, 6), 2)
    z12, add12 = z_add(12)
    same("multiples of 9 in Z_12", as_list(powers(z12, add12, 9)), "9, 6, 3, 0")
    same("orders of 0, 1, ..., 11 in Z_12", as_list(order(z12, add12, a) for a in z12), "1, 12, 6, 4, 3, 12, 2, 12, 3, 4, 6, 12")
    same("generators of Z_12", as_list(a for a in z12 if order(z12, add12, a) == 12), "1, 5, 7, 11")
    require("the order of a in Z_n is n / gcd(a, n), n = 1..40", all(order(*z_add(n), a) == n // gcd(a, n) for n in range(1, 41) for a in range(n)))
    groups = test_groups()
    require("<a> has exactly ord(a) distinct elements in every test group", all(len(set(powers(e, op, a))) == order(e, op, a) for e, op in groups for a in e))
    require("a^k = e exactly when ord(a) divides k, k = 0..36", all((power(e, op, a, k) == identity(e, op)) == (k % order(e, op, a) == 0) for e, op in groups for a in e for k in range(37)))
    ok = True
    for n in range(1, 13):
        elems, op = z_add(n)
        found = sorted(tuple(h) for h in subgroups(elems, op))
        expected = sorted(tuple(sorted(set(powers(elems, op, (n // d) % n)))) for d in range(1, n + 1) if n % d == 0)
        ok &= found == expected
    require("the subgroups of Z_n are the <n/d> for d dividing n, n = 1..12", ok)
    found = subgroups(z12, add12)
    same("number of subgroups of Z_12", len(found), 6)
    stated = {1: "{0}", 2: "{0, 6}", 3: "{0, 4, 8}", 4: "{0, 3, 6, 9}", 6: "{0, 2, 4, 6, 8, 10}"}
    for h in sorted(found, key=len):
        if len(h) in stated:
            same(f"the subgroup of Z_12 of order {len(h)}", as_set(h), stated[len(h)])
        else:
            same("the subgroup of Z_12 of order 12", len(h), 12)


def check_cosets() -> None:
    section("Cosets")
    z12, add12 = z_add(12)
    found = cosets(z12, add12, [0, 4, 8])
    for c, stated in zip(found, ["{0, 4, 8}", "{1, 5, 9}", "{2, 6, 10}", "{3, 7, 11}"]):
        same(f"coset {c[0]} + {{0, 4, 8}} in Z_12", as_set(c), stated)
    same("number of cosets of {0, 4, 8}", len(found), 4)
    same("4 + {0, 4, 8}", as_set(coset(add12, 4, [0, 4, 8])), "{0, 4, 8}")
    found = cosets(z12, add12, [0, 3, 6, 9])
    for c, stated in zip(found, ["{0, 3, 6, 9}", "{1, 4, 7, 10}", "{2, 5, 8, 11}"]):
        same(f"coset {c[0]} + {{0, 3, 6, 9}} in Z_12", as_set(c), stated)
    b3, xor = bit_strings(3)
    r3 = [0b000, 0b111]
    for a, stated in zip([0b000, 0b001, 0b010, 0b100], ["{000, 111}", "{001, 110}", "{010, 101}", "{011, 100}"]):
        same(f"{show(a, 3)} XOR R_3", as_set(coset(xor, a, r3), 3), stated)
    same("number of cosets of R_3", len(cosets(b3, xor, r3)), 4)
    cases = [(z12, add12, h) for h in subgroups(z12, add12)]
    cases += [(b3, xor, h) for h in subgroups(b3, xor)]
    u7, mul7 = z_units(7)
    cases += [(u7, mul7, h) for h in subgroups(u7, mul7)]
    require("aH = bH exactly when a^-1 b lies in H, every case", all((coset(op, a, h) == coset(op, b, h)) == (op(inverse(g, op, a), b) in h) for g, op, h in cases for a in g for b in g))
    require("two cosets are equal or share no element, every case", all(coset(op, a, h) == coset(op, b, h) or not set(coset(op, a, h)) & set(coset(op, b, h)) for g, op, h in cases for a in g for b in g))
    require("every coset has as many elements as H, every case", all(len(coset(op, a, h)) == len(h) for g, op, h in cases for a in g))
    require("a lies in its own coset aH, every case", all(a in coset(op, a, h) for g, op, h in cases for a in g))


def check_lagrange() -> None:
    section("Lagrange's theorem")
    z12, add12 = z_add(12)
    found = sorted(subgroups(z12, add12), key=len)
    same("orders of the subgroups of Z_12", as_list(len(h) for h in found), "1, 2, 3, 4, 6, 12")
    same("number of cosets of each, 12 / |H|", as_list(len(cosets(z12, add12, h)) for h in found), "12, 6, 4, 3, 2, 1")
    u7, mul7 = z_units(7)
    same("orders of 1, 2, 3, 4, 5, 6 mod 7", as_list(order(u7, mul7, a) for a in u7), "1, 3, 6, 3, 6, 2")
    same("powers of 5 mod 7", as_list(powers(u7, mul7, 5)), "5, 4, 6, 2, 3, 1")
    same("powers of 1 mod 7", as_list(powers(u7, mul7, 1)), "1")
    same("powers of 4 mod 7", as_list(powers(u7, mul7, 4)), "4, 2, 1")
    same("powers of 6 mod 7", as_list(powers(u7, mul7, 6)), "6, 1")
    same("coset 3{1, 2, 4} mod 7", as_set(coset(mul7, 3, powers(u7, mul7, 2))), "{3, 5, 6}")
    same("elements of order 4 mod 7", len([a for a in u7 if order(u7, mul7, a) == 4]), 0)
    groups = test_groups()
    require("every element order divides the group order, every test group", all(len(e) % order(e, op, a) == 0 for e, op in groups for a in e))
    require("a^|G| = e in every test group", all(power(e, op, a, len(e)) == identity(e, op) for e, op in groups for a in e))
    require("a^(p - 1) = 1 mod p for every prime p < 60 and 0 < a < p", all(pow(a, p - 1, p) == 1 for p in primes_below(60) for a in range(1, p)))
    require("a^phi(n) = 1 mod n for every n from 2 to 40 and a coprime to n", all(pow(a, len(z_units(n)[0]), n) == 1 for n in range(2, 41) for a in z_units(n)[0]))
    ok = True
    for p in primes_below(14):
        elems, op = z_add(p)
        ok &= [len(h) for h in subgroups(elems, op)] == [1, p]
        ok &= all(order(elems, op, a) == p for a in elems if a)
    require("Z_p has only {0} and Z_p, and every a != 0 generates it, p < 14", ok)
    b3, xor = bit_strings(3)
    same("sizes of the subgroups of the 3-bit strings", as_list(sorted({len(h) for h in subgroups(b3, xor)})), "1, 2, 4, 8")
    same("64 mod 12, so no 12-word linear code of length 6", 64 % 12, 4)
    same("32 mod 6, so no 6-word linear code of length 5", 32 % 6, 2)
    same("12 mod 5, so no subgroup of order 5 when |G| = 12", 12 % 5, 2)
    same("cosets of a 16-word linear code of length 7", 2 ** 7 // 16, 8)
    odd = [p for p in primes_below(100) if p > 2]
    orders = {p: order_of_two(p) for p in odd}
    same("odd primes below 100", len(odd), 24)
    require("the order of 2 mod p divides p - 1 for each of them", all((p - 1) % orders[p] == 0 for p in odd))
    full = [p for p in odd if orders[p] == p - 1]
    half = [p for p in odd if orders[p] == (p - 1) // 2]
    rest = [p for p in odd if p not in full and p not in half]
    same("primes where 2 has order p - 1", as_list(full), "3, 5, 11, 13, 19, 29, 37, 53, 59, 61, 67, 83")
    same("how many", len(full), 12)
    same("primes where 2 has order (p - 1) / 2", as_list(half), "7, 17, 23, 41, 47, 71, 79, 97")
    same("how many", len(half), 8)
    same("the other primes and the order of 2", ", ".join(f"{p}: {orders[p]}" for p in rest), "31: 5, 43: 14, 73: 9, 89: 11")


def check_homomorphisms() -> None:
    section("Homomorphisms and kernels")
    z12, add12 = z_add(12)
    z4, add4 = z_add(4)
    mod4 = lambda x: x % 4
    require("x -> x mod 4 from Z_12 to Z_4 is a homomorphism", is_hom(z12, add12, add4, mod4))
    same("its kernel", as_set(kernel(z12, mod4, 0)), "{0, 4, 8}")
    fibers = [as_set([x for x in z12 if mod4(x) == y]) for y in z4]
    same("the words sent to 1", fibers[1], "{1, 5, 9}")
    require("the words sent to each value are a coset of the kernel", all(f == as_set(c) for f, c in zip(fibers, cosets(z12, add12, [0, 4, 8]))))
    b3, xor = bit_strings(3)
    parity = lambda x: weight(x) % 2
    require("parity from 3-bit strings to Z_2 is a homomorphism", is_hom(b3, xor, lambda a, b: (a + b) % 2, parity))
    same("kernel of parity on 3-bit strings", as_set(kernel(b3, parity, 0), 3), "{000, 011, 101, 110}")
    same("the 3-bit strings of odd parity", as_set([x for x in b3 if parity(x)], 3), "{001, 010, 100, 111}")
    maps = []
    for k in range(12):
        f = (lambda k: lambda x: (k * x) % 12)(k)
        maps.append(f)
    require("every x -> kx mod 12 is a homomorphism of Z_12", all(is_hom(z12, add12, add12, f) for f in maps))
    require("f(0) = 0 and f(-a) = -f(a) for each of them", all(f(0) == 0 and all(f(inverse(z12, add12, a)) == inverse(z12, add12, f(a)) for a in z12) for f in maps))
    require("kernel and image are subgroups for each of them", all(is_group(kernel(z12, f, 0), add12) and is_group(image(z12, f), add12) for f in maps))
    require("f(a) = f(b) exactly when b - a is in the kernel, each map", all((f(a) == f(b)) == (add12(b, inverse(z12, add12, a)) in kernel(z12, f, 0)) for f in maps for a in z12 for b in z12))
    require("one-to-one exactly when the kernel is {0}, each map", all((len(image(z12, f)) == 12) == (kernel(z12, f, 0) == [0]) for f in maps))
    u5, mul5 = z_units(5)
    two_to = lambda k: pow(2, k, 5)
    require("k -> 2^k mod 5 is a homomorphism from Z_4", is_hom(z4, add4, mul5, two_to))
    same("2^4 mod 5", pow(2, 4, 5), 1)
    require("k -> 2^k mod 5 is one-to-one and onto", sorted(two_to(k) for k in z4) == u5)
    b2, xor2 = bit_strings(2)
    homs = [p for p in permutations(b2) if is_hom(z4, add4, xor2, lambda k, p=p: p[k])]
    same("bijections from Z_4 to 2-bit strings that are homomorphisms", len(homs), 0)
    same("bijections tried", len(list(permutations(b2))), 24)


def check_quotients() -> None:
    section("Quotient groups")
    z12, add12 = z_add(12)
    h = [0, 4, 8]
    reps = [c[0] for c in cosets(z12, add12, h)]
    name = lambda x: coset(add12, x, h)[0]
    rows = [as_row([name(add12(a, b)) for b in reps]) for a in reps]
    for row, stated in zip(rows, ["0 1 2 3", "1 2 3 0", "2 3 0 1", "3 0 1 2"]):
        same(f"cosets of {{0, 4, 8}} added, row {row[0]}", row, stated)
    same("1 + 2 in Z_12", add12(1, 2), 3)
    same("5 + 10 in Z_12", add12(5, 10), 3)
    same("{1, 5, 9} + {2, 6, 10}", as_set(coset(add12, 3, h)), "{3, 7, 11}")
    b3, xor = bit_strings(3)
    cases = [(z12, add12, k) for k in subgroups(z12, add12)] + [(b3, xor, k) for k in subgroups(b3, xor)]
    ok = True
    for g, op, k in cases:
        for a in g:
            for a2 in coset(op, a, k):
                for b in g:
                    for b2 in coset(op, b, k):
                        ok &= coset(op, op(a, b), k) == coset(op, op(a2, b2), k)
    require("the coset sum does not depend on the representatives", ok)
    require("the cosets form a group under that sum, every case", all(is_group([c[0] for c in cosets(g, op, k)], lambda a, b, op=op, k=k: coset(op, op(a, b), k)[0]) for g, op, k in cases))
    require("|G/ker f| = |image f| for every x -> kx mod 12", all(12 // len(kernel(z12, f, 0)) == len(image(z12, f)) for f in [(lambda k: lambda x: (k * x) % 12)(k) for k in range(12)]))
    parity = lambda x: weight(x) % 2
    same("cosets of the even-weight words in 3-bit strings", len(cosets(b3, xor, kernel(b3, parity, 0))), 2)


def check_syndromes() -> None:
    section("Cosets and syndromes")
    b5, xor = bit_strings(5)
    require("all three parity checks are 0 on every codeword of C_5", all(syndrome(c) == 0 for c in C5))
    require("the syndrome of x XOR y is the XOR of their syndromes", is_hom(b5, xor, xor, syndrome))
    same("kernel of the syndrome map", as_set(kernel(b5, syndrome, 0), 5), "{00000, 01101, 10110, 11011}")
    same("number of syndromes reached", len(image(b5, syndrome)), 8)
    same("number of cosets of C_5, 32 / 4", len(cosets(b5, xor, C5)), 8)
    require("same syndrome exactly when x XOR y lies in C_5", all((syndrome(x) == syndrome(y)) == ((x ^ y) in C5) for x in b5 for y in b5))
    rows = standard_array(C5, 5)
    stated = [("00000 01101 10110 11011", "000"), ("10000 11101 00110 01011", "110"),
              ("01000 00101 11110 10011", "101"), ("00100 01001 10010 11111", "100"),
              ("00010 01111 10100 11001", "010"), ("00001 01100 10111 11010", "001"),
              ("01010 00111 11100 10001", "111"), ("00011 01110 10101 11000", "011")]
    for row, (words, syn) in zip(rows, stated):
        same(f"coset with leader {show(row[0], 5)}", as_row(row, 5), words)
        same(f"its syndrome", show(syndrome(row[0]), 3), syn)
        require("every word of that coset has that syndrome", all(syndrome(x) == syndrome(row[0]) for x in row))
    same("leader weights", as_list(weight(r[0]) for r in rows), "0, 1, 1, 1, 1, 1, 2, 2")
    same("8 cosets of 4 words", sum(len(r) for r in rows), 32)
    same("syndrome of 11111", show(syndrome(word("11111")), 3), "100")
    same("11111 XOR its leader 00100", show(word("11111") ^ word("00100"), 5), "11011")
    ties = [as_list(sorted((x for x in r if weight(x) == 2), reverse=True), 5) for r in rows[6:]]
    same("least-weight words in the coset of syndrome 111", ties[0], "10001, 01010")
    same("least-weight words in the coset of syndrome 011", ties[1], "11000, 00011")
    doubles = [x for x in b5 if weight(x) == 2]
    in_low = [x for x in doubles if weight(leader(next(r for r in rows if x in r))) == 1]
    same("weight-2 words that lie in a coset with a weight-1 leader", len(in_low), 6)
    same("weight-2 words in all", len(doubles), 10)
    require("the leader of r's coset is r XOR a nearest codeword, all r", all(weight(leader(next(r for r in rows if x in r))) == min(weight(x ^ c) for c in C5) for x in b5))


def check_examples() -> None:
    section("Worked examples")
    u8, mul8 = z_units(8)
    for row, stated in zip(table_rows(u8, mul8), ["1 3 5 7", "3 1 7 5", "5 7 1 3", "7 5 3 1"]):
        same(f"residues 1, 3, 5, 7 mod 8 multiplied, row {row[0]}", row, stated)
    require("Z_8* is a group", is_group(u8, mul8))
    same("squares of 1, 3, 5, 7 mod 8", as_list(mul8(a, a) for a in u8), "1, 1, 1, 1")
    b2, xor2 = bit_strings(2)
    rename = dict(zip(u8, b2))
    require("1, 3, 5, 7 -> 00, 01, 10, 11 is an isomorphism onto 2-bit strings", is_hom(u8, mul8, xor2, lambda a: rename[a]))
    u9, mul9 = z_units(9)
    same("residues mod 9 that share no factor with 9", as_list(u9), "1, 2, 4, 5, 7, 8")
    same("powers of 2 mod 9", as_list(powers(u9, mul9, 2)), "2, 4, 8, 7, 5, 1")
    same("powers of 4 mod 9", as_list(powers(u9, mul9, 4)), "4, 7, 1")
    same("powers of 8 mod 9", as_list(powers(u9, mul9, 8)), "8, 1")
    same("orders of 1, 2, 4, 5, 7, 8 mod 9", as_list(order(u9, mul9, a) for a in u9), "1, 6, 3, 6, 3, 2")
    same("subgroups of Z_9*", ", ".join(as_set(h) for h in sorted(subgroups(u9, mul9), key=len)), "{1}, {1, 8}, {1, 4, 7}, {1, 2, 4, 5, 7, 8}")
    same("coset 2{1, 4, 7} mod 9", as_set(coset(mul9, 2, [1, 4, 7])), "{2, 5, 8}")
    r = word("10011")
    same("syndrome of 10011", show(syndrome(r), 3), "101")
    ld = next(row[0] for row in standard_array(C5, 5) if syndrome(row[0]) == syndrome(r))
    same("its coset leader", show(ld, 5), "01000")
    same("10011 XOR 01000", show(r ^ ld, 5), "11011")
    same("distances from 10011 to the four codewords", as_list(weight(r ^ c) for c in C5), "3, 4, 2, 1")


def check_exercises() -> None:
    section("Exercises")
    u5, mul5 = z_units(5)
    same("order of 2 mod 5", order(u5, mul5, 2), 4)
    same("order of 4 mod 5", order(u5, mul5, 4), 2)
    z12, add12 = z_add(12)
    for h, stated in [([0, 6], True), ([0, 5], False), ([0, 3, 6, 9], True), ([0, 2, 4, 8], False)]:
        same(f"{as_set(h)} is a subgroup of Z_12", is_group(h, add12), stated)
    same("5 + 5 in Z_12", add12(5, 5), 10)
    same("2 + 4 in Z_12", add12(2, 4), 6)
    same("order of 8 in Z_12", order(z12, add12, 8), 3)
    same("order of 10 in Z_12", order(z12, add12, 10), 6)
    same("order of 6 in Z_20", order(*z_add(20), 6), 10)
    b5, xor = bit_strings(5)
    r5 = [0, word("11111")]
    found = cosets(b5, xor, r5)
    same("number of cosets of R_5", len(found), 16)
    same("leaders of weight 0, 1, 2", as_list(sum(1 for c in found if weight(leader(c)) == w) for w in range(3)), "1, 5, 10")
    require("each coset of R_5 holds exactly one word of weight at most 2", all(sum(1 for x in c if weight(x) <= 2) == 1 for c in found))
    triple = lambda x: (3 * x) % 12
    same("kernel of x -> 3x mod 12", as_set(kernel(z12, triple, 0)), "{0, 4, 8}")
    same("image of x -> 3x mod 12", as_set(image(z12, triple)), "{0, 3, 6, 9}")
    r = word("10101")
    same("syndrome of 10101", show(syndrome(r), 3), "011")
    row = next(row for row in standard_array(C5, 5) if r in row)
    same("the coset of 10101", as_set(sorted(row), 5), "{00011, 01110, 10101, 11000}")
    same("distances from 10101 to the four codewords", as_list(weight(r ^ c) for c in C5), "3, 2, 2, 3")


def main() -> None:
    check_tables()
    check_axioms()
    check_subgroups()
    check_cyclic()
    check_cosets()
    check_lagrange()
    check_homomorphisms()
    check_quotients()
    check_syndromes()
    check_examples()
    check_exercises()
    print("\nAll chapter 14 numbers reproduced.")


if __name__ == "__main__":
    main()
