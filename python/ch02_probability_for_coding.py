"""Reproduce every number in chapter 2 of Error Correction Coding, "Probability for Coding".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch02_probability_for_coding.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. A run that reaches
the final line has reproduced every number in the chapter.
"""

import random
import sys
from fractions import Fraction
from itertools import product
from math import ceil, comb, erfc, exp, log, pi, sqrt


# Events on three bits ---------------------------------------------------------

def pattern_prob(pattern: str, p):
    """Probability that a BSC(p) turns 000... into `pattern`, flips being independent."""
    w = pattern.count("1")
    return p**w * (1 - p) ** (len(pattern) - w)


def prob(event, p):
    """Probability of a set of received patterns, by additivity over its outcomes."""
    return sum((pattern_prob(x, p) for x in event), Fraction(0))


def patterns(n: int) -> list:
    return ["".join(bits) for bits in product("01", repeat=n)]


# The binomial distribution ------------------------------------------------------

def binom_pmf(n: int, k: int, p):
    """P(X = k) for X the number of flips among n bits of a BSC(p)."""
    return comb(n, k) * p**k * (1 - p) ** (n - k)


def binom_tail(n: int, t: int, p: Fraction) -> Fraction:
    """Exact P(X >= t), computed as one minus the head with integer arithmetic."""
    num, den = p.numerator, p.denominator
    head = sum(comb(n, k) * num**k * (den - num) ** (n - k) for k in range(t))
    return Fraction(den**n - head, den**n)


def rep_error(n: int, p: Fraction) -> Fraction:
    """Decoded error of the repetition code R_n: more than half of n copies flipped."""
    return binom_tail(n, (n + 1) // 2, p)


# Bounds on a tail ---------------------------------------------------------------

def markov(mean: float, a: float) -> float:
    return mean / a


def chebyshev(var: float, d: float) -> float:
    return var / d**2


def kl(a: float, p: float) -> float:
    """D(a||p) = a ln(a/p) + (1-a) ln((1-a)/(1-p)), in nats."""
    return a * log(a / p) + (1 - a) * log((1 - a) / (1 - p))


def chernoff(n: int, a: float, p: float) -> float:
    """The Chernoff bound exp(-n D(a||p)) on P(X >= an), valid for p < a < 1."""
    return exp(-n * kl(a, p))


# The Gaussian and the Q-function ---------------------------------------------

def phi(x: float) -> float:
    """The standard Gaussian density."""
    return exp(-x * x / 2) / sqrt(2 * pi)


def Q(x: float) -> float:
    """P(Z > x) for a standard Gaussian Z, through the complementary error function."""
    return 0.5 * erfc(x / sqrt(2))


def q_lower(x: float) -> float:
    return x * phi(x) / (1 + x * x)


def q_upper(x: float) -> float:
    return phi(x) / x


def q_half(x: float) -> float:
    return 0.5 * exp(-x * x / 2)


# The law of large numbers, simulated -------------------------------------------

def running_average(seed: int, n: int, p: float) -> list:
    """X_k / k for k = 1..n, where X_k counts flips among the first k simulated bits.

    random.Random(seed).random() gives the same sequence on every Python version
    for the same integer seed, so the chart drawn from it never changes.
    """
    rng = random.Random(seed)
    flips, out = 0, []
    for k in range(1, n + 1):
        flips += rng.random() < p
        out.append(flips / k)
    return out


# Checks ----------------------------------------------------------------------

def half_unit(stated: str) -> float:
    """Half a unit in the last digit of a number written as the chapter writes it."""
    text = stated.replace(",", "")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    """Print a computed value beside the chapter's, and stop at the first disagreement."""
    value = float(stated.replace(",", ""))
    print(f"  {label:<60} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def main() -> None:
    p = Fraction(1, 10)

    print("Events on three bits")
    space = patterns(3)
    check("sum of all eight pattern probabilities", prob(space, p), "1")
    check("P(000)", pattern_prob("000", p), "0.729")
    check("P(100), one flip", pattern_prob("100", p), "0.081")
    check("P(110), two flips", pattern_prob("110", p), "0.009")
    check("P(111)", pattern_prob("111", p), "0.001")
    A = {x for x in space if x[0] == "1"}
    B = {x for x in space if x.count("1") >= 2}
    C = {x for x in space if x[1] == "1"}
    check("P(A), first bit flipped", prob(A, p), "0.1")
    check("P(B), decoder fails", prob(B, p), "0.028")
    check("P(A and B)", prob(A & B, p), "0.019")
    check("P(A or B)", prob(A | B, p), "0.109")
    check("P(A) + P(B) - P(A and B)", prob(A, p) + prob(B, p) - prob(A & B, p), "0.109")
    check("P(not B)", 1 - prob(B, p), "0.972")
    check("P(B | A)", prob(A & B, p) / prob(A, p), "0.19")
    check("1 - 0.9^2", 1 - (1 - p) ** 2, "0.19")
    notA = set(space) - A
    check("P(B | not A)", prob(B & notA, p) / prob(notA, p), "0.01")
    check("total probability 0.19(0.1) + 0.01(0.9)", Fraction(19, 100) * p + Fraction(1, 100) * (1 - p), "0.028")
    bayes = prob(A & B, p) / prob(B, p)
    assert bayes == Fraction(19, 28)
    check("P(A | B) = 19/28", bayes, "0.6786")
    check("P(A)P(B)", prob(A, p) * prob(B, p), "0.0028")
    check("P(A and B) / (P(A)P(B))", prob(A & B, p) / (prob(A, p) * prob(B, p)), "6.79")
    check("P(A and C), first and second bits flipped", prob(A & C, p), "0.01")
    assert prob(A & C, p) == prob(A, p) * prob(C, p)

    print("The binomial distribution")
    for k, stated in enumerate(["0.729", "0.243", "0.027", "0.001"]):
        check(f"n = 3, p = 0.1: P(X = {k})", binom_pmf(3, k, p), stated)
    q = Fraction(1, 100)
    eight = ["0.922745", "0.074565", "0.0026361", "5.3255e-5", "6.7242e-7"]
    for k, stated in enumerate(eight):
        check(f"n = 8, p = 0.01: P(X = {k})", binom_pmf(8, k, q), stated)
    check("n = 8, p = 0.01: P(X >= 5)", binom_tail(8, 5, q), "5.4612e-9")
    check("n = 8, p = 0.01: P(X >= 1)", binom_tail(8, 1, q), "0.077255")
    odd = sum(binom_pmf(8, k, q) for k in range(1, 9, 2))
    check("n = 8, p = 0.01: P(X odd)", odd, "0.0746185")
    check("(1 - 0.98^8) / 2", (1 - (1 - 2 * q) ** 8) / 2, "0.0746185")
    even2 = sum(binom_pmf(8, k, q) for k in range(2, 9, 2))
    check("n = 8, p = 0.01: P(X even and at least 2)", even2, "0.0026368")
    twenty = ["0.1216", "0.2702", "0.2852", "0.1901", "0.0898", "0.0319", "0.0089", "0.0020", "0.0004"]
    for k, stated in enumerate(twenty):
        check(f"n = 20, p = 0.1: P(X = {k})", binom_pmf(20, k, p), stated)
    check("n = 20, p = 0.1: (n + 1)p", 21 * p, "2.1")
    check("n = 20, p = 0.1: P(X >= 5)", binom_tail(20, 5, p), "0.0432")
    ratio = [binom_pmf(20, k + 1, p) / binom_pmf(20, k, p) for k in range(4)]
    check("n = 20: P(1)/P(0) = (20/1)(1/9)", ratio[0], "2.22")
    check("n = 20: P(2)/P(1) = (19/2)(1/9)", ratio[1], "1.06")
    check("n = 20: P(3)/P(2) = (18/3)(1/9)", ratio[2], "0.667")

    print("Expectation and variance")
    mean3 = sum(k * binom_pmf(3, k, p) for k in range(4))
    sq3 = sum(k * k * binom_pmf(3, k, p) for k in range(4))
    check("n = 3: E[X]", mean3, "0.3")
    check("n = 3: E[X^2]", sq3, "0.36")
    check("n = 3: Var X = E[X^2] - E[X]^2", sq3 - mean3**2, "0.27")
    check("n = 3: np(1 - p)", 3 * p * (1 - p), "0.27")
    direct = sum((k - mean3) ** 2 * binom_pmf(3, k, p) for k in range(4))
    check("n = 3: E[(X - 0.3)^2] term by term", direct, "0.27")
    for k, stated in enumerate(["0.09", "0.49", "2.89", "7.29"]):
        check(f"n = 3: (k - 0.3)^2 at k = {k}", (k - mean3) ** 2, stated)
    check("n = 3: standard deviation", sqrt(0.27), "0.5196")
    ber = 1e-12
    hour = 8e9 * 3600
    check("bits in an hour at 8 Gb/s", hour, "2.88e13")
    check("errors per second at 8 Gb/s and p = 1e-12", 8e9 * ber, "0.008")
    check("hour: mean errors", hour * ber, "28.8")
    check("hour: standard deviation", sqrt(hour * ber * (1 - ber)), "5.367")
    day = hour * 24
    check("bits in a day at 8 Gb/s", day, "6.912e14")
    check("day: mean errors", day * ber, "691.2")
    check("day: standard deviation", sqrt(day * ber * (1 - ber)), "26.29")
    check("terabyte: mean errors, 8e12 bits", 8e12 * ber, "8")
    check("terabyte: standard deviation", sqrt(8e12 * ber * (1 - ber)), "2.828")
    check("hour: sd / mean", sqrt(28.8) / 28.8, "0.186")
    check("day: sd / mean", sqrt(691.2) / 691.2, "0.038")
    check("terabyte: sd / mean", sqrt(8) / 8, "0.354")

    print("Markov, Chebyshev and the law of large numbers")
    exact20 = binom_tail(1000, 20, q)
    check("n = 1000, p = 0.01: mean", 1000 * 0.01, "10")
    check("n = 1000, p = 0.01: variance", 1000 * 0.01 * 0.99, "9.9")
    check("n = 1000: Markov bound on P(X >= 20)", markov(10, 20), "0.5")
    check("n = 1000: Chebyshev bound on P(X >= 20)", chebyshev(9.9, 10), "0.099")
    check("n = 1000: exact P(X >= 20)", exact20, "0.0032884")
    sd_day = sqrt(691.2 * (1 - ber))
    check("day: four standard deviations", 4 * sd_day, "105.2")
    check("day: mean - 4 sd", 691.2 - 4 * sd_day, "586.0")
    check("day: mean + 4 sd", 691.2 + 4 * sd_day, "796.4")
    check("Chebyshev at four sd: 1/16", Fraction(1, 16), "0.0625")
    eps, delta, pm = 1e-4, 0.05, 1e-3
    check("BER p = 1e-3 to 1e-4 at 95%: p(1-p)/(eps^2 delta)", pm * (1 - pm) / (eps**2 * delta), "1,998,000")
    check("n = 1e6, p = 1e-3: sd of the measured fraction", sqrt(pm * (1 - pm) / 1e6), "3.16e-5")
    check("n = 1e6, p = 1e-3: two sd of the fraction", 2 * sqrt(pm * (1 - pm) / 1e6), "6.3e-5")
    for seed, stated in ((1, "0.1"), (2, "0.1"), (3, "0.1")):
        final = running_average(seed, 10000, 0.1)[-1]
        assert abs(final - 0.1) <= 2 * sqrt(0.1 * 0.9 / 10000), (seed, final)
        check(f"simulated run {seed}: X/n at n = 10000, to one digit", round(final, 1), stated)
    check("band half-width 2 sqrt(p(1-p)/n) at n = 10000", 2 * sqrt(0.1 * 0.9 / 10000), "0.006")

    print("The Gaussian and the Q-function")
    for x, stated in ((1, "0.1587"), (2, "0.02275"), (3, "0.00135"), (4, "3.167e-5"), (5, "2.867e-7"), (6, "9.866e-10")):
        check(f"Q({x})", Q(x), stated)
    check("peak density of N(1, 0.25): 1/(0.5 sqrt(2 pi))", 1 / (0.5 * sqrt(2 * pi)), "0.798")
    check("phi(0) = 1/sqrt(2 pi)", phi(0), "0.3989")
    check("sigma = 0.4: Q(2.5)", Q(2.5), "0.00621")
    check("sigma = 0.25: Q(4)", Q(4), "3.167e-5")
    check("Q(2) / Q(4)", Q(2) / Q(4), "718")
    table = {
        1: ("1.210e-1", "1.587e-1", "2.420e-1", "3.033e-1"),
        2: ("2.160e-2", "2.275e-2", "2.700e-2", "6.767e-2"),
        3: ("1.330e-3", "1.350e-3", "1.477e-3", "5.554e-3"),
        4: ("3.149e-5", "3.167e-5", "3.346e-5", "1.677e-4"),
        5: ("2.859e-7", "2.867e-7", "2.973e-7", "1.863e-6"),
        6: ("9.853e-10", "9.866e-10", "1.013e-9", "7.615e-9"),
    }
    for x, (lo, mid, up, half) in table.items():
        assert q_lower(x) < Q(x) < q_upper(x) and Q(x) <= q_half(x)
        check(f"x = {x}: lower bound x phi/(1 + x^2)", q_lower(x), lo)
        check(f"x = {x}: Q(x)", Q(x), mid)
        check(f"x = {x}: upper bound phi/x", q_upper(x), up)
        check(f"x = {x}: bound exp(-x^2/2)/2", q_half(x), half)
    check("x = 6: upper / Q", q_upper(6) / Q(6), "1.026")
    check("x = 6: Q / lower", Q(6) / q_lower(6), "1.0013")
    n_mean, n_sd = 10.0, sqrt(9.9)
    check("Gaussian stand-in, n = 1000: (20 - 10)/sqrt(9.9)", (20 - n_mean) / n_sd, "3.178")
    check("Gaussian stand-in, n = 1000: Q(3.178)", Q((20 - n_mean) / n_sd), "7.41e-4")
    z2 = (200 - 100) / sqrt(99)
    check("Gaussian stand-in, n = 10000: (200 - 100)/sqrt(99)", z2, "10.05")
    check("Gaussian stand-in, n = 10000: Q(10.05)", Q(z2), "4.6e-24")
    check("n = 10000, p = 0.01: exact P(X >= 200)", binom_tail(10000, 200, q), "5.630e-19")

    print("The union bound")
    m = Fraction(1, 10**6)
    check("72 bits, q = 1e-6: union bound 72q", 72 * m, "7.2e-5")
    check("72 bits: exact P(X >= 1)", binom_tail(72, 1, m), "7.19974e-5")
    check("72 bits: C(72,2)", comb(72, 2), "2556")
    check("72 bits: union bound C(72,2) q^2", comb(72, 2) * m**2, "2.556e-9")
    check("72 bits: exact P(X >= 2)", binom_tail(72, 2, m), "2.55588e-9")
    check("R3 at p = 0.1: 3p^2", 3 * p**2, "0.03")
    check("R3 at p = 0.1: exact", rep_error(3, p), "0.028")
    check("R3 at p = 0.1: overcount 2p^3", 2 * p**3, "0.002")
    check("R5 at p = 0.01: 10p^3", 10 * q**3, "1e-5")
    check("R5 at p = 0.01: exact", rep_error(5, q), "9.8506e-6")
    check("DIMM: 16 x 0.08", 16 * 0.08, "1.28")
    check("DIMM: 1 - 0.92^16", 1 - 0.92**16, "0.7366")

    print("The Chernoff bound")
    a = 0.02
    es = a * (1 - 0.01) / (0.01 * (1 - a))
    check("n = 1000, a = 0.02: e^s = a(1-p)/(p(1-a))", es, "2.0204")
    check("n = 1000, a = 0.02: s", log(es), "0.7033")
    check("n = 1000, a = 0.02: (1 - p)/(1 - a)", 0.99 / 0.98, "1.0102")
    check("n = 1000, a = 0.02: D(a||p) in nats", kl(a, 0.01), "0.0039136")
    check("n = 1000, a = 0.02: n D", 1000 * kl(a, 0.01), "3.9136")
    check("n = 1000: Chernoff bound on P(X >= 20)", chernoff(1000, a, 0.01), "0.0200")
    check("n = 10000: Markov bound on P(X >= 200)", markov(100, 200), "0.5")
    check("n = 10000: Chebyshev bound on P(X >= 200)", chebyshev(99, 100), "0.0099")
    check("n = 10000: Chernoff bound on P(X >= 200)", chernoff(10000, a, 0.01), "1.008e-17")
    check("D(1/2 || 0.1)", kl(0.5, 0.1), "0.5108")
    check("-ln(0.36)/2", -log(0.36) / 2, "0.5108")
    check("R25 at p = 0.1: bound 0.36^12.5", 0.36**12.5, "2.843e-6")
    check("R25 at p = 0.1: exact", rep_error(25, p), "1.621e-7")
    need = 2 * log(1e9) / log(1 / 0.36)
    check("R_n target 1e-9: 2 ln(1e9)/ln(1/0.36)", need, "40.57")
    check("R41 bound 0.36^20.5", 0.36**20.5, "8.02e-10")
    for n, stated in ((33, "2.396e-9"), (35, "8.392e-10"), (41, "3.636e-11")):
        check(f"R{n} at p = 0.1: exact", rep_error(n, p), stated)
    check("R39 bound 0.36^19.5", 0.36**19.5, "2.23e-9")

    print("Worked examples")
    r = Fraction(1, 10**4)
    check("sector, 4096 bits at 1e-4: mean", 4096 * r, "0.4096")
    check("sector: variance", 4096 * r * (1 - r), "0.40956")
    check("sector: P(X = 0)", binom_pmf(4096, 0, r), "0.66390")
    check("sector: P(X = 1)", binom_pmf(4096, 1, r), "0.27196")
    check("sector: P(X >= 2)", binom_tail(4096, 2, r), "0.064136")
    check("sector: Markov bound on P(X >= 2)", markov(0.4096, 2), "0.2048")
    check("sector: C(4096,2)", comb(4096, 2), "8,386,560")
    check("sector: union bound C(4096,2) p^2", comb(4096, 2) * r**2, "0.08387")
    check("phi(3.5)", phi(3.5), "8.727e-4")
    check("x = 3.5: upper bound phi/x", q_upper(3.5), "2.493e-4")
    check("x = 3.5: lower bound x phi/(1 + x^2)", q_lower(3.5), "2.305e-4")
    check("x = 3.5: 3.5/13.25", 3.5 / 13.25, "0.2642")
    check("Q(3.5)", Q(3.5), "2.326e-4")
    s = Fraction(1, 1000)
    check("31-bit block, corrects 2, p = 1e-3: C(31,3)", comb(31, 3), "4495")
    check("31-bit block: union bound C(31,3) p^3", comb(31, 3) * s**3, "4.495e-6")
    check("31-bit block: exact P(X >= 3)", binom_tail(31, 3, s), "4.4016e-6")
    check("31-bit block: ratio bound / exact", Fraction(comb(31, 3)) * s**3 / binom_tail(31, 3, s), "1.021")

    print("Exercises")
    even = {x for x in space if x.count("1") % 2 == 0}
    check("P(even number of flips) at p = 0.1", prob(even, p), "0.756")
    check("P(A and even)", prob(A & even, p), "0.018")
    check("P(A)P(even)", prob(A, p) * prob(even, p), "0.0756")
    half = Fraction(1, 2)
    assert prob(A & even, half) == prob(A, half) * prob(even, half) == Fraction(1, 4)
    q1000 = (1 - Fraction(1, 10**4)) ** 1000
    check("ARQ: q = (1 - 1e-4)^1000", q1000, "0.904833")
    check("ARQ: E = 1/q", 1 / q1000, "1.105176")
    pe, eps6 = 1e-6, 1e-8
    nbits = pe * (1 - pe) / (eps6**2 * 0.05)
    check("BER p = 1e-6 to 1% at 95%: bits", nbits, "199,999,800,000")
    check("  seconds at 8 Gb/s", nbits / 8e9, "25.0")
    check("R7 at p = 0.1: union bound C(7,4) p^4", comb(7, 4) * p**4, "0.0035")
    check("R7 at p = 0.1: exact", rep_error(7, p), "0.002728")
    check("R7 at p = 0.1: Chernoff 0.36^3.5", 0.36**3.5, "0.02799")
    mu_tight, a_tight = 0.4, 2
    check("Markov tight: P(X = a) = mu/a with mu = 0.4, a = 2", mu_tight / a_tight, "0.2")
    # Ratios and small counts the prose states in words.
    check("patterns with 5 or more flips in 8 bits", sum(comb(8, k) for k in range(5, 9)), "93")
    check("one decision error in about this many bits at Q(3.5)", 1 / Q(3.5), "4.30e3")
    exact_1000 = float(binom_tail(1000, 20, q))
    exact_10000 = float(binom_tail(10000, 200, q))
    check("Markov over exact, n = 1000 (more than a hundred)", 0.5 / exact_1000, "152")
    check("Chebyshev over Markov improvement, n = 1000", 0.5 / 0.099, "5")
    check("Chebyshev over exact, n = 1000 (about thirty)", 0.099 / exact_1000, "30")
    check("exact over Gaussian stand-in, n = 1000 (about four)", exact_1000 / Q(10 / sqrt(9.9)), "4")
    check("exact over Gaussian stand-in, n = 10000", exact_10000 / Q(100 / sqrt(99)), "1.2e5")
    check("exact tail ratio between the two columns", exact_1000 / exact_10000, "5.8e15")
    check("Chernoff over exact, n = 10000 (about 18)", chernoff(10000, 0.02, 0.01) / exact_10000, "18")
    check("Chebyshev over Chernoff, n = 1000 (five times)", 0.099 / chernoff(1000, 0.02, 0.01), "5")
    check("n = 10000: exponent nD", 10000 * kl(0.02, 0.01), "39.136")
    check("R25 bound over exact (about eighteen)", 0.36**12.5 / float(rep_error(25, p)), "18")
    print("All chapter 2 numbers reproduced.")


if __name__ == "__main__":
    main()
