"""Reproduce every number in chapter 7 of Error Correction Coding, "Capacity and Shannon's Theorem".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch07_capacity_and_shannons_theorem.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. Binomial sums are
computed exactly with integers and fractions, the small random code is drawn
from a fixed seed and decoded exactly over every received word, and the BPSK
capacity is a Simpson-rule integral with 4000 intervals. A run that reaches the
final line has reproduced every number in the chapter.

The module also holds the capacity, typical-set and random-coding functions
that ch07_charts.py plots, so the charts and this check share one source.
"""

import random
import sys
from fractions import Fraction
from itertools import combinations
from math import ceil, comb, e, erfc, exp, lgamma, log, log1p, log2, log10, pi, sqrt
from statistics import NormalDist


# Helpers ------------------------------------------------------------------------

def half_unit(stated: str) -> float:
    text = stated.replace(",", "")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", ""))
    print(f"  {label:<62} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


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


# Information measures ----------------------------------------------------------------

def binary_entropy(p: float) -> float:
    if p <= 0 or p >= 1:
        return 0.0
    return -p * log2(p) - (1 - p) * log2(1 - p)


def db(x: float) -> float:
    return 10 * log10(x)


def q_function(x: float) -> float:
    return 0.5 * erfc(x / sqrt(2))


def q_inverse(prob: float) -> float:
    return -NormalDist().inv_cdf(prob)


def log2_comb(n: int, k: int) -> float:
    return log2(comb(n, k))


# Capacities of the binary channels ------------------------------------------------

def bsc_output_one(q: float, p: float) -> float:
    """P(Y = 1) for a BSC(p) whose input is 1 with probability q."""
    return q * (1 - p) + (1 - q) * p


def bsc_mutual_information(q: float, p: float) -> float:
    return binary_entropy(bsc_output_one(q, p)) - binary_entropy(p)


def bec_mutual_information(q: float, eps: float) -> float:
    return (1 - eps) * binary_entropy(q)


def z_mutual_information(q: float, p: float) -> float:
    """Z-channel: a 1 turns into 0 with probability p, a 0 is never disturbed."""
    return binary_entropy(q * (1 - p)) - q * binary_entropy(p)


def z_capacity(p: float) -> float:
    return log2(1 + 2 ** (-binary_entropy(p) / (1 - p)))


def z_best_output(p: float) -> float:
    return 1 / (1 + 2 ** (binary_entropy(p) / (1 - p)))


def bsc_capacity(p: float) -> float:
    return 1 - binary_entropy(p)


def bisect(f, lo: float, hi: float, target: float, steps: int = 80) -> float:
    """The x in [lo, hi] where the increasing function f reaches target."""
    for _ in range(steps):
        mid = (lo + hi) / 2
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# The Gaussian channel ------------------------------------------------------------

def gaussian_capacity(snr: float) -> float:
    """Bits per real dimension, snr = P/N."""
    return 0.5 * log2(1 + snr)


def bpsk_soft_capacity(snr: float, intervals: int = 4000) -> float:
    """1 - E[log2(1 + exp(-2Y/s2))] for Y ~ N(1, s2) and s2 = 1/snr, by Simpson's rule on +-12 sigma."""
    s2 = 1 / snr
    s = sqrt(s2)
    a, b = 1 - 12 * s, 1 + 12 * s
    step = (b - a) / intervals
    total = 0.0
    for i in range(intervals + 1):
        y = a + i * step
        weight = 1 if i in (0, intervals) else (4 if i % 2 else 2)
        z = -2 * y / s2
        f = (z + log1p(exp(-z))) / log(2) if z > 0 else log1p(exp(z)) / log(2)
        total += weight * f * exp(-(y - 1) ** 2 / (2 * s2)) / sqrt(2 * pi * s2)
    return 1 - total * step / 3


def bpsk_hard_capacity(snr: float) -> float:
    return bsc_capacity(q_function(sqrt(snr)))


def min_ebn0_db(capacity, rate: float) -> float:
    """Smallest Eb/N0 in dB at which capacity(snr) reaches rate, with Eb/N0 = snr / (2 rate)."""
    snr = bisect(capacity, 1e-6, 1e3, rate, 70)
    return db(snr / (2 * rate))


def shannon_hartley(bandwidth: float, snr: float) -> float:
    return bandwidth * log2(1 + snr)


def ebn0_boundary(eta: float) -> float:
    """Smallest Eb/N0 (as a ratio) at spectral efficiency eta bits per second per hertz."""
    return (2 ** eta - 1) / eta


# Typical sequences ---------------------------------------------------------------

def typical_counts(n: int, p: Fraction, window: Fraction) -> tuple:
    """Ones counts k with |k/n - p| <= window, the number of such sequences and their total probability."""
    ks = [k for k in range(n + 1) if abs(Fraction(k, n) - p) <= window]
    size = sum(comb(n, k) for k in ks)
    prob = sum(Fraction(comb(n, k)) * p ** k * (1 - p) ** (n - k) for k in ks)
    return ks, size, prob


def typical_probability_float(n: int, p: float, window: float) -> float:
    """The same probability through log-gamma, fast enough for the chart."""
    total = 0.0
    for k in range(n + 1):
        if abs(k / n - p) <= window + 1e-12:
            total += exp(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + k * log(p) + (n - k) * log(1 - p))
    return total


# Random coding on the BSC --------------------------------------------------------

def ball_size(n: int, r: int) -> int:
    return sum(comb(n, k) for k in range(r + 1))


def binomial_tail(n: int, p: Fraction, r: int) -> Fraction:
    """P(K > r) for K ~ Binomial(n, p), exactly."""
    return 1 - sum(Fraction(comb(n, k)) * p ** k * (1 - p) ** (n - k) for k in range(r + 1))


def log_add(a: float, b: float) -> float:
    if a == float("-inf"):
        return b
    if b == float("-inf"):
        return a
    m = max(a, b)
    return m + log1p(exp(-abs(a - b)))


def sphere_decoder_bound(n: int, rate: float, p: float) -> float:
    """Average error of a random code of 2^(nR) words with the radius-r decoder, minimized over r.

    The bound is P(K > r) + (2^(nR) - 1) |ball(r)| / 2^n, each term computed in the natural-log domain.
    """
    log_pmf = [lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + k * log(p) + (n - k) * log(1 - p)
               for k in range(n + 1)]
    log_binom = [lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) for k in range(n + 1)]
    tails = [float("-inf")] * (n + 2)
    for k in range(n, -1, -1):
        tails[k] = log_add(tails[k + 1], log_pmf[k])
    best = 1.0
    log_ball = float("-inf")
    log_words = n * rate * log(2) + log1p(-2.0 ** (-n * rate))
    for r in range(0, n // 2 + 1):
        log_ball = log_add(log_ball, log_binom[r])
        if r < int(n * p):
            continue
        union = log_words + log_ball - n * log(2)
        tail = tails[r + 1]
        value = exp(tail) + (exp(union) if union < 0 else 1.0)
        best = min(best, value)
    return best


def random_code(seed: int, words: int, n: int) -> list:
    rng = random.Random(seed)
    return [rng.getrandbits(n) for _ in range(words)]


def weight(x: int) -> int:
    return bin(x).count("1")


def ml_block_error(code: list, n: int, p: float) -> float:
    """Exact block error of maximum-likelihood (minimum-distance) decoding on a BSC(p), ties split evenly."""
    pw = [p ** w * (1 - p) ** (n - w) for w in range(n + 1)]
    correct = 0.0
    for y in range(1 << n):
        dists = [weight(c ^ y) for c in code]
        # Each of the w codewords tied at the minimum distance d wins y with probability 1/w, and
        # each has likelihood pw[d], so y adds w * pw[d] / w = pw[d] to the summed correct probability.
        correct += pw[min(dists)]
    return 1 - correct / len(code)


# Fano and the converse ------------------------------------------------------------

def fano_bound(pe: float, m: int) -> float:
    return binary_entropy(pe) + pe * log2(m - 1)


def converse_bound(capacity: float, rate: float, n: float) -> float:
    return 1 - capacity / rate - 1 / (n * rate)


def dispersion_bsc(p: float) -> float:
    return p * (1 - p) * log2((1 - p) / p) ** 2


def normal_approximation(n: int, eps: float, p: float) -> float:
    return bsc_capacity(p) - sqrt(dispersion_bsc(p) / n) * q_inverse(eps)


# The checks ---------------------------------------------------------------------

def check_counting() -> None:
    section("Counting noise patterns")
    h = binary_entropy(0.1)
    check("h(0.1)", h, "0.4690")
    check("1 - h(0.1)", 1 - h, "0.5310")
    for n, stated, per, nh in ((100, "43.98", "0.4398", "46.90"), (1000, "464.42", "0.4644", "469.00"),
                               (10000, "4683.72", "0.4684", "4689.96")):
        check(f"log2 C({n}, {n // 10})", log2_comb(n, n // 10), stated)
        check(f"log2 C({n}, {n // 10}) / n", log2_comb(n, n // 10) / n, per)
        check(f"{n} h(0.1)", n * h, nh)
    check("1000 - 1000 h(0.1), exponent of the codeword count", 1000 - 1000 * h, "531.00")
    check("chapter 6 source over BSC(0.1), uses per letter", h / (1 - h), "0.8832")


def check_binary_channels() -> None:
    section("The binary symmetric and erasure channels")
    check("h(0.11)", binary_entropy(0.11), "0.4999")
    check("capacity of BSC(0.11)", bsc_capacity(0.11), "0.5001")
    check("crossover with capacity exactly 0.5", bisect(lambda p: binary_entropy(p), 0.01, 0.5, 0.5), "0.110028")
    rows = ((0.1, "0.188", "0.6973", "0.1974", "0.2345"), (0.2, "0.266", "0.8357", "0.3358", "0.3610"),
            (0.3, "0.344", "0.9286", "0.4287", "0.4406"), (0.4, "0.422", "0.9824", "0.4825", "0.4855"),
            (0.5, "0.500", "1.0000", "0.5001", "0.5000"))
    for q, r, hr, i_bsc, i_bec in rows:
        check(f"q = {q}: P(Y = 1) for BSC(0.11)", bsc_output_one(q, 0.11), r)
        check(f"q = {q}: h(P(Y = 1))", binary_entropy(bsc_output_one(q, 0.11)), hr)
        check(f"q = {q}: I(X;Y) for BSC(0.11)", bsc_mutual_information(q, 0.11), i_bsc)
        check(f"q = {q}: I(X;Y) for BEC(0.5)", bec_mutual_information(q, 0.5), i_bec)
    check("BSC(0.1) at q = 0.5", bsc_mutual_information(0.5, 0.1), "0.5310")


def check_gaussian() -> None:
    section("The Gaussian channel")
    check("uniform, unit variance, h in bits", log2(2 * sqrt(3)), "1.7925")
    check("Laplace, unit variance, h in bits", log2(sqrt(2) * e), "1.9427")
    check("Gaussian, unit variance, h in bits", 0.5 * log2(2 * pi * e), "2.0471")
    check("Gaussian minus uniform", 0.5 * log2(2 * pi * e) - log2(2 * sqrt(3)), "0.2546")
    check("half log2(pi e / 6)", 0.5 * log2(pi * e / 6), "0.2546")
    check("C at P/N = 3", gaussian_capacity(3), "1")
    check("C at P/N = 15", gaussian_capacity(15), "2")
    check("C at P/N = 63", gaussian_capacity(63), "3")
    check("P/N = 3 in dB", db(3), "4.77")
    check("P/N = 15 in dB", db(15), "11.76")
    check("P/N = 63 in dB", db(63), "17.99")
    check("factor 4 in dB", db(4), "6.02")
    rows = ((-10, "0.0688", "0.0687", "0.0449"), (-5, "0.1982", "0.1977", "0.1353"), (0, "0.5000", "0.4859", "0.3689"),
            (5, "1.0287", "0.8592", "0.7685"), (10, "1.7297", "0.9968", "0.9908"))
    for snr_db, cg, cs, ch in rows:
        snr = 10 ** (snr_db / 10)
        check(f"{snr_db} dB: Gaussian input", gaussian_capacity(snr), cg)
        check(f"{snr_db} dB: BPSK, soft", bpsk_soft_capacity(snr), cs)
        check(f"{snr_db} dB: BPSK, hard", bpsk_hard_capacity(snr), ch)
    require("Gaussian >= soft >= hard at every listed SNR",
            all(gaussian_capacity(10 ** (d / 10)) >= bpsk_soft_capacity(10 ** (d / 10)) >= bpsk_hard_capacity(10 ** (d / 10))
                for d in range(-10, 11)))


def check_minimum_ebn0() -> None:
    section("Minimum Eb/N0 by rate")
    gauss = {0.5: "0.00", 1 / 3: "-0.55", 0.1: "-1.29", 0.01: "-1.56"}
    soft = {0.5: "0.19", 1 / 3: "-0.50", 0.1: "-1.29", 0.01: "-1.56"}
    hard = {0.5: "1.77", 1 / 3: "1.21", 0.1: "0.59", 0.01: "0.39"}
    for rate in (0.5, 1 / 3, 0.1, 0.01):
        check(f"R = {rate:.3g}: Gaussian input, (2^(2R) - 1)/(2R) in dB", db(ebn0_boundary(2 * rate)), gauss[rate])
        check(f"R = {rate:.3g}: BPSK soft", min_ebn0_db(bpsk_soft_capacity, rate), soft[rate])
        check(f"R = {rate:.3g}: BPSK hard", min_ebn0_db(bpsk_hard_capacity, rate), hard[rate])
    soft_half = min_ebn0_db(bpsk_soft_capacity, 0.5)
    hard_half = min_ebn0_db(bpsk_hard_capacity, 0.5)
    check("rate 1/2, soft, three places", soft_half, "0.187")
    check("rate 1/2, hard, three places", hard_half, "1.773")
    check("rate 1/2, hard minus soft", hard_half - soft_half, "1.585")
    snr_hard = bisect(bpsk_hard_capacity, 1e-6, 1e3, 0.5, 70)
    check("rate 1/2, crossover of the hard BSC", q_function(sqrt(snr_hard)), "0.1100")
    check("low-rate hard-decision penalty, 10 log10(pi/2)", db(pi / 2), "1.96")
    check("rate 0.01, hard minus soft", min_ebn0_db(bpsk_hard_capacity, 0.01) - min_ebn0_db(bpsk_soft_capacity, 0.01),
          "1.95")


def check_bandwidth() -> None:
    section("Bandwidth and power")
    snr = 10 ** 2.5
    check("25 dB as a ratio", snr, "316.23")
    check("log2(1 + 316.23)", log2(1 + snr), "8.3094")
    check("20 MHz at 25 dB, Mbit/s", shannon_hartley(20, snr), "166.19")
    check("half the S/N", snr / 2, "158.11")
    check("log2(1 + 158.11)", log2(1 + snr / 2), "7.3139")
    check("40 MHz at the same power, Mbit/s", shannon_hartley(40, snr / 2), "292.56")
    check("double the S/N", snr * 2, "632.46")
    check("log2(1 + 632.46)", log2(1 + 2 * snr), "9.3071")
    check("20 MHz at double power, Mbit/s", shannon_hartley(20, 2 * snr), "186.14")
    for b, stated in ((1e5, "0.3459"), (1e6, "1.0000"), (1e7, "1.3750"), (1e8, "1.4355")):
        check(f"S/N0 = 1e6, B = {b:.0e}, Mbit/s", shannon_hartley(b, 1e6 / b) / 1e6, stated)
    check("log2 e, the wideband limit in Mbit/s", log2(e), "1.4427")
    rows = ((8, "31.875", "15.03"), (6, "10.500", "10.21"), (4, "3.750", "5.74"), (2, "1.500", "1.76"),
            (1, "1.000", "0.00"), (0.5, "0.828", "-0.82"), (0.25, "0.757", "-1.21"), (0.1, "0.718", "-1.44"),
            (0.01, "0.696", "-1.58"), (0.001, "0.693", "-1.59"))
    for eta, ratio, in_db in rows:
        check(f"eta = {eta}: (2^eta - 1)/eta", ebn0_boundary(eta), ratio)
        check(f"eta = {eta}: in dB", db(ebn0_boundary(eta)), in_db)
    check("eta = 3, in dB", db(ebn0_boundary(3)), "3.68")
    check("uncoded BPSK at 1e-5, dB (chapter 3)", db(q_inverse(1e-5) ** 2 / 2), "9.59")
    check("uncoded BPSK minus the bound at eta = 2", db(q_inverse(1e-5) ** 2 / 2) - db(ebn0_boundary(2)), "7.83")
    check("uncoded BPSK minus the Shannon limit", db(q_inverse(1e-5) ** 2 / 2) - db(log(2)), "11.18")


def check_shannon_limit() -> None:
    section("The Shannon limit")
    check("ln 2", log(2), "0.6931")
    check("ln 2, six places", log(2), "0.693147")
    check("log10(ln 2)", log10(log(2)), "-0.15917")
    check("10 log10(ln 2)", db(log(2)), "-1.5917")
    check("10 log10(ln 2), two places", db(log(2)), "-1.59")
    check("(ln 2)^2 / 2, slope of the series", log(2) ** 2 / 2, "0.2402")
    check("eta = 0.1 from the series, ln 2 + 0.1 (ln 2)^2/2", log(2) + 0.1 * log(2) ** 2 / 2, "0.7172")
    check("eta = 0.1 exactly", ebn0_boundary(0.1), "0.7177")
    require("the boundary decreases as eta falls, on a grid from 10 to 0.001",
            all(ebn0_boundary(10 ** (-3 + i / 100)) < ebn0_boundary(10 ** (-3 + (i + 1) / 100)) for i in range(400)))


def check_typical() -> None:
    section("Typical sequences")
    p = Fraction(1, 10)
    h = binary_entropy(0.1)
    check("log2 9", log2(9), "3.1699")
    check("delta = 0.05 log2 9", 0.05 * log2(9), "0.1585")
    rows = ((0, "1", "0.12158", "0.12158", "0.1520"), (1, "20", "0.013509", "0.27017", "0.3105"),
            (2, "190", "0.0015009", "0.28518", "0.4690"), (3, "1140", "0.00016677", "0.19012", "0.6275"),
            (4, "4845", "0.000018530", "0.08978", "0.7860"), (5, "15504", "0.0000020589", "0.03192", "0.9445"))
    for k, count, each, total, rate in rows:
        each_exact = p ** k * (1 - p) ** (20 - k)
        same(f"n = 20, k = {k}: C(20, k)", comb(20, k), int(count))
        check(f"n = 20, k = {k}: probability of each", each_exact, each)
        check(f"n = 20, k = {k}: total probability", comb(20, k) * each_exact, total)
        check(f"n = 20, k = {k}: -(1/n) log2 p(x)", -log2(each_exact) / 20, rate)
        check(f"n = 20, k = {k}: h + (k/n - p) log2 9", h + (k / 20 - 0.1) * log2(9), rate)
    check("2^(20 h)", 2 ** (20 * h), "666.25")
    check("2^(20 (h - delta))", 2 ** (20 * (h - 0.05 * log2(9))), "74.03")
    check("2^(20 (h + delta))", 2 ** (20 * (h + 0.05 * log2(9))), "5996.2")
    ks, size, prob = typical_counts(20, p, Fraction(1, 20))
    same("n = 20: typical ones counts", ks, [1, 2, 3])
    same("n = 20: |T|", size, 1350)
    check("n = 20: log2 |T|", log2(size), "10.399")
    check("n = 20: P(T)", prob, "0.74547")
    check("n = 20: n (h - delta)", 20 * (h - 0.05 * log2(9)), "6.21")
    check("n = 20: n (h + delta)", 20 * (h + 0.05 * log2(9)), "12.55")
    for n, lo, hi, pt, lt, b_lo, b_hi in ((100, 5, 15, "0.93640", "58.086", "31.05", "62.75"),
                                          (1000, 50, 150, "0.99999972", "605.30", "310.50", "627.49")):
        ks, size, prob = typical_counts(n, p, Fraction(1, 20))
        same(f"n = {n}: range of typical k", (ks[0], ks[-1]), (lo, hi))
        check(f"n = {n}: P(T)", prob, pt)
        check(f"n = {n}: log2 |T|", log2(size), lt)
        check(f"n = {n}: n (h - delta)", n * (h - 0.05 * log2(9)), b_lo)
        check(f"n = {n}: n (h + delta)", n * (h + 0.05 * log2(9)), b_hi)
    ks, size, prob = typical_counts(10000, p, Fraction(1, 20))
    check("n = 10000: log2 |T|", log2(size), "6092.2")
    require("n = 10000: P(T) above 1 - 1e-40", prob > 1 - Fraction(1, 10 ** 40))
    ks, size, prob = typical_counts(10, p, Fraction(1, 20))
    same("n = 10: typical ones counts", ks, [1])
    check("n = 10: probability of each typical sequence", p * (1 - p) ** 9, "0.038742")
    check("n = 10: P(T)", prob, "0.38742")
    check("n = 10: 2^(10 h)", 2 ** (10 * h), "25.81")


def check_random_coding() -> None:
    section("The random coding argument")
    n, r, p = 1000, 140, Fraction(1, 10)
    check("2^400", 2.0 ** 400, "2.58e120")
    tail = binomial_tail(n, p, r)
    check("P(K > 140), n = 1000, p = 0.1", tail, "2.424e-5")
    lg = log2(ball_size(n, r)) - n
    check("log2 P(a wrong codeword within 140)", lg, "-420.29")
    check("1 - h(0.14)", 1 - binary_entropy(0.14), "0.4158")
    check("bound -1000 (1 - h(0.14))", -n * (1 - binary_entropy(0.14)), "-415.76")
    check("union exponent, exact", 400 + lg, "-20.29")
    check("union term, exact", 2.0 ** (400 + lg), "7.8e-7")
    check("union exponent from the h bound", 400 - n * (1 - binary_entropy(0.14)), "-15.76")
    check("union term from the h bound", 2.0 ** (400 - n * (1 - binary_entropy(0.14))), "1.8e-5")
    total = float(tail) + 2.0 ** (400 + lg)
    check("average error bound, both terms", total, "2.50e-5")
    check("P(K > 125)", binomial_tail(n, p, 125), "0.00453")
    check("log2 P(wrong within 125) + 400", 400 + log2(ball_size(n, 125)) - n, "-60.93")
    require("ball bound sum C(n,k), k <= rho n, <= 2^(n h(rho)) for n = 1000, rho = 0.14",
            log2(ball_size(1000, 140)) <= 1000 * binary_entropy(0.14))
    require("same bound for n = 1000 and every rho = r/n up to 1/2",
            all(log2(ball_size(1000, rr)) <= 1000 * binary_entropy(rr / 1000) + 1e-9 for rr in range(1, 501)))
    check("sphere-decoder bound, n = 1000, R = 0.4", sphere_decoder_bound(1000, 0.4, 0.1), "2.1e-5")


def check_small_code() -> None:
    section("A small random code decoded exactly")
    code = random_code(1, 16, 12)
    same("seed 1, first codeword", format(code[0], "012b"), "001000100110")
    same("seed 1, last codeword", format(code[-1], "012b"), "001101011011")
    same("distinct codewords", len(set(code)), 16)
    dists = [weight(a ^ b) for a, b in combinations(code, 2)]
    same("minimum distance", min(dists), 2)
    same("pairs at distance 2", dists.count(2), 2)
    pairs = sorted(tuple(sorted((format(a, "012b"), format(b, "012b")))) for a, b in combinations(code, 2)
                   if weight(a ^ b) == 2)
    same("the two close pairs", pairs, [("110000101100", "110000111000"), ("110010011110", "110011010110")])
    check("block error of ML decoding, n = 12, M = 16", ml_block_error(code, 12, 0.1), "0.10539")
    check("uncoded 4 bits, block error", 1 - 0.9 ** 4, "0.3439")
    check("R3 per bit, 3p^2(1-p) + p^3", 3 * 0.01 * 0.9 + 0.001, "0.028")
    check("R3 on 4 bits, block error", 1 - 0.972 ** 4, "0.10738")
    code16 = random_code(1, 16, 16)
    same("n = 16: minimum distance", min(weight(a ^ b) for a, b in combinations(code16, 2)), 3)
    check("n = 16, R = 1/4: block error", ml_block_error(code16, 16, 0.1), "0.05220")


def check_fano() -> None:
    section("Fano's inequality and the converse")
    check("log2 3", log2(3), "1.5850")
    check("Fano, M = 4, Pe = 0.1", fano_bound(0.1, 4), "0.6275")
    check("h(0.05)", binary_entropy(0.05), "0.2864")
    check("0.05 log2 7", 0.05 * log2(7), "0.1404")
    check("Fano, M = 8, Pe = 0.05", fano_bound(0.05, 8), "0.4268")
    c = bsc_capacity(0.1)
    check("1 - C/R, BSC(0.1), R = 0.6", 1 - c / 0.6, "0.1150")
    check("converse, n = 100", converse_bound(c, 0.6, 100), "0.0983")
    check("converse, n = 1000", converse_bound(c, 0.6, 1000), "0.1133")
    check("BSC(0.2) capacity", bsc_capacity(0.2), "0.2781")
    check("converse limit, BSC(0.2), R = 1/3", 1 - bsc_capacity(0.2) / (1 / 3), "0.1658")
    check("R3 bit error on BSC(0.2)", 3 * 0.04 * 0.8 + 0.008, "0.104")
    check("converse limit, BSC(0.11), R = 0.75", 1 - bsc_capacity(0.11) / 0.75, "0.3332")


def check_sources() -> None:
    section("Sources over channels")
    check("h(0.2)", binary_entropy(0.2), "0.7219")
    check("source p = 0.2 over BSC(0.1), uses per letter", binary_entropy(0.2) / bsc_capacity(0.1), "1.3596")
    check("D / sigma^2 at P/N = 15, one use", 1 / 16, "0.0625")
    check("in dB", db(1 / 16), "-12.04")
    check("D / sigma^2 at P/N = 15, two uses", 1 / 256, "0.00390625")
    check("in dB", db(1 / 256), "-24.08")
    check("uncoded scaling with MMSE: N/(P+N) at P/N = 15", 1 / (1 + 15), "0.0625")


def check_open_questions() -> None:
    section("What the theorem leaves open")
    v = dispersion_bsc(0.1)
    check("dispersion V of BSC(0.1)", v, "0.90436")
    check("sqrt V", sqrt(v), "0.95098")
    check("Q^-1(1e-3)", q_inverse(1e-3), "3.0902")
    check("Q^-1(1e-6)", q_inverse(1e-6), "4.7534")
    for n, a, b in ((100, "0.2371", "0.0790"), (1000, "0.4381", "0.3881"), (10000, "0.5016", "0.4858")):
        check(f"normal approximation, n = {n}, eps = 1e-3", normal_approximation(n, 1e-3, 0.1), a)
        check(f"normal approximation, n = {n}, eps = 1e-6", normal_approximation(n, 1e-6, 0.1), b)
    check("dropped term log2(1000) / 2000", log2(1000) / 2000, "0.0050")
    check("fraction of capacity at n = 1000, eps = 1e-3", normal_approximation(1000, 1e-3, 0.1) / bsc_capacity(0.1),
          "0.825")


def check_examples() -> None:
    section("Worked examples")
    p = 0.1
    check("Z-channel: h(p)/(1 - p)", binary_entropy(p) / (1 - p), "0.52111")
    check("Z-channel: 2^(-h/(1-p))", 2 ** (-binary_entropy(p) / (1 - p)), "0.69684")
    check("Z-channel: best P(Y = 1)", z_best_output(p), "0.41067")
    check("Z-channel: best q", z_best_output(p) / (1 - p), "0.45630")
    check("Z-channel: capacity, closed form", z_capacity(p), "0.76285")
    check("Z-channel: I at the best q", z_mutual_information(z_best_output(p) / (1 - p), p), "0.76285")
    check("Z-channel: I at q = 0.5", z_mutual_information(0.5, p), "0.7583")
    best_grid = max(range(1, 1000), key=lambda i: z_mutual_information(i / 1000, p))
    same("Z-channel: best q on a grid of 0.001", best_grid, 456)
    check("BEC(0.25) capacity", 0.75, "0.75")
    check("uses for 1000 bits on BEC(0.25)", 1000 / 0.75, "1333.3")


def check_prose_details() -> None:
    """Numbers the prose states in passing, outside the tables and the main derivations."""
    section("Numbers stated in passing")
    h = binary_entropy(0.1)
    # Counting noise patterns: the shortfall n h - log2 C(n, n/10) grows far more slowly than n.
    for n, stated in ((100, "2.92"), (1000, "4.57"), (10000, "6.23")):
        check(f"shortfall n h - log2 C(n, n/10), n = {n}", n * h - log2_comb(n, n // 10), stated)
    check("BSC crossover for capacity one half, mirrored", 1 - bisect(binary_entropy, 0.0, 0.5, 0.5), "0.890")
    check("BSC(0.11): P(Y = 1) term 0.89", 1 - 0.11, "0.89")
    check("flip cost h(0.11)/0.11", binary_entropy(0.11) / 0.11, "4.54")
    # Gaussian and BPSK.
    check("Gaussian capacity at 20 dB", gaussian_capacity(10 ** 2), "3.33")
    check("pi/2", pi / 2, "1.571")
    check("soft over hard BPSK capacity at -10 dB", bpsk_soft_capacity(0.1) / bpsk_hard_capacity(0.1), "1.53")
    check("table ratio 0.0687/0.0449", 0.0687 / 0.0449, "1.53")
    # Bandwidth against power at 25 dB.
    snr = 10 ** 2.5
    base = shannon_hartley(20, snr)
    check("doubled bandwidth: S/N + 1", snr / 2 + 1, "159.11")
    check("doubled power: S/N + 1", 2 * snr + 1, "633.46")
    check("bandwidth gain, percent", 100 * (shannon_hartley(40, snr / 2) / base - 1), "76")
    check("power gain, percent", 100 * (shannon_hartley(20, 2 * snr) / base - 1), "12")
    # Typical sequences.
    check("delta = 0.05 log2 9", 0.05 * log2(9), "0.1585")
    check("variance p(1 - p)", 0.1 * 0.9, "0.09")
    check("window squared", 0.05 ** 2, "0.0025")
    check("Chebyshev constant 0.09/0.0025", 0.09 / 0.0025, "36")
    check("Chebyshev bound at n = 1000", 36 / 1000, "0.036")
    _, _, prob1000 = typical_counts(1000, Fraction(1, 10), Fraction(1, 20))
    check("exact atypical probability at n = 1000", float(1 - prob1000), "2.8e-7")
    check("n = 10: 2^(10 h) over the 10 typical blocks", 2 ** (10 * h) / 10, "2.6")
    check("typical probability at n = 10, window 0.05", typical_probability_float(10, 0.1, 0.05), "0.39")
    check("typical probability at n = 100, window 0.05", typical_probability_float(100, 0.1, 0.05), "0.94")
    require("typical probability above 0.999 by n = 400", typical_probability_float(400, 0.1, 0.05) > 0.999)
    check("typical probability at n = 1000, window 0.02", typical_probability_float(1000, 0.1, 0.02), "0.97")
    # Random coding.
    check("codebook size in bits, 2^400 x 1000", 2.0 ** 400 * 1000, "2.58e123")
    pe = sphere_decoder_bound(1000, 0.4, 0.1)
    check("expurgated rate 0.4 - 1/n", 0.4 - 1 / 1000, "0.399")
    check("expurgated maximum error 2 x 2.50e-5", 2 * 2.50e-5, "5.0e-5")
    require("rate 0.30 bound below 1e-16 near n = 1400", sphere_decoder_bound(1400, 0.30, 0.1) < 1e-16
            and sphere_decoder_bound(1200, 0.30, 0.1) > 1e-16)
    check("rate 0.45 bound at n = 4000", sphere_decoder_bound(4000, 0.45, 0.1), "1.1e-7")
    require("rate 0.40 bound at n = 1000 is the figure's 2.1e-5", abs(pe - 2.1e-5) < 0.05e-5)
    check("small code: 2^12 received words", 2 ** 12, "4096")
    # Fano and the converse.
    c2 = bsc_capacity(0.2)
    check("converse on BSC(0.2), R = 1/3, n = 3", converse_bound(c2, 1 / 3, 3), "-0.8342")
    check("1/(nR) at n = 3, R = 1/3", 1 / (3 * (1 / 3)), "1")
    check("R3 bit success on BSC(0.2)", 1 - (3 * 0.2 ** 2 * 0.8 + 0.2 ** 3), "0.896")
    c1 = bsc_capacity(0.1)
    require("converse bound rises with n toward its limit",
            converse_bound(c1, 0.6, 1000) < converse_bound(c1, 0.6, 2000) < 1 - c1 / 0.6)
    # Sources over channels.
    check("SNR with two repetitions of SNR 15", 2 * 15, "30")
    check("uncoded error with two uses", 1 / 31, "0.0323")
    # What the theorem leaves open.
    require("n = 100 at error 1e-3 reaches less than half of capacity",
            normal_approximation(100, 1e-3, 0.1) < 0.5 * c1)
    # Worked examples.
    check("Z-channel: 1 + 2^(-s)", 1 + 2 ** (-binary_entropy(0.1) / 0.9), "1.69684")
    check("Z-channel: loss of the uniform input", z_capacity(0.1) - z_mutual_information(0.5, 0.1), "0.0046")
    same("uses for 1000 bits on BEC(0.25), rounded up", ceil(1000 / 0.75), 1334)


def main() -> None:
    check_counting()
    check_binary_channels()
    check_gaussian()
    check_minimum_ebn0()
    check_bandwidth()
    check_shannon_limit()
    check_typical()
    check_random_coding()
    check_small_code()
    check_fano()
    check_sources()
    check_open_questions()
    check_examples()
    check_prose_details()
    print("All chapter 7 numbers reproduced.")


if __name__ == "__main__":
    main()
