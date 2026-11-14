"""Reproduce every number in chapter 3 of Error Correction Coding, "Channels, Signals and Likelihoods".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch03_channels_signals_and_likelihoods.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. A run that reaches
the final line has reproduced every number in the chapter.

The module also holds the seeded channel simulators that ch03_charts.py draws
in the burst figure, so the chart and this check share one source.
"""

import random
import sys
from itertools import product
from math import comb, erfc, exp, log, log10, pi, sqrt


# Helpers ------------------------------------------------------------------------

def half_unit(stated: str) -> float:
    text = stated.replace(",", "")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", ""))
    print(f"  {label:<60} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<60} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH in '{label}'")
        sys.exit(1)


def Q(x: float) -> float:
    """Upper tail of the standard Gaussian, P(Z > x)."""
    return 0.5 * erfc(x / sqrt(2))


def Q_inv(target: float) -> float:
    """The x with Q(x) = target, by bisection, since Q falls steadily."""
    lo, hi = -10.0, 40.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if Q(mid) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def db(ratio: float) -> float:
    return 10 * log10(ratio)


def undb(value_db: float) -> float:
    return 10 ** (value_db / 10)


def bpsk_ber(ebn0_db: float) -> float:
    """Uncoded BPSK on AWGN, P_b = Q(sqrt(2 Eb/N0)), which soft-decision R3 equals."""
    return Q(sqrt(2 * undb(ebn0_db)))


def r3_copy_p(ebn0_db: float) -> float:
    """Per-copy flip probability of R3 when each copy carries Es = Eb/3."""
    return Q(sqrt(2 * undb(ebn0_db) / 3))


def majority3(p: float) -> float:
    """Probability that two or more of three independent copies are wrong."""
    return 3 * p**2 - 2 * p**3


def hard_r3(ebn0_db: float) -> float:
    return majority3(r3_copy_p(ebn0_db))


def solve_db(f, target: float, lo: float = -5.0, hi: float = 20.0) -> float:
    """The Eb/N0 in dB at which the falling error curve f reaches target."""
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def gauss_pdf(y: float, mean: float, sigma: float) -> float:
    return exp(-((y - mean) ** 2) / (2 * sigma**2)) / (sqrt(2 * pi) * sigma)


def prob_of_llr(L: float) -> float:
    """P(b = 0 | y) from the LLR, 1/(1 + e^-L)."""
    return 1 / (1 + exp(-L))


# Gilbert-Elliott parameters and simulators ----------------------------------------

GE_B, GE_G, GE_PG, GE_PB = 0.01, 0.09, 0.001, 0.1
BURST_BITS = 2000
GE_SEED = 18
BSC_SEED = 8


def ge_stationary(b: float = GE_B, g: float = GE_G) -> tuple:
    return g / (b + g), b / (b + g)


def ge_simulate(n: int = BURST_BITS, seed: int = GE_SEED, b: float = GE_B, g: float = GE_G,
                pg: float = GE_PG, pb: float = GE_PB) -> tuple:
    """States (True for bad) and error flags of n bits, the first state drawn from the stationary law."""
    r = random.Random(seed)
    bad = r.random() < b / (b + g)
    states, errors = [], []
    for _ in range(n):
        states.append(bad)
        errors.append(r.random() < (pb if bad else pg))
        bad = (r.random() >= g) if bad else (r.random() < b)
    return states, errors


def bsc_simulate(n: int = BURST_BITS, seed: int = BSC_SEED, p: float | None = None) -> list:
    if p is None:
        pi_g, pi_b = ge_stationary()
        p = pi_g * GE_PG + pi_b * GE_PB
    r = random.Random(seed)
    return [r.random() < p for _ in range(n)]


def runs(flags: list) -> list:
    """(start, length) of each maximal run of True values."""
    out, start = [], None
    for i, f in enumerate(flags + [False]):
        if f and start is None:
            start = i
        elif not f and start is not None:
            out.append((start, i - start))
            start = None
    return out


def ge_r3_back_to_back(b: float = GE_B, g: float = GE_G, pg: float = GE_PG, pb: float = GE_PB) -> float:
    """P(two or more of three consecutive bits wrong), summed over the eight state paths."""
    pi_g, pi_b = ge_stationary(b, g)
    start = {"G": pi_g, "B": pi_b}
    step = {("G", "G"): 1 - b, ("G", "B"): b, ("B", "G"): g, ("B", "B"): 1 - g}
    flip = {"G": pg, "B": pb}
    total = 0.0
    for path in product("GB", repeat=3):
        weight = start[path[0]] * step[path[0], path[1]] * step[path[1], path[2]]
        p1, p2, p3 = (flip[s] for s in path)
        two_or_more = p1 * p2 * (1 - p3) + p1 * (1 - p2) * p3 + (1 - p1) * p2 * p3 + p1 * p2 * p3
        total += weight * two_or_more
    return total


# The checks ---------------------------------------------------------------------

def main() -> None:
    print("Three discrete channels")
    p = 0.1
    check("BSC(0.1): 1011 received as 1001", 0.9 * 0.9 * 0.1 * 0.9, "0.0729")
    check("BEC(0.1): 1011 received as 1?11", 0.9**3 * 0.1, "0.0729")
    check("Z(0.1): 1011 received as 1001", 0.9 * 1 * 0.1 * 0.9, "0.081")
    check("R3 majority on BSC(0.1)", majority3(p), "0.028")
    check("R3 on BEC(0.1), all three erased", p**3, "0.001")
    check("R3 on BEC(0.1), bit error after a guess", p**3 / 2, "0.0005")
    check("R3 any-1 rule on Z(0.1), sent 1 fails", p**3, "0.001")
    check("R3 any-1 rule on Z(0.1), averaged", p**3 / 2, "0.0005")
    z_major = 3 * p**2 * (1 - p) + p**3
    check("R3 majority on Z(0.1), sent 1 fails", z_major, "0.028")
    check("R3 majority on Z(0.1), averaged", z_major / 2, "0.014")
    check("majority over any-1 on Z(0.1)", (z_major / 2) / (p**3 / 2), "28")
    check("Z(0.1): P(sent 1 and received 0)", 0.5 * p, "0.05")
    check("Z(0.1): P(received 0)", 0.5 + 0.5 * p, "0.55")
    check("Z(0.1): P(X = 1 | Y = 0)", 0.5 * p / (0.5 + 0.5 * p), "0.0909")

    print("\nDecisions by Bayes' rule")
    q = 0.05
    check("alarm: P(Y=1 | X=1) P(X=1)", 0.9 * q, "0.045")
    check("alarm: P(Y=1 | X=0) P(X=0)", 0.1 * (1 - q), "0.095")
    check("alarm: P(Y=1)", 0.9 * q + 0.1 * (1 - q), "0.14")
    check("alarm: P(X=1 | Y=1)", 0.9 * q / (0.9 * q + 0.1 * (1 - q)), "0.3214")
    check("alarm: received 0, alarm product", 0.1 * q, "0.005")
    check("alarm: received 0, no-alarm product", 0.9 * (1 - q), "0.855")
    check("alarm: always-0 receiver error", q, "0.05")
    check("alarm: believe-the-channel receiver error", 0.1, "0.1")
    check("BPSK a=1, sigma=0.5, ML error Q(2)", Q(2), "0.02275")

    print("\nSignal space and the matched filter")
    s = (0.1, 0.7, 0.7, 0.1)
    energy = sum(x * x for x in s)
    check("pulse energy", energy, "1")
    check("second sample, signal", s[1], "0.7")
    check("second sample, error Q(0.7)", Q(0.7), "0.2420")
    check("plain sum, signal", sum(s), "1.6")
    check("plain sum, noise variance", 4 * 1.0, "4")
    check("plain sum, noise standard deviation", sqrt(4), "2")
    check("plain sum, signal over noise", sum(s) / 2, "0.8")
    check("plain sum, error Q(0.8)", Q(sum(s) / 2), "0.2119")
    check("matched, signal", energy, "1")
    check("matched, noise variance", energy * 1.0, "1")
    check("matched, error Q(1)", Q(1), "0.1587")
    sv, yv = (0.6, 0.8), (0.9, -0.1)
    inner = sum(a * c for a, c in zip(sv, yv))
    d_plus = sum((a - c) ** 2 for a, c in zip(yv, sv))
    d_minus = sum((a + c) ** 2 for a, c in zip(yv, sv))
    check("signal space: 0.6 x 0.9", sv[0] * yv[0], "0.54")
    check("signal space: 0.8 x 0.1", -sv[1] * yv[1], "0.08")
    check("signal space: (0.9 - 0.6)^2", (yv[0] - sv[0]) ** 2, "0.09")
    check("signal space: (-0.1 - 0.8)^2", (yv[1] - sv[1]) ** 2, "0.81")
    check("signal space: (0.9 + 0.6)^2", (yv[0] + sv[0]) ** 2, "2.25")
    check("signal space: (-0.1 + 0.8)^2", (yv[1] + sv[1]) ** 2, "0.49")
    check("signal space: inner product", inner, "0.46")
    check("signal space: squared distance to s", d_plus, "0.90")
    check("signal space: squared distance to -s", d_minus, "2.74")
    check("signal space: difference of squared distances", d_minus - d_plus, "1.84")
    require("difference equals four times the inner product", abs((d_minus - d_plus) - 4 * inner) < 1e-12)
    n0 = 2 * 1.0
    check("four-sample example: N0", n0, "2")
    check("four-sample example: Eb/N0", energy / n0, "0.5")
    check("four-sample example: Q(sqrt(2 Eb/N0))", Q(sqrt(2 * energy / n0)), "0.1587")

    print("\nEnergy per bit and the decibel")
    check("ratio 2 in dB", db(2), "3.010")
    check("ratio 1/2 in dB", db(0.5), "-3.010")
    check("ratio 10 in dB", db(10), "10")
    check("ratio 100 in dB", db(100), "20")
    check("ratio 1 in dB", db(1), "0")
    check("plain sum, squared ratio", 0.8**2, "0.64")
    check("single sample, squared ratio", 0.7**2, "0.49")
    check("plain sum loss, ratio 0.64 in dB", -db(0.8**2), "1.94")
    check("single sample loss, ratio 0.49 in dB", -db(0.7**2), "3.10")
    check("chapter 2 example: N0", 2 * 0.25, "0.5")
    check("chapter 2 example: Eb/N0", 1 / 0.5, "2")
    check("chapter 2 example: Eb/N0 in dB", db(2), "3.01")
    check("chapter 2 example: Q(sqrt(4))", Q(sqrt(4)), "0.02275")
    table = [
        (0, "1", "1.4142", "0.07865"),
        (2, "1.5849", "1.7804", "0.03751"),
        (4, "2.5119", "2.2414", "0.01250"),
        (6, "3.9811", "2.8217", "0.002388"),
        (8, "6.3096", "3.5523", "1.909e-4"),
        (10, "10", "4.4721", "3.872e-6"),
        (12, "15.849", "5.6301", "9.006e-9"),
    ]
    for d, ratio, root, pb in table:
        check(f"{d} dB: ratio", undb(d), ratio)
        check(f"{d} dB: sqrt(2 Eb/N0)", sqrt(2 * undb(d)), root)
        check(f"{d} dB: P_b", bpsk_ber(d), pb)
    check("0 dB over 2 dB, about half", bpsk_ber(0) / bpsk_ber(2), "2.1")
    check("10 dB over 12 dB, about 430", bpsk_ber(10) / bpsk_ber(12), "430")
    x5 = Q_inv(1e-5)
    check("Q(x) = 1e-5: x", x5, "4.2649")
    check("Q(x) = 1e-5: x squared", x5**2, "18.189")
    check("Q(x) = 1e-5: Eb/N0", x5**2 / 2, "9.0946")
    check("Q(x) = 1e-5: Eb/N0 in dB", db(x5**2 / 2), "9.588")
    check("about 9.6 dB", db(x5**2 / 2), "9.6")
    check("rate 1/3 in dB", db(1 / 3), "-4.771")
    check("SNR: Rb/B = 1/2 in dB", db(0.5), "-3.010")
    check("SNR for 1e-5 at Rb/B = 1/2", db(x5**2 / 2) + db(0.5), "6.578")

    print("\nHard and soft decisions")
    check("4 dB: per-copy Es/N0 in dB", 4 + db(1 / 3), "-0.771")
    check("4 dB: per-copy Es/N0", undb(4) / 3, "0.8373")
    rows = [
        (4, "0.09782", "0.02684", "0.01250", "2.15"),
        (6, "0.05164", "0.007726", "0.002388", "3.23"),
        (8, "0.02014", "0.0012001", "1.909e-4", "6.29"),
        (10, "0.004912", "7.214e-5", "3.872e-6", "18.6"),
    ]
    for d, pc, hard, soft, ratio in rows:
        check(f"{d} dB: p per copy", r3_copy_p(d), pc)
        check(f"{d} dB: hard R3", hard_r3(d), hard)
        check(f"{d} dB: soft R3 = uncoded", bpsk_ber(d), soft)
        check(f"{d} dB: hard over soft", hard_r3(d) / bpsk_ber(d), ratio)
    # Soft R3 sums three copies: signal 3 sqrt(Es), noise sd sqrt(3) sigma, Es = Eb/3.
    ebn0 = undb(7.3)
    es, sigma = ebn0 / 3, sqrt(0.5)
    require("soft R3 equals uncoded BPSK at 7.3 dB",
            abs(Q(3 * sqrt(es) / (sqrt(3) * sigma)) - Q(sqrt(2 * ebn0))) < 1e-15)
    hard_db = solve_db(hard_r3, 1e-5)
    soft_db = solve_db(bpsk_ber, 1e-5)
    check("hard R3 reaches 1e-5 at (dB)", hard_db, "11.029")
    check("soft R3 reaches 1e-5 at (dB)", soft_db, "9.588")
    check("hard-decision loss (dB)", hard_db - soft_db, "1.441")
    check("loss over ten, the exponent", (hard_db - soft_db) / 10, "0.1441")
    check("hard-decision loss as a power ratio", undb(hard_db - soft_db), "1.393")
    check("extra power in percent", 100 * (undb(hard_db - soft_db) - 1), "39")
    check("figure caption: hard crossing", hard_db, "11.03")
    check("figure caption: soft crossing", soft_db, "9.59")
    check("figure caption: loss", hard_db - soft_db, "1.44")

    print("\nLikelihoods and log-likelihood ratios")
    a, sig = 1.0, 0.5
    c = 2 * a / sig**2
    check("2a / sigma^2 at a=1, sigma=0.5", c, "8")
    check("likelihood ratio at y=0.05", exp(c * 0.05), "1.49")
    check("likelihood ratio at y=1.9", exp(c * 1.9), "3.99e6")
    check("LLR at y=0.05", c * 0.05, "0.4")
    check("LLR at y=1.9", c * 1.9, "15.2")
    check("P(0 | y=0.05)", prob_of_llr(c * 0.05), "0.5987")
    check("P(0 | y=1.9)", prob_of_llr(c * 1.9), "0.99999975")
    p_hard = Q(a / sig)
    check("hard decision: p = Q(2)", p_hard, "0.02275")
    hard_llr = log((1 - p_hard) / p_hard)
    check("hard decision: LLR of a received 0", hard_llr, "3.760")
    check("figure: hard LLR", hard_llr, "3.76")
    check("soft and hard LLRs agree at |y|", hard_llr / c, "0.470")
    check("figure: crossing", hard_llr / c, "0.47")
    check("BSC(0.1) channel LLR", log(0.9 / 0.1), "2.197")
    check("Z(0.1) LLR of a received 0", log(1 / 0.1), "2.303")
    # The Z entry from the definition: P(0|0)/P(0|1) = 1/p.
    check("Z(0.1) LLR from transition probabilities", log(1.0 / 0.1), "2.303")
    c8 = 2 * 1.0 / 0.8**2
    ys = (0.9, -0.2, -0.3)
    llrs = [c8 * y for y in ys]
    check("2a / sigma^2 at sigma=0.8", c8, "3.125")
    check("LLR of 0.9", llrs[0], "2.8125")
    check("LLR of -0.2", llrs[1], "-0.625")
    check("LLR of -0.3", llrs[2], "-0.9375")
    check("sum of LLRs", sum(llrs), "1.25")
    check("P(0 | three copies)", prob_of_llr(sum(llrs)), "0.7773")
    check("voltage sum", sum(ys), "0.4")
    require("hard majority votes 1", sum(1 for y in ys if y < 0) >= 2)

    print("\nBursty channels")
    pi_g, pi_b = ge_stationary()
    check("pi_B", pi_b, "0.1")
    check("pi_G", pi_g, "0.9")
    check("convergence factor 1 - b - g", 1 - GE_B - GE_G, "0.9")
    check("0.9^50", (1 - GE_B - GE_G) ** 50, "0.00515")
    fixed = (1 - GE_B - GE_G) * pi_b + GE_B
    require("b/(b+g) is a fixed point of the recursion", abs(fixed - pi_b) < 1e-15)
    x = 0.7
    for _ in range(400):
        x = (1 - GE_B - GE_G) * x + GE_B
    require("the recursion converges to pi_B from 0.7", abs(x - pi_b) < 1e-12)
    avg = pi_g * GE_PG + pi_b * GE_PB
    check("average error rate", avg, "0.0109")
    check("mean bad spell 1/g", 1 / GE_G, "11.11")
    check("mean good spell 1/b", 1 / GE_B, "100")
    tail = sum((1 - GE_G) ** (k - 1) for k in range(1, 5000))
    check("tail sum of (1-g)^(k-1)", tail, "11.11")
    t1, t2 = pi_g * GE_PG, (1 - GE_B) * GE_PG + GE_B * GE_PB
    t3, t4 = pi_b * GE_PB, GE_G * GE_PG + (1 - GE_G) * GE_PB
    check("pi_G p_G", t1, "0.0009")
    check("(1-b) p_G + b p_B", t2, "0.00199")
    check("pi_B p_B", t3, "0.01")
    check("g p_G + (1-g) p_B", t4, "0.09109")
    pair = t1 * t2 + t3 * t4
    check("P(two consecutive errors)", pair, "0.000912691")
    check("memoryless value 0.0109^2", avg**2, "0.00011881")
    check("ratio to memoryless", pair / avg**2, "7.68")
    check("P(next bit wrong | error)", pair / avg, "0.0837")
    check("that over the average rate, nearly eight", (pair / avg) / avg, "7.7")
    back = ge_r3_back_to_back()
    spread = majority3(avg)
    check("R3 copies back to back, failure", back, "0.0024927")
    check("R3 copies spread apart, failure", spread, "0.00035384")
    check("back to back over spread", back / spread, "7.04")
    states, errors = ge_simulate()
    bsc_errors = bsc_simulate()
    print(f"  (figure data, not stated in the text: GE seed {GE_SEED} has {sum(errors)} errors and "
          f"{len(runs(states))} bad spells in {BURST_BITS} bits, BSC seed {BSC_SEED} has {sum(bsc_errors)} errors)")

    print("\nFault models for memory")
    check("C(72, 2)", comb(72, 2), "2556")
    word = comb(72, 2) * 1e-9**2
    check("SEC word failure bound per interval", word, "2.556e-15")
    check("across 2^30 words", word * 2**30, "2.744e-6")

    print("\nWorked examples")
    check("parity: no erasure in nine bits", 0.9**9, "0.38742")
    check("parity: exactly one erasure", 9 * 0.1 * 0.9**8, "0.38742")
    check("parity: recoverable", 0.9**9 + 9 * 0.1 * 0.9**8, "0.77484")
    data = "10110010"
    check("parity bit of 10110010", data.count("1") % 2, "0")
    known = data[:4] + data[5:] + "0"
    check("recovered fifth bit", known.count("1") % 2, "0")
    require("recovered bit equals the sent bit", int(data[4]) == known.count("1") % 2)
    x6 = Q_inv(1e-6)
    check("Q(x) = 1e-6: x", x6, "4.7534")
    check("Q(x) = 1e-6: Eb/N0", x6**2 / 2, "11.2975")
    check("Q(x) = 1e-6: Eb/N0 in dB", db(x6**2 / 2), "10.530")
    halved_db = db(x6**2 / 2) - db(2)
    check("halved power (dB)", halved_db, "7.520")
    check("halved power, ratio", undb(halved_db), "5.6488")
    check("halved power, P_b", Q(sqrt(2 * undb(halved_db))), "3.881e-4")
    check("that over the target, almost four hundred", Q(sqrt(2 * undb(halved_db))) / 1e-6, "388")
    check("figure caption: halved power (dB)", halved_db, "7.52")
    lnr = log(0.1 / 0.9)
    tau = (sig**2 / (2 * a)) * lnr
    check("ln(0.1/0.9)", lnr, "-2.1972")
    check("sigma^2 / 2a", sig**2 / (2 * a), "0.125")
    check("MAP threshold tau", tau, "-0.2747")
    check("figure: MAP threshold", tau, "-0.275")
    check("(a - tau) / sigma", (a - tau) / sig, "2.5493")
    check("P(error | 0 sent)", Q((a - tau) / sig), "0.005397")
    check("(a + tau) / sigma", (a + tau) / sig, "1.4507")
    check("P(error | 1 sent)", Q((a + tau) / sig), "0.07343")
    map_err = 0.9 * Q((a - tau) / sig) + 0.1 * Q((a + tau) / sig)
    check("MAP total error", map_err, "0.01220")
    check("ML error", Q(2), "0.02275")
    require("the weighted curves cross at tau",
            abs(0.9 * gauss_pdf(tau, 1, sig) - 0.1 * gauss_pdf(tau, -1, sig)) < 1e-12)
    check("peak of 0.9 f(y | +1)", 0.9 * gauss_pdf(1, 1, sig), "0.7181")
    bsc_part, awgn_part = -log(9), c * 0.3
    check("BSC(0.1) received 1", bsc_part, "-2.197")
    check("AWGN y = 0.3", awgn_part, "2.4")
    check("total LLR", awgn_part + bsc_part, "0.203")
    check("P(b = 0)", prob_of_llr(awgn_part + bsc_part), "0.5505")

    print("\nExercises")
    qz, pz = 0.9, 0.1
    check("Z prior: qp", qz * pz, "0.09")
    check("Z prior: (1-q) + qp", (1 - qz) + qz * pz, "0.19")
    check("Z prior: P(X=1 | Y=0)", qz * pz / ((1 - qz) + qz * pz), "0.47368")
    require("alarm example: q = 0.05 below p = 0.1, so a received 1 decodes as 0",
            (1 - 0.1) * 0.05 < 0.1 * (1 - 0.05))
    check("7 dB as a ratio", undb(7), "5.012")
    check("-3 dB as a ratio", undb(-3), "0.5012")
    check("13 dB as a ratio", undb(13), "19.95")
    check("ratio 50 in dB", db(50), "16.99")
    s2 = (0.6, 0.8)
    norm2 = sqrt(sum(v * v for v in s2))
    check("two-sample matched: Q(|s| / sigma)", Q(norm2 / 0.5), "0.02275")
    check("second sample alone: Q(1.6)", Q(0.8 / 0.5), "0.05480")
    check("plain sum: signal", sum(s2), "1.4")
    check("plain sum: signal over noise", sum(s2) / (0.5 * sqrt(2)), "1.9799")
    check("plain sum: error", Q(sum(s2) / (0.5 * sqrt(2))), "0.02386")
    l99 = log(0.99 / 0.01)
    check("BSC(0.01) LLR", l99, "4.595")
    check("received 0, 0, 1: LLR sum", l99 + l99 - l99, "4.595")
    check("received 0, 0, 1: posterior", prob_of_llr(l99), "0.99")
    check("Gilbert: average rate", ge_stationary()[1] * GE_PB, "0.01")
    gil_pair = ge_stationary()[1] * GE_PB * (1 - GE_G) * GE_PB
    check("Gilbert: two consecutive errors", gil_pair, "0.00091")
    check("Gilbert: ratio to memoryless", gil_pair / 0.01**2, "9.1")
    check("Gilbert: P(next | error)", gil_pair / 0.01, "0.091")
    es7 = undb(7) / 3
    check("7 dB: per-copy Es/N0", es7, "1.6706")
    check("7 dB: 2 Es/N0", 2 * es7, "3.3412")
    check("7 dB: p per copy", r3_copy_p(7), "0.03378")
    check("7 dB: hard R3", hard_r3(7), "0.003347")
    check("7 dB: Eb/N0 ratio", undb(7), "5.0119")
    check("7 dB: soft R3", bpsk_ber(7), "7.727e-4")
    check("7 dB: hard over soft", hard_r3(7) / bpsk_ber(7), "4.33")
    words = list(product((0, 1), repeat=8))
    corrupted = sum(1 for w in words if w[3] == 1) / len(words)
    check("stuck-at-0 cell: share of uniform 8-bit words it corrupts", corrupted, "0.5")

    print()
    print("All chapter 3 numbers reproduced.")


if __name__ == "__main__":
    main()
