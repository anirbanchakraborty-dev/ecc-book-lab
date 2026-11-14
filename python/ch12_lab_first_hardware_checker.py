"""Reproduce every number in chapter 12 of Error Correction Coding, "Lab: A First Hardware Checker".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch12_lab_first_hardware_checker.py

The module is the software model of the parity bench in hdl/parity/. It holds the
generator, the channel and the checker as a few lines of integer arithmetic, the
gate and depth counts of the recursive XOR tree in hdl/parity/xor_tree.sv, the
closed forms for the alarm, undetected-error and coverage rates of single parity
on a binary symmetric channel, and the Wilson interval, zero-failure bound and
trial count that chapter 11 builds.

Run as a script, it replays the seeded streams that tb/parity/test_parity.py drives
through the Verilator simulation, recounts everything the cocotb tests log, and
prints each value computed here beside the value the chapter states. It stops
with exit status 1 at the first pair that disagrees by more than half a unit in
the last digit the chapter prints, and at the first count that differs at all.
Every random number comes from random.Random(seed).random(), the one method whose
sequence for a given seed the Python documentation guarantees across Python
versions, and every draw is compared with a fixed threshold, so the counts
reproduce on every machine. The gate counts and depths Yosys reported, which this
script cannot rerun without Yosys, are stored below as measured by
synth/synth_parity.py and checked against the formulas the chapter proves.
"""

import random
import sys
from math import ceil, comb, log2, sqrt

SEED = 12                    # the channel generator; the bench restarts it for every p
SOURCE_OFFSET = 10**9        # the data generator is seeded with SEED + SOURCE_OFFSET
Z95 = 1.96
P_VALUES = (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
WORDS_PER_P = 100_000
WAVE_P = 0.2
WAVE_WORDS = 12
STEP_NS = 10

# What the cocotb run logged for K = 7, seed 12, 100,000 words at each p:
# p: (damaged, flagged, slipped, coverage, Wilson low, Wilson high, formula, verdict)
MEASURED = {
    0.01: (7640, 7379, 261, "0.96584", "0.96152", "0.96968", "0.96587", "inside"),
    0.02: (14829, 13806, 1023, "0.93101", "0.92682", "0.93498", "0.93345", "inside"),
    0.05: (33558, 28362, 5196, "0.84516", "0.84125", "0.84899", "0.84606", "inside"),
    0.1: (56906, 41303, 15603, "0.72581", "0.72213", "0.72946", "0.73062", "OUTSIDE"),
    0.2: (83373, 49271, 34102, "0.59097", "0.58763", "0.59430", "0.59071", "inside"),
    0.3: (94208, 49611, 44597, "0.52661", "0.52342", "0.52980", "0.53024", "OUTSIDE"),
    0.4: (98315, 49962, 48353, "0.50818", "0.50506", "0.51131", "0.50854", "inside"),
    0.5: (99641, 50215, 49426, "0.50396", "0.50085", "0.50706", "0.50196", "inside"),
}

# The twelve words of the waveform window, as test_waveform logs them:
# (time in ns, data, error mask, received word, error flag)
WAVEFORM = (
    (0, "0000101", "00011000", "00011101", 0),
    (10, "0101010", "01010000", "11111010", 0),
    (20, "1101100", "10000101", "11101001", 1),
    (30, "1011110", "01001001", "10010111", 1),
    (40, "1010000", "01101001", "00111001", 0),
    (50, "0110001", "00100000", "10010001", 1),
    (60, "0011010", "01000000", "11011010", 1),
    (70, "0110101", "00010000", "00100101", 1),
    (80, "0001110", "01010000", "11011110", 0),
    (90, "0001110", "00010001", "10011111", 0),
    (100, "1011110", "10010000", "01001110", 0),
    (110, "0100001", "10000100", "10100101", 0),
)

# What synth/synth_parity.py measured with Yosys 0.71:
# (block, n, width W of its tree, gates, XOR cells, XNOR cells, longest path in cells)
YOSYS = (
    ("parity_gen", 8, 7, 6, 0, 6, 3),
    ("parity_gen", 16, 15, 14, 4, 10, 4),
    ("parity_gen", 32, 31, 30, 8, 22, 5),
    ("parity_gen", 64, 63, 62, 20, 42, 6),
    ("parity_check", 8, 8, 7, 3, 4, 3),
    ("parity_check", 16, 16, 15, 5, 10, 4),
    ("parity_check", 32, 32, 31, 11, 20, 5),
    ("parity_check", 64, 64, 63, 21, 42, 6),
)


# Helpers ------------------------------------------------------------------------

def half_unit(stated: str) -> float:
    text = stated.replace(",", "").rstrip("%")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", "").rstrip("%"))
    if stated.endswith("%"):
        computed = 100 * float(computed)
    print(f"  {label:<64} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<64} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def holds(label: str, condition: bool) -> None:
    print(f"  {label:<64} {'yes' if condition else 'NO'}")
    if not condition:
        print(f"MISMATCH: '{label}' does not hold")
        sys.exit(1)


# The bench in software ----------------------------------------------------------

def weight(x: int) -> int:
    """The number of 1 bits in x."""
    return bin(x).count("1")


def codeword(data: int, k: int) -> int:
    """What parity_gen.sv sends: the k data bits with the even-parity bit on top, at bit k."""
    return data | ((weight(data) & 1) << k)


def flag(word: int) -> int:
    """What parity_check.sv raises: 1 when the word has odd weight."""
    return weight(word) & 1


def or_root_flag(word: int, n: int) -> int:
    """A broken checker whose top gate is OR where it should be XOR."""
    low = n // 2
    return (weight(word & ((1 << low) - 1)) & 1) | (weight(word >> low) & 1)


def tree_gates(w: int) -> int:
    """Gates in xor_tree.sv over w inputs, counted the way the module builds them."""
    return 0 if w == 1 else tree_gates(w // 2) + tree_gates(w - w // 2) + 1


def tree_depth(w: int) -> int:
    """Gates on the deepest input-to-output path of the same tree."""
    return 0 if w == 1 else 1 + max(tree_depth(w // 2), tree_depth(w - w // 2))


def bsc_stream(p: float, k: int, seed: int = SEED):
    """The (data, mask) pairs test_parity.py drives, from the same two generators."""
    source = random.Random(seed + SOURCE_OFFSET)
    channel = random.Random(seed)
    while True:
        data = sum(int(source.random() < 0.5) << i for i in range(k))
        mask = sum(int(channel.random() < p) << i for i in range(k + 1))
        yield data, mask


def run_bsc(p: float, k: int = 7, words: int = WORDS_PER_P, seed: int = SEED) -> tuple:
    """Damaged, flagged and slipped counts for one crossover probability."""
    stream = bsc_stream(p, k, seed)
    damaged = flagged = slipped = 0
    for _ in range(words):
        data, mask = next(stream)
        f = flag(codeword(data, k) ^ mask)
        flagged += f
        if mask:
            damaged += 1
            slipped += 1 - f
    return damaged, flagged, slipped


def exhaustive(k: int, checker=None) -> tuple:
    """Every data word against every mask: vectors, flagged, slipped, checker disagreements."""
    n = k + 1
    flagged = slipped = wrong = 0
    first_wrong = None
    for data in range(2**k):
        for mask in range(2**n):
            received = codeword(data, k) ^ mask
            f = flag(received)
            flagged += f
            slipped += int(mask != 0 and f == 0)
            if checker is not None and checker(received, n) != f:
                wrong += 1
                if first_wrong is None:
                    first_wrong = (data, mask)
    return 2**k * 2**n, flagged, slipped, wrong, first_wrong


# Statistics and closed forms ----------------------------------------------------

def wilson(k: int, n: int, z: float = Z95) -> tuple:
    ph = k / n
    center = ph + z * z / (2 * n)
    half = z * sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    scale = 1 + z * z / n
    return (center - half) / scale, (center + half) / scale


def zero_failure_bound(n: int, confidence: float = 0.95) -> float:
    return 1 - (1 - confidence) ** (1 / n)


def trials_needed(p: float, r: float, z: float = Z95) -> int:
    return ceil((z / r) ** 2 * (1 - p) / p)


def p_alarm(p: float, n: int) -> float:
    """Probability that an odd number of the n bits flip, so the checker raises its flag."""
    return (1 - (1 - 2 * p) ** n) / 2


def p_damaged(p: float, n: int) -> float:
    return 1 - (1 - p) ** n


def p_undetected(p: float, n: int) -> float:
    """Probability of an even, nonzero number of flips."""
    return (1 + (1 - 2 * p) ** n) / 2 - (1 - p) ** n


def coverage(p: float, n: int) -> float:
    """Fraction of damaged words the checker flags."""
    return p_alarm(p, n) / p_damaged(p, n)


def z_score(caught: int, damaged: int, c: float) -> float:
    return (caught / damaged - c) / sqrt(c * (1 - c) / damaged)


def pearson(a: list, b: list) -> float:
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return cov / sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))


# Checks -------------------------------------------------------------------------

def check_tree() -> None:
    print("The XOR tree")
    holds("gates = W - 1 for every W from 1 to 1,024",
          all(tree_gates(w) == w - 1 for w in range(1, 1025)))
    holds("depth = ceil(log2 W) for every W from 1 to 1,024",
          all(tree_depth(w) == ceil(log2(w)) for w in range(1, 1025)))
    same("generator tree for K = 7: gates", tree_gates(7), 6)
    same("generator tree for K = 7: depth", tree_depth(7), 3)
    same("checker tree for n = 8: gates", tree_gates(8), 7)
    same("checker tree for n = 8: depth", tree_depth(8), 3)
    same("W = 13 splits into", (13 // 2, 13 - 13 // 2), (6, 7))
    same("W = 13: gates", tree_gates(13), 12)
    same("W = 13: depth", tree_depth(13), 4)
    same("W = 6: depth", tree_depth(6), 3)
    same("W = 7: depth", tree_depth(7), 3)
    same("n = 64 checker: gates", tree_gates(64), 63)
    same("n = 64 checker: depth", tree_depth(64), 6)
    same("W = 2, 3, 4: depth", (tree_depth(2), tree_depth(3), tree_depth(4)), (1, 2, 2))
    d = [(0b1101100 >> i) & 1 for i in range(7)]
    same("1101100: d6 down to d0", d[::-1], [1, 1, 0, 1, 1, 0, 0])
    same("low group: d1 xor d2, then d0 xor that", (d[1] ^ d[2], d[0] ^ d[1] ^ d[2]), (1, 1))
    same("high group: d3 xor d4, d5 xor d6, then both", (d[3] ^ d[4], d[5] ^ d[6], d[3] ^ d[4] ^ d[5] ^ d[6]), (1, 0, 1))
    same("root: parity bit of 1101100", d[0] ^ d[1] ^ d[2] ^ d[3] ^ d[4] ^ d[5] ^ d[6], 0)
    print("Yosys netlists against the formulas")
    for block, n, w, gates, xor, xnor, depth in YOSYS:
        holds(f"{block} n = {n}: {gates} gates = W - 1, depth {depth} = ceil(log2 {w})",
              gates == w - 1 == tree_gates(w) and xor + xnor == gates and depth == ceil(log2(w)) == tree_depth(w))
        holds(f"{block} n = {n}: {xnor} XNOR cells, an even number", xnor % 2 == 0)


def check_exhaustive() -> None:
    print("Every error pattern")
    for k, vectors, flagged, slipped in ((7, 32768, 16384, 16256), (3, 128, 64, 56)):
        v, f, s, _, _ = exhaustive(k)
        same(f"k = {k}: vectors", v, vectors)
        same(f"k = {k}: flagged", f, flagged)
        same(f"k = {k}: slipped", s, slipped)
        same(f"k = {k}: flagged per data word", f // 2**k, 2**k)
        same(f"k = {k}: slipped per data word", s // 2**k, 2**k - 1)
    holds("every k = 7 vector: flag = weight(mask) mod 2",
          all(flag(codeword(d, 7) ^ m) == weight(m) % 2 for d in range(128) for m in range(256)))
    holds("every clean k = 7 and k = 3 word is unflagged",
          all(flag(codeword(d, 7)) == 0 for d in range(128)) and all(flag(codeword(d, 3)) == 0 for d in range(8)))
    check("share of nonzero masks that slip, 127/255", 127 / 255, "0.498")
    same("clean words sent by test_clean_words for k = 7", 2**7, 128)
    same("vectors for 32-bit codewords, 2^(2k + 1) with k = 31", 2 ** (2 * 31 + 1), 2**63)
    v, _, _, wrong, first = exhaustive(7, or_root_flag)
    same("OR at the root: vectors the mutant gets wrong", wrong, 8192)
    same("OR at the root: first failing (data, mask)", first, (0, 17))
    same("mask 17 in binary", f"{17:08b}", "00010001")
    same("words of 8 bits with both halves odd", sum(1 for w in range(256) if weight(w & 15) % 2 and weight(w >> 4) % 2), 64)
    check("share of vectors the mutant gets wrong", wrong / v, "0.25")


def check_waveform() -> None:
    print("The waveform window, p = 0.2, seed 12")
    stream = bsc_stream(WAVE_P, 7)
    caught = slipped = 0
    weights = []
    for t, data_s, mask_s, recv_s, f_s in WAVEFORM:
        data, mask = next(stream)
        received = codeword(data, 7) ^ mask
        same(f"t = {t:>3} ns", (f"{data:07b}", f"{mask:08b}", f"{received:08b}", flag(received)),
             (data_s, mask_s, recv_s, f_s))
        weights.append(weight(mask))
        caught += flag(received)
        slipped += int(mask != 0 and flag(received) == 0)
    same("mask weights in order", weights, [2, 2, 3, 3, 4, 1, 1, 1, 2, 2, 2, 2])
    same("caught in the window", caught, 5)
    same("slipped in the window", slipped, 7)
    same("window length in ns", WAVE_WORDS * STEP_NS, 120)
    first = codeword(0b1101100, 7)
    same("t = 20 ns: codeword of 1101100", f"{first:08b}", "01101100")
    same("t = 20 ns: weight of the received word", weight(first ^ 0b10000101), 5)
    same("t = 20 ns: bits the mask flips", [i for i in range(7, -1, -1) if 0b10000101 >> i & 1], [7, 2, 0])
    same("times at which error is high", [row[0] for row in WAVEFORM if row[4]], [20, 30, 50, 60, 70])
    slip = p_undetected(WAVE_P, 8)
    check("probability a word slips at p = 0.2", slip, "0.341")
    check("slips expected among the 12 words", WAVE_WORDS * slip, "4.1")
    check("seven or more slips among 12", sum(comb(12, k) * slip**k * (1 - slip) ** (12 - k) for k in range(7, 13)), "7%")


def check_fault_injection() -> None:
    print("Fault injection, 100,000 words at each p, seed 12 (replaying the bench)")
    damaged_counts = []
    for p in P_VALUES:
        damaged, flagged, slipped = run_bsc(p)
        d_s, f_s, s_s, cov_s, lo_s, hi_s, form_s, verdict_s = MEASURED[p]
        caught = damaged - slipped
        lo, hi = wilson(caught, damaged)
        c = coverage(p, 8)
        same(f"p = {p}: damaged, flagged, slipped", (damaged, flagged, slipped), (d_s, f_s, s_s))
        holds(f"p = {p}: every flagged word was damaged", flagged == caught)
        check(f"p = {p}: measured coverage", caught / damaged, cov_s)
        check(f"p = {p}: Wilson low", lo, lo_s)
        check(f"p = {p}: Wilson high", hi, hi_s)
        check(f"p = {p}: formula", c, form_s)
        same(f"p = {p}: formula inside the interval", "inside" if lo <= c <= hi else "OUTSIDE", verdict_s)
        damaged_counts.append(damaged)
    holds("damaged counts rise with p", all(a < b for a, b in zip(damaged_counts, damaged_counts[1:])))
    holds("every p = 0.1 mask sits inside the p = 0.3 mask, word by word",
          all((a[1] & b[1]) == a[1] for _, a, b in zip(range(WORDS_PER_P), bsc_stream(0.1, 7), bsc_stream(0.3, 7))))
    same("words sent through the channel in all", len(P_VALUES) * WORDS_PER_P, 800_000)
    same("simulated time at the end of the eight runs, ns", len(P_VALUES) * WORDS_PER_P * STEP_NS, 8_000_000)
    same("seed of the data generator", SEED + SOURCE_OFFSET, 1_000_000_012)
    check("misses expected among eight 95% intervals", 8 * 0.05, "0.4")
    check("z at p = 0.1", z_score(41303, 56906, coverage(0.1, 8)), "-2.59")
    check("z at p = 0.3", z_score(49611, 94208, coverage(0.3, 8)), "-2.23")
    check("two or more misses of 8 if they were independent", 1 - 0.95**8 - 8 * 0.05 * 0.95**7, "0.057")
    check("that is one chance in", 1 / (1 - 0.95**8 - 8 * 0.05 * 0.95**7), "17")
    print("How strongly the shared draws couple the runs, seeds 0 to 399, 5,000 words each")
    coupled = (0.01, 0.02, 0.1, 0.3)
    zs_by_p = {p: [] for p in coupled}
    for seed in range(400):
        for p in coupled:
            damaged, _, slipped = run_bsc(p, words=5_000, seed=seed)
            zs_by_p[p].append(z_score(damaged - slipped, damaged, coverage(p, 8)))
    check("correlation of z at p = 0.01 and p = 0.02", pearson(zs_by_p[0.01], zs_by_p[0.02]), "0.43")
    check("correlation of z at p = 0.1 and p = 0.3", pearson(zs_by_p[0.1], zs_by_p[0.3]), "0.02")
    print("The same experiment at p = 0.1 over seeds 0 to 299, 20,000 words each")
    misses = 0
    zs = []
    c = coverage(0.1, 8)
    for seed in range(300):
        damaged, _, slipped = run_bsc(0.1, words=20_000, seed=seed)
        lo, hi = wilson(damaged - slipped, damaged)
        misses += int(not lo <= c <= hi)
        zs.append(z_score(damaged - slipped, damaged, c))
    mean = sum(zs) / len(zs)
    sd = sqrt(sum((z - mean) ** 2 for z in zs) / (len(zs) - 1))
    same("intervals that miss the formula", misses, 12)
    check("share that miss", misses / 300, "4%")
    lo, hi = wilson(misses, 300)
    check("Wilson interval on 12 of 300: low", lo, "0.023")
    check("Wilson interval on 12 of 300: high", hi, "0.069")
    check("mean z over 300 seeds", mean, "-0.09")
    check("standard deviation of z over 300 seeds", sd, "0.97")


def check_closed_forms() -> None:
    print("Closed forms at n = 8, p = 0.01")
    n, p = 8, 0.01
    check("(1 - 2p)^n", (1 - 2 * p) ** n, "0.850763")
    check("probability of an alarm", p_alarm(p, n), "0.074618")
    check("probability the word is damaged", p_damaged(p, n), "0.077255")
    check("probability of an undetected error", p_undetected(p, n), "0.0026368")
    check("alarm + undetected = damaged", p_alarm(p, n) + p_undetected(p, n), "0.077255")
    check("coverage", coverage(p, n), "0.96587")
    check("approximation 1 - (n - 1)p/2", 1 - (n - 1) * p / 2, "0.965")
    check("coverage at p = 0.5, 128/255", coverage(0.5, n), "0.50196")
    holds("coverage at p = 0.5 is exactly 128/255", abs(coverage(0.5, n) - 128 / 255) < 1e-15)
    check("expected damaged words in 100,000 at p = 0.01", 100_000 * p_damaged(p, n), "7,725.5")
    check("(1 - p)^n", (1 - p) ** n, "0.922745")
    check("second flip against first, (n - 1)p/2", (n - 1) * p / 2, "3.5%")
    two_or_more = 1 - (1 - p) ** n - n * p * (1 - p) ** (n - 1)
    check("share of damaged words that carry a second flip", two_or_more / p_damaged(p, n), "3.5%")
    sd = sqrt(100_000 * p_damaged(p, n) * (1 - p_damaged(p, n)))
    check("standard deviation of the damaged count", sd, "84.4")
    check("7,640 damaged in standard deviations from 7,725.5", (7640 - 100_000 * p_damaged(p, n)) / sd, "-1.01")
    print("The Wilson interval for 7,379 of 7,640, step by step")
    c_hat, D, z2 = 7379 / 7640, 7640, Z95**2
    check("estimate", c_hat, "0.965838")
    check("z^2", z2, "3.8416")
    check("center before scaling", c_hat + z2 / (2 * D), "0.966089")
    check("spread", Z95 * sqrt(c_hat * (1 - c_hat) / D + z2 / (4 * D * D)), "0.0040809")
    check("scale", 1 + z2 / D, "1.000503")
    print("Precision of the estimates")
    check("Wilson half-width scale at p = 0.01, 1.96 sqrt(c(1-c)/7,640)",
          Z95 * sqrt(coverage(p, n) * (1 - coverage(p, n)) / 7640), "0.0041")
    check("measured undetected rate at p = 0.01, 261/100,000", 261 / 100_000, "0.00261")
    lo, hi = wilson(261, 100_000)
    check("its Wilson low", lo, "0.00231")
    check("its Wilson high", hi, "0.00295")
    holds("the closed form 0.0026368 lies inside", lo <= p_undetected(p, n) <= hi)
    P = p_undetected(p, n)
    check("(1 - P)/P for the undetected rate", (1 - P) / P, "378.245")
    check("(1.96/0.1)^2", (Z95 / 0.1) ** 2, "384.16")
    same("words for 10% precision on the undetected rate", trials_needed(P, 0.1), 145_307)
    check("precision 100,000 words give, 1.96 sqrt((1 - P)/(N P))", Z95 * sqrt((1 - P) / (100_000 * P)), "12%")
    check("extra words 10% precision needs beyond 100,000", trials_needed(P, 0.1) / 100_000 - 1, "45%")


def check_exercise_answers() -> None:
    print("Exercise answers")
    n, p = 8, 0.001
    P = p_undetected(p, n)
    check("undetected rate at p = 0.001", P, "2.78e-5")
    check("expected slips in 20,000 words", 20_000 * P, "0.56")
    check("probability of no slip in 20,000 words", (1 - P) ** 20_000, "0.57")
    check("zero-failure bound at N = 20,000", zero_failure_bound(20_000), "1.4978e-4")
    check("bound over the true rate", zero_failure_bound(20_000) / P, "5.4")
    same("odd masks among the 16 of length 4", sum(weight(m) & 1 for m in range(16)), 8)
    same("even nonzero masks among the 16 of length 4", sum(1 for m in range(1, 16) if weight(m) % 2 == 0), 7)
    check("coverage at n = 4, p = 0.5, 8/15", coverage(0.5, 4), "0.5333")
    same("three XNOR gates invert the parity, 3 mod 2", 3 % 2, 1)


def main() -> None:
    check_tree()
    check_exhaustive()
    check_waveform()
    check_closed_forms()
    check_fault_injection()
    check_exercise_answers()
    print("\nAll chapter 12 numbers reproduced.")


if __name__ == "__main__":
    main()
