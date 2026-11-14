"""Reproduce every number in chapter 11 of Error Correction Coding, "Project: A Channel Simulator".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch11_project_channel_simulator.py

The module is the chapter's channel simulator and its reference model at once.
Importing it gives the three channels (binary symmetric, binary erasure and BPSK
over Gaussian noise, with the Gaussian samples made by the Box-Muller transform),
the repetition and single parity check codes with their decoders, the Monte Carlo
estimator with its Wald and Wilson intervals, the zero-failure bound, the number
of trials a target precision needs, and the closed forms each simulation is
checked against.

Run as a script, it repeats every seeded experiment the chapter reports and
prints each value computed here beside the value the chapter states. It stops
with exit status 1 at the first pair that disagrees by more than half a unit in
the last digit the chapter prints, and at the first count that differs at all.
Every random number comes from random.Random(seed).random(), the one method whose
sequence for a given seed the Python documentation guarantees across Python
versions. The counts on the discrete channels compare those uniforms with fixed
thresholds, so in practice they reproduce on every machine. The Gaussian counts
also pass the uniforms through log, cos and sin from the platform's math library
and scale them by a sigma computed with pow, so they reproduce on every machine
unless two libraries round one of those functions differently in the last bit.
"""

import random
import sys
from math import ceil, comb, cos, degrees, erfc, exp, log, pi, sin, sqrt

Z95 = 1.96
SOURCE_OFFSET = 10**9


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


# Channels -----------------------------------------------------------------------

def bsc(word: list, p: float, rng: random.Random) -> list:
    """Binary symmetric channel: each bit flips when a uniform draw falls below p."""
    return [b ^ int(rng.random() < p) for b in word]


def bec(word: list, eps: float, rng: random.Random) -> list:
    """Binary erasure channel: each bit becomes None, an erasure, when a draw falls below eps."""
    return [None if rng.random() < eps else b for b in word]


class GaussianSource:
    """Standard normal samples made from uniform draws by the Box-Muller transform.

    Each pair of uniform draws gives two independent samples, so the second is
    kept and handed out on the next call.
    """

    def __init__(self, rng: random.Random):
        self.rng = rng
        self.spare = None

    def next(self) -> float:
        if self.spare is not None:
            z, self.spare = self.spare, None
            return z
        u1 = 1.0 - self.rng.random()  # in (0, 1], so log(u1) is finite
        u2 = self.rng.random()        # in [0, 1)
        r = sqrt(-2.0 * log(u1))
        theta = 2.0 * pi * u2
        self.spare = r * sin(theta)
        return r * cos(theta)


def awgn_sigma(ebn0_db: float, rate: float) -> float:
    """Noise standard deviation for pulses of energy Es = 1 carrying `rate` data bits each."""
    ebn0 = 10 ** (ebn0_db / 10)
    return sqrt(1.0 / (2.0 * rate * ebn0))


def bpsk_awgn(word: list, sigma: float, noise: GaussianSource) -> list:
    """Send 0 as +1 and 1 as -1, then add Gaussian noise of standard deviation sigma."""
    return [(1.0 - 2.0 * b) + sigma * noise.next() for b in word]


# Codes --------------------------------------------------------------------------

def repetition_encode(bit: int, n: int) -> list:
    return [bit] * n


def majority(word: list) -> int:
    return int(sum(word) > len(word) // 2)


def hard_repetition(received: list) -> int:
    """Decide each pulse by its sign, then take the majority."""
    return majority([int(y < 0) for y in received])


def soft_repetition(received: list) -> int:
    """Add the received values and decide by the sign of the sum."""
    return int(sum(received) < 0)


def repetition_erasures(received: list):
    """On the erasure channel any copy that arrives is correct. None when every copy is erased."""
    for b in received:
        if b is not None:
            return b
    return None


def parity_encode(data: list) -> list:
    return data + [sum(data) % 2]


def parity_ok(word: list) -> bool:
    return sum(word) % 2 == 0


def parity_fill(received: list):
    """Fill a single erasure from the parity rule. None when two or more positions are erased."""
    missing = [i for i, b in enumerate(received) if b is None]
    if len(missing) > 1:
        return None
    word = list(received)
    if missing:
        word[missing[0]] = sum(b for b in received if b is not None) % 2
    return word


# The simulation loop ------------------------------------------------------------

def count_failures(trials: int, seed: int, trial, gaussian: bool = False) -> int:
    """Run `trial` `trials` times and count the runs that return True.

    The message bits come from their own generator, seeded with seed + 10**9, so
    the channel sees the same uniform draws whatever the messages are.
    """
    source = random.Random(seed + SOURCE_OFFSET)
    channel = random.Random(seed)
    noise = GaussianSource(channel) if gaussian else channel
    failures = 0
    for _ in range(trials):
        failures += trial(source, noise)
    return failures


def random_bits(source: random.Random, k: int) -> list:
    return [int(source.random() < 0.5) for _ in range(k)]


def repetition_bsc_trial(n: int, p: float):
    def trial(source, channel):
        b = random_bits(source, 1)[0]
        return majority(bsc(repetition_encode(b, n), p, channel)) != b
    return trial


def repetition_bec_trial(n: int, eps: float):
    def trial(source, channel):
        b = random_bits(source, 1)[0]
        return repetition_erasures(bec(repetition_encode(b, n), eps, channel)) != b
    return trial


def parity_undetected_trial(n: int, p: float):
    """True when the received word passes the parity check but differs from the word sent."""
    def trial(source, channel):
        word = parity_encode(random_bits(source, n - 1))
        received = bsc(word, p, channel)
        return parity_ok(received) and received != word
    return trial


def parity_detected_trial(n: int, p: float):
    def trial(source, channel):
        word = parity_encode(random_bits(source, n - 1))
        return not parity_ok(bsc(word, p, channel))
    return trial


def parity_bec_trial(n: int, eps: float):
    def trial(source, channel):
        word = parity_encode(random_bits(source, n - 1))
        filled = parity_fill(bec(word, eps, channel))
        if filled is not None and filled != word:
            raise AssertionError("parity filled an erasure wrongly")
        return filled is None
    return trial


def uncoded_bpsk_trial(ebn0_db: float):
    sigma = awgn_sigma(ebn0_db, 1.0)

    def trial(source, noise):
        b = random_bits(source, 1)[0]
        return int(bpsk_awgn([b], sigma, noise)[0] < 0) != b
    return trial


def repetition_awgn_trial(n: int, ebn0_db: float, decoder):
    sigma = awgn_sigma(ebn0_db, 1.0 / n)

    def trial(source, noise):
        b = random_bits(source, 1)[0]
        return decoder(bpsk_awgn(repetition_encode(b, n), sigma, noise)) != b
    return trial


# Estimation ---------------------------------------------------------------------

def wald(k: int, n: int, z: float = Z95) -> tuple:
    ph = k / n
    half = z * sqrt(ph * (1 - ph) / n)
    return ph - half, ph + half


def wilson(k: int, n: int, z: float = Z95) -> tuple:
    ph = k / n
    center = ph + z * z / (2 * n)
    half = z * sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    scale = 1 + z * z / n
    return (center - half) / scale, (center + half) / scale


def exact_coverage(p: float, n: int, interval=wilson) -> float:
    """Probability that the interval after n trials contains the true rate p."""
    total = 0.0
    for k in range(n + 1):
        lo, hi = interval(k, n)
        if lo <= p <= hi:
            total += comb(n, k) * p**k * (1 - p) ** (n - k)
    return total


def zero_failure_bound(n: int, confidence: float = 0.95) -> float:
    """Largest rate that still shows no failure in n trials with probability 1 - confidence."""
    return 1 - (1 - confidence) ** (1 / n)


def trials_needed(p: float, r: float, z: float = Z95) -> int:
    """Trials for the interval's half-width to be at most r times p."""
    return ceil((z / r) ** 2 * (1 - p) / p)


# Closed forms -------------------------------------------------------------------

def q(x: float) -> float:
    return 0.5 * erfc(x / sqrt(2))


def repetition_bsc_exact(n: int, p: float) -> float:
    t = (n - 1) // 2
    return sum(comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(t + 1, n + 1))


def repetition_bec_exact(n: int, eps: float) -> float:
    return eps**n


def uncoded_bpsk_exact(ebn0_db: float) -> float:
    return q(sqrt(2 * 10 ** (ebn0_db / 10)))


def hard_repetition_awgn_exact(n: int, ebn0_db: float) -> float:
    return repetition_bsc_exact(n, q(sqrt(2 * 10 ** (ebn0_db / 10) / n)))


def soft_repetition_awgn_exact(n: int, ebn0_db: float) -> float:
    return uncoded_bpsk_exact(ebn0_db)


def parity_undetected_exact(n: int, p: float) -> float:
    return (1 + (1 - 2 * p) ** n) / 2 - (1 - p) ** n


def parity_detected_exact(n: int, p: float) -> float:
    return (1 - (1 - 2 * p) ** n) / 2


def parity_bec_exact(n: int, eps: float) -> float:
    return 1 - (1 - eps) ** n - n * eps * (1 - eps) ** (n - 1)


# The experiments ----------------------------------------------------------------
# Each entry: name, seed, trials, the trial factory, Gaussian or not, exact value.

def experiments() -> list:
    return [
        ("R3 on BSC(0.1)", 2026, 10_000, repetition_bsc_trial(3, 0.1), False, repetition_bsc_exact(3, 0.1)),
        ("R3 on BEC(0.2)", 12, 100_000, repetition_bec_trial(3, 0.2), False, repetition_bec_exact(3, 0.2)),
        ("parity (8,7) on BSC(0.01), undetected", 13, 100_000, parity_undetected_trial(8, 0.01), False,
         parity_undetected_exact(8, 0.01)),
        ("parity (8,7) on BEC(0.1), unfilled", 14, 20_000, parity_bec_trial(8, 0.1), False,
         parity_bec_exact(8, 0.1)),
        ("uncoded BPSK at 4 dB", 15, 100_000, uncoded_bpsk_trial(4.0), True, uncoded_bpsk_exact(4.0)),
        ("hard R3 at 4 dB", 16, 100_000, repetition_awgn_trial(3, 4.0, hard_repetition), True,
         hard_repetition_awgn_exact(3, 4.0)),
        ("soft R3 at 4 dB", 16, 100_000, repetition_awgn_trial(3, 4.0, soft_repetition), True,
         soft_repetition_awgn_exact(3, 4.0)),
    ]


def run_experiments() -> list:
    rows = []
    for name, seed, trials, trial, gaussian, exact in experiments():
        k = count_failures(trials, seed, trial, gaussian)
        lo, hi = wilson(k, trials)
        rows.append((name, seed, trials, k, lo, hi, exact))
    return rows


COVERAGE_RUNS = 1000
COVERAGE_LONG = 1000
COVERAGE_SHORT = 100
COVERAGE_EXACT = repetition_bsc_exact(3, 0.1)


def coverage_counts(trials: int, runs: int = COVERAGE_RUNS) -> list:
    """Failures of R3 on BSC(0.1) in each of `runs` runs, seeds 1 to runs."""
    trial = repetition_bsc_trial(3, 0.1)
    return [count_failures(trials, seed, trial) for seed in range(1, runs + 1)]


def box_muller_samples(count: int, seed: int) -> list:
    noise = GaussianSource(random.Random(seed))
    return [noise.next() for _ in range(count)]


BOX_MULLER_SEED = 17
BOX_MULLER_COUNT = 100_000

SWEEP_DB = list(range(0, 8))
SWEEP_TRIALS = 100_000


def sweep_seed(db: int) -> int:
    return 100 + db


def awgn_sweep() -> list:
    """Uncoded, hard R3 and soft R3 at 0 to 7 dB. Each point uses seed 100 + dB and 100,000 trials."""
    rows = []
    for db in SWEEP_DB:
        seed = sweep_seed(db)
        ku = count_failures(SWEEP_TRIALS, seed, uncoded_bpsk_trial(db), True)
        kh = count_failures(SWEEP_TRIALS, seed, repetition_awgn_trial(3, db, hard_repetition), True)
        ks = count_failures(SWEEP_TRIALS, seed, repetition_awgn_trial(3, db, soft_repetition), True)
        rows.append((db, ku, kh, ks))
    return rows


ZERO_SEED = 18
ZERO_TRIALS = 20_000


def paired_hard_soft(trials: int, seed: int, ebn0_db: float, n: int = 3) -> tuple:
    """Decode the same received words both ways and count the words both get wrong, only hard, only soft."""
    source = random.Random(seed + SOURCE_OFFSET)
    noise = GaussianSource(random.Random(seed))
    sigma = awgn_sigma(ebn0_db, 1.0 / n)
    both = hard_only = soft_only = 0
    for _ in range(trials):
        b = random_bits(source, 1)[0]
        received = bpsk_awgn(repetition_encode(b, n), sigma, noise)
        hard_wrong = hard_repetition(received) != b
        soft_wrong = soft_repetition(received) != b
        both += hard_wrong and soft_wrong
        hard_only += hard_wrong and not soft_wrong
        soft_only += soft_wrong and not hard_wrong
    return both, hard_only, soft_only


# Checks -------------------------------------------------------------------------

def check_channels_and_codes() -> None:
    print("Channels and decoders")
    rng = random.Random(1)
    same("BSC with p = 0 changes nothing", bsc([0, 1, 1, 0], 0.0, rng), [0, 1, 1, 0])
    same("BSC with p = 1 flips every bit", bsc([0, 1, 1, 0], 1.0, rng), [1, 0, 0, 1])
    same("BEC with eps = 1 erases every bit", bec([0, 1], 1.0, rng), [None, None])
    same("parity of 1011001", parity_encode([1, 0, 1, 1, 0, 0, 1]), [1, 0, 1, 1, 0, 0, 1, 0])
    same("parity check passes 10110010", parity_ok([1, 0, 1, 1, 0, 0, 1, 0]), True)
    same("parity check fails 10100010", parity_ok([1, 0, 1, 0, 0, 0, 1, 0]), False)
    same("parity check passes 10100011, two errors", parity_ok([1, 0, 1, 0, 0, 0, 1, 1]), True)
    same("parity fills 1 ? 1 1 0 0 1 0", parity_fill([1, None, 1, 1, 0, 0, 1, 0]), [1, 0, 1, 1, 0, 0, 1, 0])
    same("parity cannot fill two erasures", parity_fill([1, None, 1, None, 0, 0, 1, 0]), None)
    same("soft R3 decision on 0.9, -1.4, 0.2", soft_repetition([0.9, -1.4, 0.2]), 1)
    same("hard R3 decision on 0.9, -1.4, 0.2", hard_repetition([0.9, -1.4, 0.2]), 0)
    check("sum of 0.9, -1.4, 0.2", 0.9 - 1.4 + 0.2, "-0.3")
    check("Eb/N0 = 4 dB as a ratio", 10 ** 0.4, "2.5119")
    check("noise sigma for uncoded BPSK at 4 dB", awgn_sigma(4.0, 1.0), "0.4462")
    check("noise sigma for R3 pulses at 4 dB", awgn_sigma(4.0, 1 / 3), "0.7728")


def check_box_muller() -> None:
    print("Box-Muller")
    u1, u2 = 0.25, 0.1
    r = sqrt(-2 * log(u1))
    check("worked pair: radius for u1 = 0.25", r, "1.6651")
    check("worked pair: angle for u2 = 0.1, radians", 2 * pi * u2, "0.6283")
    check("worked pair: z0", r * cos(2 * pi * u2), "1.3471")
    check("worked pair: z1", r * sin(2 * pi * u2), "0.9787")
    check("worked pair: cos(2 pi 0.1)", cos(2 * pi * u2), "0.8090")
    check("worked pair: sin(2 pi 0.1)", sin(2 * pi * u2), "0.5878")
    check("P(R^2 > 2) = e^(-1)", exp(-1), "0.3679")
    rng = random.Random(BOX_MULLER_SEED)
    a, b = rng.random(), rng.random()
    same("first uniform from seed 17, as Python prints it", repr(a), "0.5219839097124932")
    same("second uniform from seed 17, as Python prints it", repr(b), "0.8066907771186791")
    check("first uniform from seed 17", a, "0.52198")
    check("second uniform from seed 17", b, "0.80669")
    check("u1 = 1 - first uniform", 1 - a, "0.47802")
    check("R^2 from seed 17", -2 * log(1 - a), "1.4762")
    r = sqrt(-2 * log(1 - a))
    check("radius from seed 17", r, "1.2150")
    check("angle from seed 17, radians", 2 * pi * b, "5.0686")
    check("cos of the seed-17 angle", cos(2 * pi * b), "0.3487")
    check("sin of the seed-17 angle", sin(2 * pi * b), "-0.9372")
    zs = box_muller_samples(BOX_MULLER_COUNT, BOX_MULLER_SEED)
    check("first sample from seed 17", zs[0], "0.4237")
    check("second sample from seed 17", zs[1], "-1.1387")
    mean = sum(zs) / len(zs)
    var = sum((z - mean) ** 2 for z in zs) / len(zs)
    tail = sum(abs(z) > Z95 for z in zs)
    check("mean of 100,000 samples, seed 17", mean, "0.0044")
    check("that mean in standard errors", mean * sqrt(len(zs)), "1.4")
    check("variance of 100,000 samples, seed 17", var, "0.9981")
    same("samples beyond 1.96 in size, of 100,000", tail, 5014)
    check("2Q(1.96)", 2 * q(Z95), "0.0500")
    ph = tail / len(zs)
    check("tail fraction", ph, "0.05014")
    check("one minus the tail fraction", 1 - ph, "0.94986")
    se = sqrt(ph * (1 - ph) / len(zs))
    check("estimated standard deviation of the tail fraction", se, "0.00069012")
    check("Wald half-width for the tail fraction", Z95 * se, "0.0013526")
    lo, hi = wald(tail, len(zs))
    check("Wald low for the tail count", lo, "0.04879")
    check("Wald high for the tail count", hi, "0.05149")
    lo, hi = wilson(tail, len(zs))
    check("Wilson low for the tail count", lo, "0.04880")
    check("Wilson high for the tail count", hi, "0.05151")
    counts = [0] * 32
    for z in zs:
        if -4 <= z < 4:
            counts[int((z + 4) // 0.25)] += 1
    peak = max(range(32), key=lambda i: counts[i])
    same("lower edge of the tallest histogram bar", -4 + 0.25 * peak, -0.25)
    same("samples in the tallest histogram bar", counts[peak], 9907)
    same("samples outside -4 to 4", len(zs) - sum(counts), 3)
    check("peak of the standard normal density", 1 / sqrt(2 * pi), "0.3989")


def check_estimator() -> None:
    print("Estimator and intervals")
    name, seed, trials, k, lo, hi, exact = run_experiments()[0]
    same("R3 on BSC(0.1), seed 2026, 10,000 words: failures", k, 266)
    from ch09_error_control_strategies_and_performance import simulate_repetition
    same("chapter 9's simulator on the same seed", sum(simulate_repetition(3, 0.1, 10_000, random.Random(2026))), 266)
    check("exact R3 word error at p = 0.1", exact, "0.028")
    check("Wilson low", lo, "0.02362")
    check("Wilson high", hi, "0.02994")
    wl, wh = wald(k, trials)
    check("Wald low", wl, "0.02345")
    check("Wald high", wh, "0.02975")

    print("Coverage of the intervals, 1,000 seeded runs")
    longs = coverage_counts(COVERAGE_LONG)
    shorts = coverage_counts(COVERAGE_SHORT)
    estimates = [k / COVERAGE_LONG for k in longs]
    mean = sum(estimates) / len(estimates)
    sd = sqrt(sum((e - mean) ** 2 for e in estimates) / (len(estimates) - 1))
    check("mean of 1,000 estimates at N = 1,000", mean, "0.0280")
    check("their standard deviation", sd, "0.00520")
    check("predicted standard deviation sqrt(p(1-p)/N)",
          sqrt(COVERAGE_EXACT * (1 - COVERAGE_EXACT) / COVERAGE_LONG), "0.00522")
    same("fewest failures in a run of 1,000 words", min(longs), 13)
    same("most failures in a run of 1,000 words", max(longs), 47)

    def inside(interval, k, n):
        a, b = interval(k, n)
        return a <= COVERAGE_EXACT <= b

    same("Wilson intervals covering 0.028, N = 1,000", sum(inside(wilson, k, COVERAGE_LONG) for k in longs), 957)
    same("Wald intervals covering 0.028, N = 1,000", sum(inside(wald, k, COVERAGE_LONG) for k in longs), 947)
    same("Wilson intervals covering 0.028, N = 100", sum(inside(wilson, k, COVERAGE_SHORT) for k in shorts), 981)
    same("Wald intervals covering 0.028, N = 100", sum(inside(wald, k, COVERAGE_SHORT) for k in shorts), 950)
    same("runs at N = 100 with no failure", sum(k == 0 for k in shorts), 49)
    check("probability of no failure in 100 words", (1 - COVERAGE_EXACT) ** 100, "0.0584")
    same("Wald misses at N = 100 that come from K = 0",
         sum(k == 0 and not inside(wald, k, COVERAGE_SHORT) for k in shorts), 49)
    same("Wilson misses at N = 100 with K below 7",
         sum(k < 7 and not inside(wilson, k, COVERAGE_SHORT) for k in shorts), 0)
    check("Wald low at K = 1, N = 100", wald(1, 100)[0], "-0.0095")
    same("misses among the first 50 Wilson intervals at N = 1,000",
         sum(not inside(wilson, k, COVERAGE_LONG) for k in longs[:50]), 2)
    same("those two misses as (seed, failures)",
         [(i + 1, k) for i, k in enumerate(longs[:50]) if not inside(wilson, k, COVERAGE_LONG)], [(20, 45), (50, 16)])
    same("failure counts of the Wilson misses at N = 100",
         sorted({k for k in shorts if not inside(wilson, k, COVERAGE_SHORT)}), [7, 8, 9])
    same("Wald misses at N = 100", sum(not inside(wald, k, COVERAGE_SHORT) for k in shorts), 50)
    check("exact Wilson coverage at P = 0.028, N = 100", exact_coverage(COVERAGE_EXACT, COVERAGE_SHORT), "0.9774")
    check("Wilson low at K = 1, N = 100", wilson(1, 100)[0], "0.0018")
    holds("at P = 0.0017, N = 100 only K = 0 covers",
          wilson(0, 100)[0] <= 0.0017 <= wilson(0, 100)[1] and wilson(1, 100)[0] > 0.0017)
    check("1 - 0.0017", 1 - 0.0017, "0.9983")
    check("exact Wilson coverage at P = 0.0017, N = 100", exact_coverage(0.0017, COVERAGE_SHORT), "0.8435")

    print("Zero failures and the rule of three")
    k0 = count_failures(ZERO_TRIALS, ZERO_SEED, repetition_bsc_trial(5, 0.01))
    same("R5 on BSC(0.01), seed 18, 20,000 words: failures", k0, 0)
    p5 = repetition_bsc_exact(5, 0.01)
    check("exact R5 word error at p = 0.01", p5, "9.8506e-6")
    check("expected failures in 20,000 words", ZERO_TRIALS * p5, "0.197")
    check("probability of no failure in 20,000 words", (1 - p5) ** ZERO_TRIALS, "0.821")
    check("-ln 0.05", -log(0.05), "2.9957")
    check("exact zero-failure bound for N = 20,000", zero_failure_bound(ZERO_TRIALS), "1.4978e-4")
    check("rule of three for N = 20,000", 3 / ZERO_TRIALS, "1.5e-4")
    check("rule-of-three bound over the true R5 rate", 3 / ZERO_TRIALS / p5, "15")
    check("e^-3", exp(-3), "0.0498")
    check("Wilson high at K = 0, N = 20,000", wilson(0, ZERO_TRIALS)[1], "1.9204e-4")
    check("-ln 0.025", -log(0.025), "3.6889")
    same("Wald interval at K = 0 has zero width", wald(0, ZERO_TRIALS), (0.0, 0.0))

    print("Trials for a target precision")
    check("(1.96/0.2)^2", (Z95 / 0.2) ** 2, "96.04")
    check("(1.96/0.1)^2", (Z95 / 0.1) ** 2, "384.16")
    pu = parity_undetected_exact(8, 0.01)
    same("trials for 10% precision on the parity escape rate", trials_needed(pu, 0.1), 145_307)
    same("trials for 20% precision on the parity escape rate", trials_needed(pu, 0.2), 36_327)
    check("trials for 20% at P = 1e-5, millions", trials_needed(1e-5, 0.2) / 1e6, "9.60")


# What the chapter's summary table states for each experiment: failures, Wilson
# limits, exact value, and whether the exact value lies inside the interval.
STATED = {
    "R3 on BSC(0.1)": (266, "0.02362", "0.02994", "0.028", True),
    "R3 on BEC(0.2)": (741, "0.006897", "0.007961", "0.008", False),
    "parity (8,7) on BSC(0.01), undetected": (265, "0.002350", "0.002988", "0.0026368", True),
    "parity (8,7) on BEC(0.1), unfilled": (3720, "0.18067", "0.19145", "0.18690", True),
    "uncoded BPSK at 4 dB": (1246, "0.011791", "0.013167", "0.012501", True),
    "hard R3 at 4 dB": (2691, "0.025925", "0.027931", "0.026835", True),
    "soft R3 at 4 dB": (1250, "0.011830", "0.013208", "0.012501", True),
}


def check_closed_forms() -> None:
    print("Closed-form checks")
    rows = run_experiments()
    for name, seed, trials, k, lo, hi, exact in rows:
        sk, slo, shi, sexact, sinside = STATED[name]
        same(f"{name}: failures in {trials:,}, seed {seed}", k, sk)
        check(f"{name}: Wilson low", lo, slo)
        check(f"{name}: Wilson high", hi, shi)
        check(f"{name}: exact", exact, sexact)
        same(f"{name}: exact value inside the interval", lo <= exact <= hi, sinside)
    same("intervals in the table that miss, of 7", sum(not (r[4] <= r[6] <= r[5]) for r in rows), 1)

    print("The erasure miss")
    sd = sqrt(100_000 * 0.008 * 0.992)
    check("standard deviation of K for R3 on BEC(0.2), N = 100,000", sd, "28.17")
    check("distance of 741 from 800 in standard deviations", (741 - 800) / sd, "-2.09")
    k = count_failures(1_000_000, 12, repetition_bec_trial(3, 0.2))
    same("same run continued to 1,000,000 words: failures", k, 7974)
    lo, hi = wilson(k, 1_000_000)
    check("its Wilson low", lo, "0.007802")
    check("its Wilson high", hi, "0.008150")
    trial = repetition_bec_trial(3, 0.2)
    covered = 0
    for s in range(1, 101):
        lo, hi = wilson(count_failures(100_000, s, trial), 100_000)
        covered += lo <= 0.008 <= hi
    same("runs with seeds 1 to 100 whose interval covers 0.008", covered, 98)

    print("Parity and Gaussian details")
    check("0.98^8", 0.98**8, "0.850763")
    check("0.99^8", 0.99**8, "0.922745")
    check("0.9^8", 0.9**8, "0.430467")
    check("8 x 0.1 x 0.9^7", 8 * 0.1 * 0.9**7, "0.382638")
    check("parity detected rate, n = 8, p = 0.01", parity_detected_exact(8, 0.01), "0.07462")
    same("parity words flagged, seed 13, 100,000 words", count_failures(100_000, 13, parity_detected_trial(8, 0.01)),
         7565)
    lo, hi = wilson(265, 100_000)
    check("relative half-width of the parity escape interval", (hi - lo) / 2 / (265 / 100_000), "0.120")
    check("1.96 / sqrt(265)", Z95 / sqrt(265), "0.120")
    check("argument of Q for uncoded BPSK at 4 dB", sqrt(2 * 10**0.4), "2.2414")
    check("argument of Q for each R3 pulse at 4 dB", sqrt(2 * 10**0.4 / 3), "1.2941")
    check("symbol error of each R3 pulse at 4 dB", q(sqrt(2 * 10**0.4 / 3)), "0.09782")
    both, hard_only, soft_only = paired_hard_soft(100_000, 16, 4.0)
    same("words both decoders get wrong, seed 16", both, 894)
    same("words only the hard decoder gets wrong", hard_only, 1797)
    same("words only the soft decoder gets wrong", soft_only, 356)
    same("extra failures of hard over soft", hard_only - soft_only, 1441)

    print("Gaussian sweep, seeds 100 + dB, 100,000 trials per point")
    misses = []
    sweep = awgn_sweep()
    for db, ku, kh, ks in sweep:
        for label, k, exact in (("uncoded", ku, uncoded_bpsk_exact(db)),
                                ("hard R3", kh, hard_repetition_awgn_exact(3, db)),
                                ("soft R3", ks, soft_repetition_awgn_exact(3, db))):
            lo, hi = wilson(k, SWEEP_TRIALS)
            if not lo <= exact <= hi:
                misses.append(f"{label} at {db} dB")
    same("sweep points whose interval misses the closed form, of 24", len(misses), 1)
    same("the one miss", misses, ["hard R3 at 1 dB"])
    same("hard R3 failures at 1 dB, seed 101", sweep[1][2], 8727)
    check("expected hard R3 failures at 1 dB", SWEEP_TRIALS * hard_repetition_awgn_exact(3, 1), "8,536.0")
    check("misses expected among 24 intervals", 24 * 0.05, "1.2")
    check("chance that 24 independent 95% intervals all cover", 0.95**24, "0.2920")
    same("uncoded, soft and hard failures at 7 dB, seed 107", (sweep[7][1], sweep[7][3], sweep[7][2]), (92, 63, 308))


def check_worked_numbers() -> None:
    print("Other numbers the chapter works by hand")
    check("-2 ln 0.25", -2 * log(0.25), "2.7726")
    check("standard error of a mean of 100,000 standard normal samples", 1 / sqrt(100_000), "0.0032")
    noise = box_muller_samples(2, BOX_MULLER_SEED)
    sigma = awgn_sigma(4.0, 1.0)
    check("received value for +1 plus the first seed-17 sample at 4 dB", 1 + sigma * noise[0], "1.1890")
    check("received value for +1 plus the second seed-17 sample at 4 dB", 1 + sigma * noise[1], "0.4919")
    check("(1 + 0.98^8) / 2", (1 + 0.98**8) / 2, "0.92538151")
    check("0.99^8 to eight places", 0.99**8, "0.92274469")
    check("their difference", (1 + 0.98**8) / 2 - 0.99**8, "0.00263682")
    check("0.95^7", 0.95**7, "0.6983")
    check("chance that at least one of 7 intervals misses", 1 - 0.95**7, "0.3017")
    check("misses expected among 7 intervals", 7 * 0.05, "0.35")
    check("1.96^2", Z95**2, "3.8416")
    check("-ln(0.025) / 20,000", -log(0.025) / ZERO_TRIALS, "1.8444e-4")
    check("trials for 100 failures at the R5 rate, millions", 100 / repetition_bsc_exact(5, 0.01) / 1e6, "10.2")
    holds("2^19937 - 1 has more than 6,000 digits", 2**19937 - 1 >= 10**6000)
    check("Box-Muller angle for u2 = 0.1, degrees", degrees(2 * pi * 0.1), "36")
    check("misses expected among 50 intervals", 50 * 0.05, "2.5")
    lo, hi = wilson(45, 1000)
    check("Wilson low for run 20, 45 failures in 1,000", lo, "0.0338")
    check("Wilson high for run 20", hi, "0.0597")
    lo, hi = wilson(16, 1000)
    check("Wilson low for run 50, 16 failures in 1,000", lo, "0.0099")
    check("Wilson high for run 50", hi, "0.0258")
    check("failures expected in 100 words", 100 * repetition_bsc_exact(3, 0.1), "2.8")
    check("runs of 1,000 expected to see no failure in 100 words", 1000 * (1 - repetition_bsc_exact(3, 0.1)) ** 100, "58")
    check("Wilson high at K = 0, N = 100", wilson(0, 100)[1], "0.0370")
    check("Wilson low at K = 7, N = 100", wilson(7, 100)[0], "0.0343")
    holds("Wilson at K = 6, N = 100 still reaches below 0.028", wilson(6, 100)[0] < 0.028)
    holds("Wald at N = 100 misses 0.028 from above exactly when K >= 9",
          wald(8, 100)[0] < 0.028 < wald(9, 100)[0])
    p = 0.01
    check("R5 term 10 p^3 (1-p)^2 at p = 0.01", 10 * p**3 * (1 - p) ** 2, "9.801e-6")
    check("R5 term 5 p^4 (1-p) at p = 0.01", 5 * p**4 * (1 - p), "4.95e-8")
    check("expected parity flags in 100,000 words", 100_000 * parity_detected_exact(8, 0.01), "7,461.8")
    pu = parity_undetected_exact(8, 0.01)
    check("(1 - P)/P for the parity escape rate", (1 - pu) / pu, "378.245")
    check("expected uncoded failures in 100,000 words at 4 dB", 100_000 * uncoded_bpsk_exact(4), "1,250.1")
    check("expected hard R3 failures in 100,000 words at 4 dB", 100_000 * hard_repetition_awgn_exact(3, 4), "2,683.5")
    check("expected soft R3 failures in 100,000 words at 4 dB", 100_000 * soft_repetition_awgn_exact(3, 4), "1,250.1")
    check("hard over soft at 4 dB", hard_repetition_awgn_exact(3, 4) / soft_repetition_awgn_exact(3, 4), "2.147")
    check("N + z^2 at N = 20,000", ZERO_TRIALS + Z95**2, "20,003.8416")
    print("Exercise answers")
    check("-2 ln 0.5", -2 * log(0.5), "1.3863")
    check("radius for u1 = 0.5", sqrt(-2 * log(0.5)), "1.1774")
    check("(1 - P)/P at P = 1e-4", (1 - 1e-4) / 1e-4, "9,999")
    check("384.16 x 9,999", (Z95 / 0.1) ** 2 * (1 - 1e-4) / 1e-4, "3,841,215.8")
    same("trials for 10% precision at P = 1e-4", trials_needed(1e-4, 0.1), 3_841_216)
    check("zero-failure bound for N = 1,000,000", zero_failure_bound(1_000_000), "2.9957e-6")
    check("(1 + 0.8^4) / 2", (1 + 0.8**4) / 2, "0.7048")
    check("0.9^4", 0.9**4, "0.6561")
    check("two flips in four at p = 0.1", comb(4, 2) * 0.1**2 * 0.9**2, "0.0486")
    check("four flips in four at p = 0.1", 0.1**4, "0.0001")
    check("parity undetected rate, n = 4, p = 0.1", parity_undetected_exact(4, 0.1), "0.0487")
    same("hard R3 decision on 0.3, 0.4, -1.2", hard_repetition([0.3, 0.4, -1.2]), 0)
    same("soft R3 decision on 0.3, 0.4, -1.2", soft_repetition([0.3, 0.4, -1.2]), 1)
    same("BSC(0.1) on 1011 with draws 0.31, 0.04, 0.77, 0.09",
         [b ^ int(u < 0.1) for b, u in zip([1, 0, 1, 1], [0.31, 0.04, 0.77, 0.09])], [1, 1, 1, 0])


def main() -> None:
    check_channels_and_codes()
    check_box_muller()
    check_estimator()
    check_closed_forms()
    check_worked_numbers()
    print("\nAll chapter 11 numbers reproduced.")


if __name__ == "__main__":
    main()
