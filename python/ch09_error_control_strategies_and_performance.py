"""Reproduce every number in chapter 9 of Error Correction Coding, "Error-Control Strategies and Performance".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch09_error_control_strategies_and_performance.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. The ARQ throughput
formulas are checked against a frame-by-frame simulation of each protocol, the
binomial word-error formula against exhaustive enumeration of every error
pattern, and every Monte Carlo estimate, which uses random.Random with a fixed
seed so that it reproduces exactly, against its exact formula. A run that
reaches the final line has reproduced every number in the chapter.

The module also holds the functions that ch09_charts.py plots, so the charts
and this check share one source.
"""

import heapq
import random
import sys
from itertools import product
from math import ceil, comb, erfc, exp, floor, lgamma, log, log10, sqrt


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
    print(f"  {label:<62} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<62} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def near(label: str, a: float, b: float, rel: float) -> None:
    """Two computations of one quantity that must agree to a relative tolerance."""
    print(f"  {label:<62} {a:<16.9g} {b:.9g}")
    if abs(a - b) > rel * abs(b):
        print(f"MISMATCH in '{label}': {a!r} against {b!r}")
        sys.exit(1)


def section(title: str) -> None:
    print(f"\n{title}")


# The Gaussian tail ----------------------------------------------------------------

def q_func(x: float) -> float:
    """The Q-function, the probability that a standard Gaussian exceeds x."""
    return 0.5 * erfc(x / sqrt(2.0))


def q_inv(p: float) -> float:
    """The x with Q(x) = p, found by bisection."""
    lo, hi = -40.0, 40.0
    for _ in range(300):
        mid = (lo + hi) / 2
        if q_func(mid) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def db(x: float) -> float:
    return 10 * log10(x)


# Binomial counting -----------------------------------------------------------------

def binom_pmf(n: int, i: int, p: float) -> float:
    return comb(n, i) * p ** i * (1 - p) ** (n - i)


def word_error(n: int, t: int, p: float) -> float:
    """Probability that more than t of n bits flip on BSC(p), the binomial word-error formula."""
    return sum(binom_pmf(n, i, p) for i in range(t + 1, n + 1))


def binom_log_pmf(n: int, k: int, p: float) -> float:
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + k * log(p) + (n - k) * log(1 - p)


def window_prob(n: int, p: float, rel: float) -> float:
    """Exact probability that K/n lands within rel*p of p when K is binomial(n, p)."""
    lo = ceil(n * p * (1 - rel) - 1e-9)
    hi = floor(n * p * (1 + rel) + 1e-9)
    return sum(exp(binom_log_pmf(n, k, p)) for k in range(max(lo, 0), min(hi, n) + 1))


# ARQ throughput ------------------------------------------------------------------

def frame_error(n: int, p: float) -> float:
    return 1.0 - (1.0 - p) ** n


def eta_sw(P: float, a: float, rate: float) -> float:
    return (1 - P) / (1 + a) * rate


def eta_gbn(P: float, window: float, rate: float) -> float:
    return (1 - P) / (1 + (window - 1) * P) * rate


def eta_sr(P: float, rate: float) -> float:
    return (1 - P) * rate


def simulate_protocol(kind: str, P: float, a: int, frames: int, rng: random.Random) -> float:
    """Slots, in frame times, that a protocol spends delivering a number of frames.

    Each transmission of a frame fails independently with probability P. A
    stop-and-wait attempt holds the channel for 1 + a slots. Go-back-N with a
    window of 1 + a slots loses the failed frame and the a frames sent behind
    it. Selective repeat loses only the failed frame. Because each failure is
    charged the slots the chapter's proof charges, this checks the averaging in
    that proof, not the slot model of the protocols.
    """
    slots = 0
    for _ in range(frames):
        while True:
            failed = rng.random() < P
            if kind == "sw":
                slots += 1 + a
            elif kind == "gbn":
                slots += (1 + a) if failed else 1
            else:
                slots += 1
            if not failed:
                break
    return slots / frames


def sr_release_delays(a: int, P: float, slots: int, seed: int) -> list:
    """Extra delay, in slots, of each frame under ideal selective repeat with in-order release.

    The sender fills every slot, resending a failed frame as soon as its NAK
    arrives a slots after the frame ends and otherwise sending the next new
    frame. The receiver buffers frames that arrive out of order and releases a
    frame only once every earlier frame is in, so a frame's delay counts the
    wait behind earlier damaged frames as well as its own retries. Delay runs
    from the frame's first transmission to its release.
    """
    rng = random.Random(seed)
    due, first, recv, nxt = [], {}, {}, 0
    for t in range(1, slots + 1):
        if due and due[0][0] <= t - 1:
            _, f = heapq.heappop(due)
        else:
            f, nxt = nxt, nxt + 1
            first[f] = t
        if rng.random() < P:
            heapq.heappush(due, (t + a, f))
        else:
            recv[f] = t
    released, delays = -1, []
    for f in range(nxt):
        if f not in recv:
            break
        released = max(released, recv[f])
        delays.append(released - first[f])
    return delays


def gbn_delivery_delays(a: int, P: float, slots: int, seed: int) -> list:
    """Extra delay, in slots, of each frame under go-back-N with a window of 1 + a.

    The receiver accepts only the next frame in order. When that frame fails,
    the sender learns of it a slots after the frame ends and goes back to it,
    resending it and every frame sent behind it. Delay runs from the frame's
    first transmission to its acceptance.
    """
    rng = random.Random(seed)
    first, recv, nxt, expect, goback = {}, {}, 0, 0, None
    for t in range(1, slots + 1):
        if goback and goback[0] <= t - 1:
            nxt, goback = goback[1], None
        f, nxt = nxt, nxt + 1
        first.setdefault(f, t)
        ok = rng.random() >= P
        if f == expect and ok:
            recv[f] = t
            expect += 1
        elif f == expect and goback is None:
            goback = (t + a, f)
    return [recv[f] - first[f] for f in sorted(recv)]


def hybrid_failure(n: int, t: int, p: float) -> float:
    return word_error(n, t, p)


def hybrid_crossover(n: int, rate_crc: float, rate_hybrid: float, t: int) -> float:
    """The bit error probability where type-I hybrid ARQ and plain selective repeat tie."""
    lo, hi = 1e-7, 1e-3
    for _ in range(200):
        mid = sqrt(lo * hi)
        if eta_sr(frame_error(n, mid), rate_crc) > (1 - hybrid_failure(n, t, mid)) * rate_hybrid:
            lo = mid
        else:
            hi = mid
    return sqrt(lo * hi)


# Single parity under ARQ -----------------------------------------------------------

def parity_outcomes(n: int, p: float) -> tuple:
    """Probabilities that one parity-checked transmission is clean, flagged, or wrong but accepted."""
    clean = (1 - p) ** n
    alarm = (1 - (1 - 2 * p) ** n) / 2
    slip = 1 - clean - alarm
    return clean, alarm, slip


def parity_outcomes_enumerated(n: int, p: float) -> tuple:
    clean = alarm = slip = 0.0
    for e in product((0, 1), repeat=n):
        w = sum(e)
        pr = p ** w * (1 - p) ** (n - w)
        if w == 0:
            clean += pr
        elif w % 2:
            alarm += pr
        else:
            slip += pr
    return clean, alarm, slip


def simulate_arq_parity(n: int, p: float, frames: int, rng: random.Random) -> tuple:
    """Deliver frames under ARQ with single parity, counting wrong deliveries and transmissions."""
    wrong = sends = 0
    for _ in range(frames):
        while True:
            sends += 1
            w = sum(rng.random() < p for _ in range(n))
            if w % 2 == 0:
                wrong += w > 0
                break
    return wrong / frames, sends / frames


# Coding gain ---------------------------------------------------------------------

def required_uncoded(target: float) -> tuple:
    """The Gaussian-channel argument x and the Eb/N0 in dB at which uncoded BPSK meets a target."""
    x = q_inv(target)
    return x, db(x * x / 2)


def required_repetition(n: int, target: float) -> tuple:
    """Channel flip probability, Q argument, Es/N0 and Eb/N0 in dB for hard-decision R_n to meet a target."""
    t = (n - 1) // 2
    lo, hi = 1e-15, 0.5
    for _ in range(300):
        mid = (lo + hi) / 2
        if word_error(n, t, mid) > target:
            hi = mid
        else:
            lo = mid
    p = (lo + hi) / 2
    x = q_inv(p)
    es = x * x / 2
    return p, x, db(es), db(es * n)


def hard_rep_ber(n: int, ebn0_db: float) -> float:
    """Decoded bit error rate of hard-decision R_n on the Gaussian channel."""
    ebn0 = 10 ** (ebn0_db / 10)
    p = q_func(sqrt(2 * ebn0 / n))
    return word_error(n, (n - 1) // 2, p)


def uncoded_ber(ebn0_db: float) -> float:
    return q_func(sqrt(2 * 10 ** (ebn0_db / 10)))


def asymptotic_hard_gain(n: int) -> float:
    t = (n - 1) // 2
    return db((t + 1) / n)


def soft_rep_ber_simulated(n: int, ebn0_db: float, trials: int, rng: random.Random) -> float:
    """Soft-decision R_n: add the n received values and take the sign."""
    ebn0 = 10 ** (ebn0_db / 10)
    es = ebn0 / n
    sigma = sqrt(0.5)
    errors = 0
    for _ in range(trials):
        total = sum(sqrt(es) + rng.gauss(0.0, sigma) for _ in range(n))
        errors += total < 0
    return errors / trials


# Monte Carlo ---------------------------------------------------------------------

Z95 = 1.96


def simulate_repetition(n: int, p: float, trials: int, rng: random.Random) -> list:
    """One entry per trial, 1 when majority voting over n copies on BSC(p) fails."""
    t = (n - 1) // 2
    out = []
    for _ in range(trials):
        flips = sum(rng.random() < p for _ in range(n))
        out.append(1 if flips > t else 0)
    return out


def wald(k: int, n: int, z: float = Z95) -> tuple:
    ph = k / n
    hw = z * sqrt(ph * (1 - ph) / n)
    return ph - hw, ph + hw


def wilson(k: int, n: int, z: float = Z95) -> tuple:
    ph = k / n
    centre = (ph + z * z / (2 * n)) / (1 + z * z / n)
    hw = z / (1 + z * z / n) * sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return centre - hw, centre + hw


def trials_needed(p: float, rel: float, z: float = Z95) -> float:
    return (z / rel) ** 2 * (1 - p) / p


MAIN_SEED = 2026
MAIN_TRIALS = 10_000
MAIN_P = 0.1


def main_run() -> list:
    return simulate_repetition(3, MAIN_P, MAIN_TRIALS, random.Random(MAIN_SEED))


POINT_PS = (0.02, 0.05, 0.1, 0.2, 0.3)
POINT_SEED = 7
POINT_TRIALS = 10_000


def mc_points() -> list:
    """Estimates of the R_3 word error rate at several p, with exact values and Wilson intervals."""
    rng = random.Random(POINT_SEED)
    rows = []
    for p in POINT_PS:
        k = sum(simulate_repetition(3, p, POINT_TRIALS, rng))
        lo, hi = wilson(k, POINT_TRIALS)
        rows.append((p, word_error(3, 1, p), k, POINT_TRIALS, lo, hi))
    return rows


# Reliability -----------------------------------------------------------------------

HOURS_PER_YEAR = 8760


def mtbf_hours(fit: float) -> float:
    return 1e9 / fit


# The checks ------------------------------------------------------------------------

def check_arq() -> None:
    section("Stop-and-wait, go-back-N and selective repeat")
    n, k, R, t_rt, p = 1000, 968, 1e6, 10e-3, 1e-4
    t_f = n / R
    a = t_rt / t_f
    window = 1 + a
    rate = k / n
    q = (1 - p) ** n
    P = frame_error(n, p)
    check("frame time in ms", t_f * 1e3, "1")
    check("a = T_rt / T_f", a, "10")
    check("window N = 1 + a", window, "11")
    check("frame efficiency k/n", rate, "0.968")
    check("clean frame probability (chapter 1)", q, "0.904833")
    check("frame error probability P", P, "0.095167")
    check("stop-and-wait slots per delivered frame", (1 + a) / q, "12.1569")
    check("go-back-N slots per delivered frame", (1 + (window - 1) * P) / q, "2.15694")
    check("selective repeat slots per delivered frame", 1 / q, "1.10518")
    check("go-back-N denominator 1 + (N - 1)P", 1 + (window - 1) * P, "1.95167")
    sw, gbn, sr = eta_sw(P, a, rate), eta_gbn(P, window, rate), eta_sr(P, rate)
    check("stop-and-wait efficiency", sw, "0.079625")
    check("go-back-N efficiency", gbn, "0.44878")
    check("selective repeat efficiency", sr, "0.87588")
    check("stop-and-wait goodput kbit/s", sw * R / 1e3, "79.6")
    check("go-back-N goodput kbit/s", gbn * R / 1e3, "448.8")
    check("selective repeat goodput kbit/s", sr * R / 1e3, "875.9")
    check("go-back-N over stop-and-wait", gbn / sw, "5.64")
    check("selective repeat over go-back-N", sr / gbn, "1.95")
    check("stop-and-wait ceiling k/(n(1 + a)) with P = 0", rate / (1 + a), "0.088")
    # Toy timeline of figure 9-2: a = 3, 16 slots, frame 2 damaged once.
    same("toy timeline, stop-and-wait frames delivered", 16 // (1 + 3) - 1, 3)
    same("toy timeline, go-back-N frames delivered (window 4 lost)", 16 - (1 + 3), 12)
    same("toy timeline, ideal selective repeat frames delivered", 16 - 1, 15)
    # Held to a window of 4, selective repeat sends frames 3 to 5, resends 2 in
    # slot 6, idles in slots 7 to 9 until that resend is acknowledged, then sends
    # frames 6 to 12 in slots 10 to 16.
    same("toy timeline, selective repeat held to a window of 4", 5 + (16 - 9), 12)
    rng = random.Random(11)
    for kind, formula in (("sw", (1 + a) / q), ("gbn", (1 + (window - 1) * P) / q), ("sr", 1 / q)):
        near(f"{kind} slots per frame, simulated over 200,000 frames", simulate_protocol(kind, P, int(a), 200_000, rng),
             formula, 0.01)
    check("mean retries P/(1 - P)", P / q, "0.10518")
    check("one retry costs (1 + a) frame times, ms", (1 + a) * t_f * 1e3, "11")
    check("mean extra delay in ms", P / q * (1 + a) * t_f * 1e3, "1.16")
    check("probability of two or more retries", P ** 2, "0.0090568")
    check("one frame in this many waits two retries", 1 / P ** 2, "110")
    check("extra wait of two retries, ms", 2 * (1 + a) * t_f * 1e3, "22")
    # In-order release: a frame also waits behind earlier damaged frames.
    sr_d = sr_release_delays(int(a), P, 400_000, 7)
    check("selective repeat, simulated mean extra delay, ms", sum(sr_d) / len(sr_d) * t_f * 1e3, "5.68")
    check("selective repeat, one frame in this many waits 22 ms or more",
          len(sr_d) / sum(d >= 22 for d in sr_d), "58")
    gbn_d = gbn_delivery_delays(int(a), P, 400_000, 7)
    near("go-back-N simulated throughput against 1/S", len(gbn_d) / 400_000, q / (1 + (window - 1) * P), 0.01)
    # Each resend of a frame trails its previous attempt by exactly 1 + a slots, and a
    # frame is sent S_GBN times on average, so the mean extra delay is (1 + a)(S_GBN - 1) T_f.
    s_gbn = (1 + (window - 1) * P) / q
    gbn_mean_ms = (1 + a) * (s_gbn - 1) * t_f * 1e3
    check("go-back-N, mean extra delay (1 + a)(S_GBN - 1) T_f, ms", gbn_mean_ms, "12.7")
    near("go-back-N simulated mean extra delay against that", sum(gbn_d) / len(gbn_d) * t_f * 1e3, gbn_mean_ms, 0.02)
    check("go-back-N, one frame in this many waits 22 ms or more",
          len(gbn_d) / sum(d >= 22 for d in gbn_d), "3")


def check_satellite() -> None:
    section("Worked example: a satellite hop")
    n, k, R, t_rt, p = 1000, 968, 1e6, 0.48, 1e-5
    a = t_rt / (n / R)
    window = 1 + a
    P = frame_error(n, p)
    rate = k / n
    check("one-way path, km", 2 * 35_786, "71,572")
    check("one-way time, s", 2 * 35_786 / 299_792.458, "0.239")
    check("round trip, s", 4 * 35_786 / 299_792.458, "0.477")
    check("a", a, "480")
    check("window", window, "481")
    check("frame error probability", P, "0.0099502")
    check("stop-and-wait efficiency", eta_sw(P, a, rate), "0.0019924")
    check("go-back-N efficiency", eta_gbn(P, window, rate), "0.16592")
    check("go-back-N denominator", 1 + (window - 1) * P, "5.7761")
    check("selective repeat efficiency", eta_sr(P, rate), "0.95837")
    check("frames held after one failure, about 1 + a, in kilobytes", window * n / 8 / 1000, "60.1")


def check_residual() -> None:
    section("The residual error of a checked link")
    for n, p, stated in ((8, 0.01, ("0.850763", "0.922745", "0.074618", "0.0026368", "1.080635", "0.0028494")),
                         (1000, 1e-4, ("0.818714", "0.904833", "0.090643", "0.0045243", "1.099678", "0.0049753"))):
        clean, alarm, slip = parity_outcomes(n, p)
        check(f"n = {n}: (1 - 2p)^n", (1 - 2 * p) ** n, stated[0])
        check(f"n = {n}: clean", clean, stated[1])
        check(f"n = {n}: alarm", alarm, stated[2])
        check(f"n = {n}: slip", slip, stated[3])
        check(f"n = {n}: expected transmissions", 1 / (1 - alarm), stated[4])
        check(f"n = {n}: residual per accepted frame", slip / (1 - alarm), stated[5])
    enumerated = parity_outcomes_enumerated(8, 0.01)
    for name, a, b in zip(("clean", "alarm", "slip"), enumerated, parity_outcomes(8, 0.01)):
        near(f"n = 8 {name}, all 256 patterns enumerated", a, b, 1e-9)
    wrong, sends = simulate_arq_parity(8, 0.01, 400_000, random.Random(5))
    near("n = 8 residual, simulated over 400,000 frames", wrong, parity_outcomes(8, 0.01)[2] / (1 - parity_outcomes(8, 0.01)[1]), 0.05)
    near("n = 8 transmissions, simulated", sends, 1 / (1 - parity_outcomes(8, 0.01)[1]), 0.005)
    check("n = 8: wrong with no check, 1 - c", 1 - 0.99 ** 8, "0.077255")
    clean, alarm, slip = parity_outcomes(1000, 1e-4)
    check("n = 1000: one accepted frame in this many is wrong", (1 - alarm) / slip, "201")
    check("n = 1000: frames accepted per second when 1000 are sent", 1000 * (1 - alarm), "909")
    check("n = 1000: wrong frames accepted per second when 1000 are sent", 1000 * slip, "4.5")
    check("the same, as 909 accepted frames over 201", 1000 * (1 - alarm) / 201, "4.5")
    check("n = 1000: share of corrupted frames that slip", slip / (1 - clean), "0.0475")
    check("2^-32", 2 ** -32, "2.33e-10")
    P = frame_error(1000, 1e-4)
    undetected = P * 2 ** -32
    check("heuristic undetected probability with 32 check bits", undetected, "2.22e-11")
    check("mean time to an undetected frame, years", 1 / (1000 * undetected) / 3600 / 24 / 365.25, "1.43")
    # The exact count behind the heuristic: 2^k - 1 of the 2^n - 1 nonzero patterns are codewords.
    nn, kk = 1000, 968
    frac = (2 ** kk - 1) / (2 ** nn - 1)
    near("escape fraction (2^k - 1)/(2^n - 1) against 2^-r", frac, 2.0 ** -(nn - kk), 1e-9)


def check_hybrid() -> None:
    section("Hybrid ARQ")
    n, t = 1000, 3
    rate_crc, rate_h = 968 / 1000, 938 / 1000
    check("type-I frame efficiency", rate_h, "0.938")
    p = 1e-3
    q = (1 - p) ** n
    check("p = 1e-3: clean frame probability", q, "0.367695")
    check("p = 1e-3: frame error probability", 1 - q, "0.632305")
    pf = hybrid_failure(n, t, p)
    check("p = 1e-3: probability of four or more errors", pf, "0.018927")
    check("p = 1e-3: hybrid efficiency", (1 - pf) * rate_h, "0.92025")
    check("p = 1e-3: selective repeat efficiency", eta_sr(1 - q, rate_crc), "0.35593")
    p = 1e-4
    pf = hybrid_failure(n, t, p)
    check("p = 1e-4: probability of four or more errors", pf, "3.8e-6")
    check("p = 1e-4: hybrid efficiency", (1 - pf) * rate_h, "0.93800")
    cross = hybrid_crossover(n, rate_crc, rate_h, t)
    check("crossover bit error probability", cross, "3.1e-5")
    check("BCH parity bits m t for m = 10, t = 3", 10 * 3, "30")
    check("data bits in the hybrid frame", 1000 - 32 - 30, "938")


def check_ber_wer() -> None:
    section("Bit and word error rates")
    p = 0.01
    dist = [binom_pmf(4, j, p) for j in range(5)]
    for j, stated in enumerate(("0.96059601", "0.03881196", "0.00058806", "3.96e-6", "1e-8")):
        check(f"uncoded 4-bit word: P(J = {j})", dist[j], stated)
    pw = 1 - dist[0]
    ej = sum(j * dist[j] for j in range(5))
    check("uncoded word error rate", pw, "0.03940399")
    check("expected wrong bits E[J]", ej, "0.04")
    check("uncoded bit error rate E[J]/4", ej / 4, "0.01")
    check("lower bound P_w/k", pw / 4, "0.0098510")
    for i, stated in ((2, "2.97e-4"), (3, "1e-6")):
        check(f"R3 term i = {i}", binom_pmf(3, i, p), stated)
    r3 = word_error(3, 1, p)
    check("R3 word error rate", r3, "2.98e-4")
    for i, stated in ((3, "9.801e-6"), (4, "4.95e-8"), (5, "1e-10")):
        check(f"R5 term i = {i}", binom_pmf(5, i, p), stated)
    r5 = word_error(5, 2, p)
    check("R5 word error rate", r5, "9.8506e-6")
    check("R3 over R5", r3 / r5, "30.3")
    # Exhaustive enumeration of every pattern for R3 and R5.
    for n, value in ((3, r3), (5, r5)):
        total = sum(p ** sum(e) * (1 - p) ** (n - sum(e)) for e in product((0, 1), repeat=n) if sum(e) > n // 2)
        near(f"R{n} word error rate, all {2 ** n} patterns enumerated", total, value, 1e-9)
    block = 1 - (1 - r3) ** 4
    check("4 bits each sent with R3: word error rate", block, "0.0011915")
    check("4 bits each sent with R3: lower bound P_w/4", block / 4, "2.9787e-4")
    check("uncoded word errors per R3 word error, 4-bit message", pw / block, "33.1")


def check_coding_gain() -> None:
    section("Coding gain")
    x, unc = required_uncoded(1e-6)
    check("uncoded: Q argument at 1e-6", x, "4.7534")
    check("uncoded: Eb/N0 ratio at 1e-6", x * x / 2, "11.2975")
    check("uncoded: Eb/N0 dB at 1e-6", unc, "10.530")
    p3, x3, es3, eb3 = required_repetition(3, 1e-6)
    check("R3: channel flip probability", p3, "5.7746e-4")
    check("R3: Q argument", x3, "3.2498")
    check("R3: Es/N0 dB", es3, "7.227")
    check("R3: Eb/N0 dB", eb3, "11.998")
    check("R3: coding gain dB", unc - eb3, "-1.468")
    p5, x5, es5, eb5 = required_repetition(5, 1e-6)
    check("R5: channel flip probability", p5, "4.6524e-3")
    check("R5: Q argument", x5, "2.6006")
    check("R5: Es/N0 dB", es5, "5.291")
    check("R5: Eb/N0 dB", eb5, "12.281")
    check("R5: coding gain dB", unc - eb5, "-1.751")
    check("10 log10 3", db(3), "4.771")
    check("10 log10 5", db(5), "6.990")
    _, unc5 = required_uncoded(1e-5)
    check("uncoded Eb/N0 dB at 1e-5 (chapter 3)", unc5, "9.588")
    check("R3 gain at 1e-5 (chapter 3)", unc5 - required_repetition(3, 1e-5)[3], "-1.441")
    check("R5 gain at 1e-5", unc5 - required_repetition(5, 1e-5)[3], "-1.718")
    check("asymptotic hard gain R3", asymptotic_hard_gain(3), "-1.761")
    check("asymptotic hard gain R5", asymptotic_hard_gain(5), "-2.218")
    check("asymptotic hard gain R7", asymptotic_hard_gain(7), "-2.430")
    check("limit 10 log10(1/2)", db(0.5), "-3.010")
    # The finite-target loss creeps toward the asymptote as the target falls.
    prev = 0.0
    for target in (1e-5, 1e-6, 1e-9, 1e-12, 1e-15):
        loss = required_uncoded(target)[1] - required_repetition(3, target)[3]
        print(f"  R3 gain at {target:g}: {loss:.4f} dB")
        if not (loss < prev and loss > asymptotic_hard_gain(3)):
            print("MISMATCH: the R3 gain does not fall monotonically toward its asymptote")
            sys.exit(1)
        prev = loss
    # Soft-decision R_n equals uncoded BPSK.
    rng = random.Random(3)
    sim = soft_rep_ber_simulated(3, 2.0, 200_000, rng)
    near("soft R3 at 2 dB, simulated over 200,000 bits", sim, uncoded_ber(2.0), 0.03)
    check("rate-1/2 soft headroom at 1e-6, dB", unc - 0.19, "10.34")
    check("rate-1/2 hard headroom at 1e-6, dB", unc - 1.77, "8.76")
    check("BER of hard R3 at 12 dB", hard_rep_ber(3, 12.0), "9.9e-7")


def check_monte_carlo() -> None:
    section("Measuring error rates by simulation")
    exact = word_error(3, 1, MAIN_P)
    check("exact R3 word error at p = 0.1", exact, "0.028")
    outcomes = main_run()
    k = sum(outcomes)
    ph = k / MAIN_TRIALS
    same("failures K in 10,000 seeded trials", k, 266)
    check("estimate p-hat", ph, "0.0266")
    se = sqrt(ph * (1 - ph) / MAIN_TRIALS)
    check("estimated standard error", se, "0.0016091")
    check("true standard error", sqrt(exact * (1 - exact) / MAIN_TRIALS), "0.0016497")
    check("1.96 times the standard error", Z95 * se, "0.0031539")
    lo, hi = wald(k, MAIN_TRIALS)
    check("Wald lower", lo, "0.023446")
    check("Wald upper", hi, "0.029754")
    wlo, whi = wilson(k, MAIN_TRIALS)
    check("Wilson lower", wlo, "0.023623")
    check("Wilson upper", whi, "0.029940")
    check("Wilson centre", (ph + Z95 ** 2 / (2 * MAIN_TRIALS)) / (1 + Z95 ** 2 / MAIN_TRIALS), "0.026782")
    check("estimate below truth, in standard errors", (exact - ph) / sqrt(exact * (1 - exact) / MAIN_TRIALS), "0.85")
    if not (lo < exact < hi and wlo < exact < whi):
        print("MISMATCH: an interval misses the exact value")
        sys.exit(1)
    # Running estimate checkpoints for the chart.
    for n, stated in ((100, "0.01"), (300, "0.02"), (1000, "0.022"), (3000, "0.026333")):
        check(f"running estimate after {n} trials", sum(outcomes[:n]) / n, stated)
    # Trials for a target precision.
    check("1 - 2Q(1.96)", 1 - 2 * q_func(Z95), "0.95")
    check("(1.96/0.2)^2", (Z95 / 0.2) ** 2, "96.04")
    check("relative half-width with 100 errors", Z95 / sqrt(100), "0.196")
    check("trials for 20% at p = 0.028", trials_needed(exact, 0.2), "3,334")
    same("ceil(100/p) at p = 0.028", ceil(100 / exact), 3572)
    check("trials for 20% at p = 1e-3", trials_needed(1e-3, 0.2), "95,944")
    check("trials for 20% at p = 1e-6", trials_needed(1e-6, 0.2), "9.604e7")
    for p, n, stated in ((0.028, 3572, "0.958"), (1e-3, 100_000, "0.960"), (1e-6, 100_000_000, "0.960")):
        check(f"exact chance of landing within 20%, p = {p:g}, N = {n:,}", window_prob(n, p, 0.2), stated)
    check("the same at p = 0.028, to five places", window_prob(3572, exact, 0.2), "0.95767")
    for p, n in ((1e-3, 100_000), (1e-6, 100_000_000)):
        same(f"N p is exactly 100 at p = {p:g}", round(n * p, 9), 100.0)
        check(f"exact chance without the boundary counts K = 80 and 120, p = {p:g}",
              window_prob(n, p, 0.2) - exp(binom_log_pmf(n, 80, p)) - exp(binom_log_pmf(n, 120, p)), "0.949")
    rng = random.Random(99)
    hits = 0
    for _ in range(400):
        kk = sum(simulate_repetition(3, MAIN_P, 3572, rng))
        hits += abs(kk / 3572 - exact) <= 0.2 * exact
    same("seeded experiments of 3,572 trials within 20%, of 400", hits, 381)
    near("that share against the exact chance", hits / 400, window_prob(3572, exact, 0.2), 0.03)
    same("Chebyshev trials for 20% at 95%, p = 0.028, rounded up",
         ceil((1 - exact) / (0.2 ** 2 * exact * 0.05)), 17358)
    check("Chebyshev factor 1/(0.05 r^2) at r = 0.2", 1 / (0.05 * 0.2 ** 2), "500")
    check("Chebyshev over the Gaussian count", (1 / (0.04 * 0.05)) / (Z95 / 0.2) ** 2, "5.2")
    # Zero failures.
    check("ln 20", log(20), "2.996")
    check("exact 95% upper bound after 10^6 clean trials", 1 - 0.05 ** (1 / 1e6), "2.9957e-6")
    # The five points of the second chart.
    for p, ex, kk, n, wl, wh in mc_points():
        print(f"  p = {p}: exact {ex:.6g}, K = {kk}, Wilson [{wl:.6g}, {wh:.6g}]")
        if not wl < ex < wh:
            print(f"MISMATCH: the interval at p = {p} misses the exact value")
            sys.exit(1)
    pts = mc_points()
    same("failures at p = 0.02", pts[0][2], 9)
    check("exact R3 word error at p = 0.02", pts[0][1], "0.001184")
    same("failures at p = 0.3", pts[4][2], 2155)
    check("exact R3 word error at p = 0.3", pts[4][1], "0.216")


def check_wilson_exercise() -> None:
    section("Exercise: two failures in 1,000 trials")
    lo, hi = wald(2, 1000)
    check("Wald lower", lo, "-0.000769")
    check("Wald upper", hi, "0.004769")
    wlo, whi = wilson(2, 1000)
    check("Wilson lower", wlo, "0.00055")
    check("Wilson upper", whi, "0.00726")


def check_reliability() -> None:
    section("Reliability in FIT, SDC and DUE")
    raw = 40_000
    uncorrected = raw * 0.0025
    due = uncorrected * 0.95
    sdc = uncorrected * 0.05
    check("raw flip rate per hour", raw * 1e-9, "4e-5")
    check("uncorrected FIT", uncorrected, "100")
    check("DUE FIT", due, "95")
    check("SDC FIT", sdc, "5")
    check("raw MTBF hours", mtbf_hours(raw), "25,000")
    check("raw MTBF years", mtbf_hours(raw) / HOURS_PER_YEAR, "2.85")
    check("DUE MTBF hours", mtbf_hours(due), "1.0526e7")
    check("DUE MTBF years", mtbf_hours(due) / HOURS_PER_YEAR, "1,202")
    check("SDC MTBF hours", mtbf_hours(sdc), "2e8")
    check("SDC MTBF years", mtbf_hours(sdc) / HOURS_PER_YEAR, "22,831")
    servers = 100_000
    check("fleet raw MTBF minutes", mtbf_hours(raw * servers) * 60, "15")
    check("fleet DUE MTBF hours", mtbf_hours(due * servers), "105.3")
    check("fleet DUE MTBF days", mtbf_hours(due * servers) / 24, "4.39")
    check("fleet SDC MTBF hours", mtbf_hours(sdc * servers), "2,000")
    check("fleet SDC MTBF days", mtbf_hours(sdc * servers) / 24, "83.3")
    lam = sdc * servers * 1e-9
    check("fleet SDC rate per hour", lam, "5e-4")
    check("fleet lambda t over a year", lam * HOURS_PER_YEAR, "4.38")
    check("fleet chance of at least one SDC in a year", 1 - exp(-lam * HOURS_PER_YEAR), "0.98747")
    check("one server's chance of an SDC in a year", 1 - exp(-sdc * 1e-9 * HOURS_PER_YEAR), "4.38e-5")
    # The survival limit (1 - lambda t/m)^m -> exp(-lambda t), checked at a large m.
    lt = 4.38
    near("(1 - 4.38/m)^m at m = 10^7 against exp(-4.38)", (1 - lt / 1e7) ** 1e7, exp(-lt), 1e-5)
    check("exercise: 500 FIT MTBF hours", mtbf_hours(500), "2e6")
    check("exercise: 500 FIT MTBF years", mtbf_hours(500) / HOURS_PER_YEAR, "228.3")
    check("exercise: 20,000 devices at 500 FIT, hours", mtbf_hours(500 * 20_000), "100")


def check_exercises() -> None:
    section("Exercises")
    n, R, t_rt, p = 1000, 10e6, 0.05e-3, 1e-5
    a = t_rt / (n / R)
    P = frame_error(n, p)
    check("short link: a", a, "0.5")
    check("short link: stop-and-wait efficiency", eta_sw(P, a, 0.968), "0.63891")
    clean, alarm, slip = parity_outcomes(4, 0.1)
    check("n = 4, p = 0.1: clean", clean, "0.6561")
    check("n = 4, p = 0.1: alarm", alarm, "0.2952")
    check("n = 4, p = 0.1: slip", slip, "0.0487")
    check("n = 4, p = 0.1: residual", slip / (1 - alarm), "0.069098")
    check("n = 4, p = 0.1: transmissions", 1 / (1 - alarm), "1.41884")
    check("trials for 10% at 95%, p = 1e-5", trials_needed(1e-5, 0.1), "3.8416e7")
    check("(1.96/0.1)^2", (Z95 / 0.1) ** 2, "384.16")
    pw = word_error(5, 2, 0.01)
    check("trials to size an R5 run at p = 0.01, 100/P_w", 100 / pw, "1.015e7")
    check("channel bits for that run", 5 * 100 / pw, "5.08e7")


def main() -> None:
    print("Chapter 9: Error-Control Strategies and Performance")
    check_arq()
    check_satellite()
    check_residual()
    check_hybrid()
    check_ber_wer()
    check_coding_gain()
    check_monte_carlo()
    check_wilson_exercise()
    check_reliability()
    check_exercises()
    print("\nAll chapter 9 numbers reproduced.")


if __name__ == "__main__":
    main()
