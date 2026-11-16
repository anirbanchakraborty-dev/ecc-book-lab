"""Reproduce every number in chapter 28 of Error Correction Coding, "Verifying ECC Hardware".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch28_verifying_ecc_hardware.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Whole numbers, bit patterns and decoder outputs are compared exactly, and every
fractional value is printed with a fixed number of places, so its comparison is
exact too.

The running design is the Hsiao (13,8) SEC-DED decoder. Its check matrix is
taken from python/ch20_secded_and_hsiao_codes.py and compared with the columns
the chapter prints. The decoder model below is written from that matrix,
separately from the cocotb model in tb/secded/secded_model.py, and it plants the
same four bugs the SystemVerilog decoder plants through its MUTANT parameter.
tb/secded/run_secded.py imports mutant_counts() from here and checks that the
two models and the hardware agree. Every count over the design envelope comes
from all 23,552 cases, not from a sample, and the random campaigns replay
random.Random(28).random(), whose sequence the Python documentation guarantees
across versions, so the Octave twin can draw the same numbers.

A column or a syndrome of r rows is held as an integer whose bit i is row i and
printed with row 0 first. A data word is printed as an 8-bit binary number with
d7 on the left. Codeword bits 0 to 7 are d0 to d7 and bits 8 to 12 are c0 to c4.
"""

import random
import sys
from itertools import combinations
from math import ceil, comb, log

from ch20_secded_and_hsiao_codes import hsiao_13_8

K = 8                        # data bits
N = 13                       # codeword bits
R = 5                        # check bits
WORDS = 2**K                 # 256 data words
RATE = 10**9                 # assumed checks per second, an assumption, not a measurement
YEAR = 31_557_600            # seconds in a year of 365.25 days
SEED = 28                    # the campaign seed, fixed before any campaign was run
T_CAMPAIGN = 100             # draws in the seeded campaign
R_CAMPAIGNS = 10_000         # campaigns in the repeated experiment
TARGET = 1e-6                # the miss probability a campaign is sized for


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


# ---------------------------------------------------------------- the code

def weight(v: int) -> int:
    return bin(v).count("1")


def bits(v: int, r: int) -> str:
    """A column or syndrome as a string, row 0 first."""
    return "".join("1" if (v >> i) & 1 else "0" for i in range(r))


def word(v: int) -> str:
    """A data word as an 8-bit binary number, d7 on the left."""
    return format(v, "08b")


COLUMNS = hsiao_13_8()       # d0 .. d7, then c0 .. c4


def encode(d: int) -> int:
    check = 0
    for j in range(K):
        if (d >> j) & 1:
            check ^= COLUMNS[j]
    return d | (check << K)


def envelope() -> list[int]:
    """The 92 patterns of weight at most 2: zero, singles by position, doubles in combinations order."""
    return [0] + [1 << i for i in range(N)] + [(1 << i) | (1 << j) for i, j in combinations(range(N), 2)]


def decoder_columns(mutant: int) -> list[int]:
    """The columns the decoder's XOR trees compute. Mutant 1 leaves received bit 6 out of row 1."""
    cols = list(COLUMNS)
    if mutant == 1:
        cols[6] &= ~(1 << 1)
    return cols


def decode(received: int, mutant: int = 0) -> tuple[int, int, int]:
    """(data, corrected, uncorrectable) for one received word, with one bug planted."""
    s = 0
    for j, c in enumerate(decoder_columns(mutant)):
        if (received >> j) & 1:
            s ^= c
    position = COLUMNS.index(s) if s in COLUMNS else None
    data = received & 0xFF
    if position is not None and position < K:
        target = position
        if mutant == 3 and position == 5:
            target = 4                                # d5's syndrome flips d4
        if mutant == 4 and position == 0 and (received >> 7) & 1:
            target = None                             # d0's flip suppressed while bit 7 is 1
        if target is not None:
            data ^= 1 << target
    corrected = int(position is not None)
    uncorrectable = int(s != 0 and position is None)
    if mutant == 2 and weight(s) != 2:
        uncorrectable = 0                             # raised only for two-1 syndromes
    return data, corrected, uncorrectable


def spec(d: int, e: int) -> tuple[int, int, int]:
    w = weight(e)
    return (d, 0, 0) if w == 0 else (d, 1, 0) if w == 1 else (d ^ (e & 0xFF), 0, 1)


def judge(d: int, e: int, mutant: int) -> tuple[bool, bool]:
    """(wrong, silent) for one case of the envelope."""
    out = decode(encode(d) ^ e, mutant)
    wrong = out != spec(d, e)
    return wrong, wrong and out[0] != d and out[2] == 0


def wrong_cases(mutant: int) -> list[bool]:
    """For every case index k = 92 d + p, whether the mutant gets case (d, pattern p) wrong."""
    pats = envelope()
    return [judge(d, e, mutant)[0] for d in range(WORDS) for e in pats]


def mutant_counts(mutant: int) -> tuple[int, int]:
    """(m, silent): the wrong cases and the silent ones among the 23,552 envelope cases."""
    m = silent = 0
    for d in range(WORDS):
        for e in envelope():
            w, s = judge(d, e, mutant)
            m += w
            silent += s
    return m, silent


def zero_word_failures(mutant: int) -> int:
    """How many of the 92 patterns the mutant gets wrong on the all-zero codeword."""
    return sum(judge(0, e, mutant)[0] for e in envelope())


def harmonic(n: int) -> float:
    total = 0.0
    for k in range(1, n + 1):
        total += 1 / k
    return total


def campaign_size(m: int, n: int, target: float = TARGET) -> int:
    """The fewest uniform draws that push the miss probability (1 - m/n)^T to the target or below."""
    return ceil(log(target) / log(1 - m / n))


def draws(rng: random.Random, count: int) -> list[int]:
    """Case indices drawn uniformly from the 23,552 envelope cases."""
    return [int(rng.random() * (WORDS * 92)) for _ in range(count)]


def envelope_size(n: int) -> int:
    return 1 + n + comb(n, 2)


# ---------------------------------------------------------------- the checks

def check_code() -> None:
    section("The Hsiao (13,8) code")
    same("data columns d0..d7, row 0 first", ",".join(bits(c, R) for c in COLUMNS[:K]),
         "11010,11001,10110,10101,01110,01101,01011,00111")
    same("check columns c0..c4, row 0 first", ",".join(bits(c, R) for c in COLUMNS[K:]),
         "10000,01000,00100,00010,00001")
    require("every column has odd weight", all(weight(c) % 2 for c in COLUMNS))
    distance = min(w for w in range(1, 5)
                   if any(not _xor(sub) for sub in combinations(COLUMNS, w)))
    same("minimum distance (fewest columns summing to zero)", distance, 4)
    d = 0b01101001
    c = encode(d)
    same("codeword of 01101001: check bits c0..c4", bits(c >> K, R), "01001")
    out = decode(c ^ (1 << 2))
    same("single error in d2: syndrome", bits(_syn(c ^ (1 << 2)), R), "10110")
    same("single error in d2: (data, corrected, uncorrectable)", f"{word(out[0])},{out[1]},{out[2]}",
         "01101001,1,0")
    r = c ^ (1 << 2) ^ (1 << 9)
    out = decode(r)
    same("double error in d2 and c1: syndrome", bits(_syn(r), R), "11110")
    same("double error in d2 and c1: (data, corrected, uncorr.)", f"{word(out[0])},{out[1]},{out[2]}",
         "01101101,0,1")
    out = decode(c ^ (1 << 10))
    same("single error in c2: (data, corrected, uncorrectable)", f"{word(out[0])},{out[1]},{out[2]}",
         "01101001,1,0")


def _xor(cols) -> int:
    s = 0
    for c in cols:
        s ^= c
    return s


def _syn(received: int) -> int:
    return _xor(COLUMNS[j] for j in range(N) if (received >> j) & 1)


def check_envelope() -> None:
    section("Guarantees and the design envelope")
    same("patterns of weight 0, 1, 2", f"{comb(N, 0)},{comb(N, 1)},{comb(N, 2)}", "1,13,78")
    same("patterns in the envelope", len(envelope()), 92)
    same("patterns outside the envelope", 2**N - 92, 8100)
    same("patterns of weight 3 to 13", ",".join(str(comb(N, w)) for w in range(3, N + 1)),
         "286,715,1287,1716,1716,1287,715,286,78,13,1")
    singles = [COLUMNS[i] for i in range(N)]
    require("the 13 single-error syndromes are distinct and odd",
            len(set(singles)) == N and all(weight(s) % 2 for s in singles))
    doubles = [COLUMNS[i] ^ COLUMNS[j] for i, j in combinations(range(N), 2)]
    same("double-error syndromes of weight 2 and of weight 4",
         f"{sum(weight(s) == 2 for s in doubles)},{sum(weight(s) == 4 for s in doubles)}", "52,26")
    for label, lo_i, hi_i, lo_j, hi_j, stated in (("check pairs", K, N, K, N, "10,0"),
                                                   ("data-and-check pairs", 0, K, K, N, "24,16"),
                                                   ("data pairs", 0, K, 0, K, "18,10")):
        split = [COLUMNS[i] ^ COLUMNS[j] for i in range(lo_i, hi_i) for j in range(lo_j, hi_j) if i < j]
        same(f"{label}: syndromes of weight 2 and of weight 4",
             f"{sum(weight(s) == 2 for s in split)},{sum(weight(s) == 4 for s in split)}", stated)
    left_out = [sum(not (COLUMNS[j] >> row) & 1 for j in range(K)) for row in range(R)]
    same("data columns leaving out rows 0 to 4", ",".join(map(str, left_out)), "4,3,3,3,3")
    same("data pairs sharing two rows, sum of C(left out, 2)", sum(comb(x, 2) for x in left_out), 18)
    require("the correct decoder meets the guarantees on all 23,552 cases",
            mutant_counts(0) == (0, 0))
    mis = flagged = 0
    for t in combinations(range(N), 3):
        out = decode(sum(1 << i for i in t))
        mis += out[1]
        flagged += out[2]
    same("triples: total, miscorrected, flagged", f"{comb(N, 3)},{mis},{flagged}", "286,220,66")
    same("triples: percent miscorrected", f"{100 * mis / comb(N, 3):.1f}", "76.9")


def check_counting() -> None:
    section("Counting the cases")
    same("(13,8): data words", WORDS, 256)
    same("(13,8): error patterns", 2**N, 8192)
    same("(13,8): envelope cases, 256 x 92", WORDS * 92, 23552)
    same("(13,8): all cases, 256 x 8192", WORDS * 2**N, 2097152)
    same("(13,8): input combinations of the formal wrapper, 2^21", 2**(K + N), 2**21)
    same("(72,64): envelope patterns 1 + 72 + 2556", f"{1 + 72 + comb(72, 2)}={1}+{72}+{comb(72, 2)}",
         "2629=1+72+2556")
    same("(72,64): data words", 2**64, 18446744073709551616)
    same("(72,64): error patterns", 2**72, 4722366482869645213696)
    same("(72,64): envelope cases, 2629 x 2^64", 2629 * 2**64, 48496490169782411198464)
    same("(72,64): all cases, 2^136", 2**136, 87112285931760246646623899502532662132736)
    section("Times at an assumed 10^9 checks per second")
    same("(13,8) envelope, microseconds", f"{WORDS * 92 * 1e6 / RATE:.1f}", "23.6")
    same("(13,8) all cases, milliseconds", f"{WORDS * 2**N * 1e3 / RATE:.1f}", "2.1")
    same("(72,64) envelope on one codeword, microseconds", f"{2629 * 1e6 / RATE:.1f}", "2.6")
    same("(72,64) all patterns on one codeword, whole years", 2**72 // (RATE * YEAR), 149642)
    same("(72,64) envelope cases, millions of years", f"{2629 * 2**64 / (RATE * YEAR) / 1e6:.2f}", "1.54")
    same("(72,64) all cases, years", f"{2**136 / (RATE * YEAR):.2e}", "2.76e+24")
    section("Envelope and all-pattern sweeps for four Hsiao codes")
    for (n, k), pats, env, full in (((13, 8), 92, 23552, 2097152),
                                    ((22, 16), 254, 16646144, 274877906944),
                                    ((39, 32), 781, 3354369458176, 2361183241434822606848),
                                    ((72, 64), 2629, 48496490169782411198464,
                                     87112285931760246646623899502532662132736)):
        same(f"({n},{k}): envelope patterns", envelope_size(n), pats)
        same(f"({n},{k}): envelope cases", 2**k * envelope_size(n), env)
        same(f"({n},{k}): all cases", 2**k * 2**n, full)
    same("(39,32) envelope sweep at 10^9 per second, minutes", f"{2**32 * 781 / RATE / 60:.1f}", "55.9")


def check_reduction() -> None:
    section("The reduction theorem, checked by brute force")
    pats = envelope()
    for mutant in (0, 2, 3):
        ok = True
        for d in range(WORDS):
            c = encode(d)
            for e in pats:
                zero = decode(e, mutant)
                if decode(c ^ e, mutant) != (d ^ zero[0], zero[1], zero[2]):
                    ok = False
        require(f"mutant {mutant}: output on c + e is d + output on e, every case", ok)
    for mutant, stated in ((1, 9856), (4, 128)):
        broken = sum(decode(encode(d) ^ e, mutant) != (d ^ decode(e, mutant)[0],) + decode(e, mutant)[1:]
                     for d in range(WORDS) for e in pats)
        same(f"mutant {mutant}: cases where that identity fails", broken, stated)
    require("mutant 1: some codeword has a nonzero syndrome",
            any(decode(encode(d), 1) != (d, 0, 0) for d in range(WORDS)))
    for mutant, stated in ((0, 0), (1, 12), (2, 26), (3, 1), (4, 0)):
        same(f"mutant {mutant}: patterns wrong on the zero codeword", zero_word_failures(mutant), stated)
    flips_d6 = [e for e in pats if (e >> 6) & 1]
    wrong_d6 = [e for e in flips_d6 if judge(0, e, 1)[0]]
    same("mutant 1, zero codeword: patterns flipping d6, wrong ones", f"{len(flips_d6)},{len(wrong_d6)}", "13,12")
    require("mutant 1, zero codeword: every wrong pattern flips d6",
            all((e >> 6) & 1 for e in pats if judge(0, e, 1)[0]))
    e = (1 << 6) | (1 << K)
    mcols = decoder_columns(1)
    out = decode(e, 1)
    same("mutant 1, zero codeword, d6 and c0: mutant syndrome, output",
         f"{bits(_xor(mcols[j] for j in range(N) if (e >> j) & 1), R)},{word(out[0])},{out[1]},{out[2]}",
         "10011,01000000,0,1")


def check_campaign_math() -> None:
    section("Coupon collector and miss probabilities")
    same("H_92", f"{harmonic(92):.4f}", "5.1044")
    same("92 H_92, expected draws to see all 92 patterns", f"{92 * harmonic(92):.1f}", "469.6")
    same("H_2629", f"{harmonic(2629):.4f}", "8.4518")
    same("2629 H_2629, expected draws to see all 2629 patterns", f"{2629 * harmonic(2629):.1f}", "22219.7")
    same("13 H_13, expected draws to see all 13 singles", f"{13 * harmonic(13):.1f}", "41.3")
    same("T for miss <= 10^-6, m = 1, N = 92", campaign_size(1, 92), 1265)
    same("T for miss <= 10^-6, m = 1, N = 2629", campaign_size(1, 2629), 36315)
    same("T for miss <= 10^-6, m = 1, N = 23552", campaign_size(1, 23552), 325376)
    same("approximation N ln(10^6), N = 92", f"{92 * log(1e6):.0f}", "1271")
    same("approximation N ln(10^6), N = 2629", f"{2629 * log(1e6):.0f}", "36321")
    same("T for miss <= 10^-6, m = 26, N = 92", campaign_size(26, 92), 42)
    same("T for miss <= 10^-6, m = 1, N = 13", campaign_size(1, 13), 173)
    same("one-case campaign over the 23,552-case sweep", f"{campaign_size(1, WORDS * 92) / (WORDS * 92):.1f}", "13.8")


def check_mutants() -> None:
    section("The planted bugs over all 23,552 envelope cases")
    stated = {0: (0, 0), 1: (9856, 6528), 2: (6656, 6656), 3: (256, 256), 4: (128, 128)}
    for mutant in range(5):
        m, silent = mutant_counts(mutant)
        same(f"mutant {mutant}: wrong cases m, silent cases", f"{m},{silent}", f"{stated[mutant][0]},{stated[mutant][1]}")
    same("mutant 1: wrong patterns per word, d6 = 0 and d6 = 1",
         f"{sum(judge(0, e, 1)[0] for e in envelope())},{sum(judge(1 << 6, e, 1)[0] for e in envelope())}",
         "12,65")
    for mutant, p_stated, miss_stated, t_stated in ((1, "0.4185", "0.000", 26), (2, "0.2826", "0.000", 42),
                                                     (3, "0.0109", "0.335", 1265), (4, "0.0054", "0.580", 2536)):
        m = mutant_counts(mutant)[0]
        p = m / (WORDS * 92)
        same(f"mutant {mutant}: chance one draw hits it", f"{p:.4f}", p_stated)
        same(f"mutant {mutant}: miss probability at T = 100", f"{(1 - p) ** T_CAMPAIGN:.3f}", miss_stated)
        same(f"mutant {mutant}: T for miss <= 10^-6", campaign_size(m, WORDS * 92), t_stated)
    m4 = mutant_counts(4)[0]
    require("mutant 4: m divides N", (WORDS * 92) % m4 == 0)
    same("mutant 4: m/N as 1/x", f"1/{WORDS * 92 // m4}", "1/184")
    same("mutant 4: sweep over campaign, N / T", f"{WORDS * 92 / campaign_size(m4, WORDS * 92):.2f}", "9.29")
    same("mutant 2: wrong cases with wt(e) < 2", sum(judge(d, e, 2)[0] for d in range(WORDS)
                                                    for e in envelope() if weight(e) < 2), 0)
    same("mutant 3: wrong cases with wt(e) < 2", sum(judge(d, e, 3)[0] for d in range(WORDS)
                                                    for e in envelope() if weight(e) < 2), 256)
    same("cases with wt(e) < 2", WORDS * 14, 3584)


def check_campaigns() -> None:
    section(f"Seeded campaigns, seed {SEED}")
    wrong = {mutant: wrong_cases(mutant) for mutant in range(1, 5)}
    rng = random.Random(SEED)
    first = draws(rng, T_CAMPAIGN)
    same("first draw as (data word, pattern index)", f"{word(first[0] // 92)},{first[0] % 92}", "00011100,84")
    hits = ",".join(str(sum(wrong[mutant][k] for k in first)) for mutant in range(1, 5))
    same("one campaign of 100: hits on mutants 1, 2, 3, 4", hits, "33,29,1,0")
    same("one campaign of 100: distinct patterns, distinct words",
         f"{len({k % 92 for k in first})},{len({k // 92 for k in first})}", "59,77")
    rng = random.Random(SEED)
    missed = {mutant: 0 for mutant in range(1, 5)}
    for _ in range(R_CAMPAIGNS):
        batch = draws(rng, T_CAMPAIGN)
        for mutant in range(1, 5):
            if not any(wrong[mutant][k] for k in batch):
                missed[mutant] += 1
    same("10,000 campaigns of 100: missed mutants 1, 2, 3, 4",
         ",".join(str(missed[mutant]) for mutant in range(1, 5)), "0,0,3401,5883")
    same("mutant 3: fraction missed", f"{missed[3] / R_CAMPAIGNS:.4f}", "0.3401")
    same("mutant 4: fraction missed", f"{missed[4] / R_CAMPAIGNS:.4f}", "0.5883")
    rng = random.Random(SEED)
    seen = set()
    all_at = 0
    for t, k in enumerate(draws(rng, 700), start=1):
        seen.add(k % 92)
        if not all_at and len(seen) == 92:
            all_at = t
    same("coupon campaign, seed 28: draw that completes all 92 patterns", all_at, 454)


def check_formal() -> None:
    section("The formal counterexamples the chapter reads")
    d, e = 0b10111100, 1 << 5
    out = decode(encode(d) ^ e, 3)
    same("mutant 3 at d = 10111100, e = d5: decoder output", f"{word(out[0])},{out[1]},{out[2]}",
         "10001100,1,0")
    same("mutant 3 at that case: wrong, silent", ",".join(str(int(x)) for x in judge(d, e, 3)), "1,1")
    same("syndrome of that received word", bits(_syn(encode(d) ^ e), R), "01101")
    for mutant, d, e, pattern, syndrome, output, guarantee in (
            (1, 0b11111100, 0b11, "0000000000011", "00011", "10111111,1,0", "11111111,0,1"),
            (2, 0b01000000, (1 << 2) | (1 << 9), "0001000000100", "11110", "01000100,0,0", "01000100,0,1"),
            (4, 0b10011001, 1, "0000000000001", "11010", "10011000,1,0", "10011001,1,0")):
        r = encode(d) ^ e
        out = decode(r, mutant)
        want = spec(d, e)
        same(f"mutant {mutant} at d = {word(d)}: pattern, position 12 first", format(e, "013b"), pattern)
        same(f"mutant {mutant} at that case: correct syndrome", bits(_syn(r), R), syndrome)
        same(f"mutant {mutant} at that case: decoder output", f"{word(out[0])},{out[1]},{out[2]}", output)
        same(f"mutant {mutant} at that case: the guarantee", f"{word(want[0])},{want[1]},{want[2]}", guarantee)
        same(f"mutant {mutant} at that case: wrong, silent", ",".join(str(int(x)) for x in judge(d, e, mutant)),
             "1,1")
    r = encode(0b11111100) ^ 0b11
    mcols = decoder_columns(1)
    same("mutant 1 at d = 11111100: mutant syndrome", bits(_xor(mcols[j] for j in range(N) if (r >> j) & 1), R),
         "01011")


def main() -> None:
    check_code()
    check_envelope()
    check_counting()
    check_reduction()
    check_campaign_math()
    check_mutants()
    check_campaigns()
    check_formal()
    print("\nAll chapter 28 numbers reproduced.")


if __name__ == "__main__":
    main()
