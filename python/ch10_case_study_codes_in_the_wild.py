"""Reproduce every number in chapter 10 of Error Correction Coding, "Case Study: Codes in the Wild".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch10_case_study_codes_in_the_wild.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. The codes are built
from their definitions rather than typed in. The Hamming, Reed-Muller, Hsiao and
5G (32,K) codes are enumerated codeword by codeword, the free distance of the
CCSDS convolutional code is found by a shortest-path search over its trellis,
the CCSDS Reed-Solomon generator is multiplied out over GF(256), the compact
disc's burst bound is checked against a frame-by-frame model of its delay
lines, the QR version 1 symbol is laid out module by module, and the 5G
lifting sizes and code block parameters follow the rules of TS 38.212. A run
that reaches the final line has reproduced every number in the chapter.

The module also holds the functions that ch10_charts.py plots, so the charts
and this check share one source.
"""

import heapq
import sys
from itertools import combinations, product
from math import ceil, comb, erfc, factorial, log2, sqrt


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


def q_func(x: float) -> float:
    """The Q-function, the probability that a standard Gaussian exceeds x."""
    return 0.5 * erfc(x / sqrt(2.0))


def tail_at_least(n: int, t: int, p: float) -> float:
    """The probability of t or more errors among n bits of BSC(p)."""
    return sum(comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(t, n + 1))


def weight(v: int) -> int:
    return bin(v).count("1")


# Parity in early memories ----------------------------------------------------------

def even_count_probability(n: int, p: float) -> float:
    """The probability that n bits of BSC(p) suffer an even number of errors, zero included."""
    return (1 + (1 - 2 * p) ** n) / 2


def byte_parity_word(p: float) -> tuple:
    """A 72-bit word stored as eight bytes, each with its own parity bit: (alarm, silent corruption)."""
    all_even = even_count_probability(9, p) ** 8
    return 1 - all_even, all_even - (1 - p) ** 72


def secded_uncorrected(p: float) -> float:
    """The probability that a 72-bit SEC-DED word holds two or more errors, so is not corrected."""
    return tail_at_least(72, 2, p)


def check_parity() -> None:
    section("Parity in Early Memories")
    same("ones in the byte 10110010", weight(0b10110010), 4)
    same("even-parity bit for 10110010", weight(0b10110010) % 2, 0)
    same("ones after flipping one bit", weight(0b10110010 ^ 0b00001000), 5)
    check("byte parity rate 8/9", 8 / 9, "0.889")
    check("check overhead 8/64", 8 / 64, "12.5%")
    same("bits in a 64-bit word stored with byte parity", 8 * 9, 72)
    # The even-count formula against direct enumeration of every 9-bit pattern.
    p = 0.01
    direct = sum(p ** weight(e) * (1 - p) ** (9 - weight(e)) for e in range(512) if weight(e) % 2 == 0)
    near("even-count formula against enumeration, n = 9", even_count_probability(9, p), direct, 1e-12)


# Hamming's code -------------------------------------------------------------------

def hamming_encode(d3: int, d5: int, d6: int, d7: int) -> list:
    """Hamming's arrangement: check bits at positions 1, 2 and 4, data at 3, 5, 6 and 7."""
    p1 = d3 ^ d5 ^ d7
    p2 = d3 ^ d6 ^ d7
    p4 = d5 ^ d6 ^ d7
    return [p1, p2, d3, p4, d5, d6, d7]


def hamming_syndrome(word: list) -> int:
    s = 0
    for position, bit in enumerate(word, start=1):
        if bit:
            s ^= position
    return s


def uncoded_block(p: float, k: int) -> float:
    return 1 - (1 - p) ** k


def hamming74_block(p: float) -> float:
    return 1 - ((1 - p) ** 7 + 7 * p * (1 - p) ** 6)


def repetition3_block(p: float) -> float:
    bit = 3 * p ** 2 * (1 - p) + p ** 3
    return 1 - (1 - bit) ** 4


def check_hamming() -> None:
    section("Hamming's Code at Bell Labs")
    word = hamming_encode(1, 0, 1, 1)
    same("codeword for data 1011", "".join(map(str, word)), "0110011")
    bad = word[:]
    bad[4] ^= 1
    same("received word after flipping position 5", "".join(map(str, bad)), "0110111")
    same("syndrome of the received word", format(hamming_syndrome(bad), "03b"), "101")
    same("syndrome of every codeword is zero",
         all(hamming_syndrome(hamming_encode(*d)) == 0 for d in product((0, 1), repeat=4)), True)
    for pos in range(1, 8):
        flipped = word[:]
        flipped[pos - 1] ^= 1
        if hamming_syndrome(flipped) != pos:
            print(f"MISMATCH: flipping position {pos} gives syndrome {hamming_syndrome(flipped)}")
            sys.exit(1)
    same("every single flip names its own position", True, True)
    weights = {}
    for d in product((0, 1), repeat=4):
        w = sum(hamming_encode(*d))
        weights[w] = weights.get(w, 0) + 1
    same("weight distribution of the (7,4) code", dict(sorted(weights.items())), {0: 1, 3: 7, 4: 7, 7: 1})
    check("rate of the (7,4) code", 4 / 7, "0.571")
    check("rate of three-fold repetition", 1 / 3, "0.333")
    p = 0.01
    check("uncoded 4 bits at p = 0.01", uncoded_block(p, 4), "0.0394")
    check("Hamming (7,4) block error at p = 0.01", hamming74_block(p), "0.00203")
    check("three-fold repetition of 4 bits at p = 0.01", repetition3_block(p), "0.00119")
    check("uncoded over Hamming at p = 0.01", uncoded_block(p, 4) / hamming74_block(p), "19.4")
    # The block error of the (7,4) code against enumeration of all 128 error patterns
    # with a syndrome decoder.
    total = 0.0
    for e in range(128):
        pattern = [(e >> i) & 1 for i in range(7)]
        received = [a ^ b for a, b in zip(word, pattern)]
        s = hamming_syndrome(received)
        if s:
            received[s - 1] ^= 1
        if received != word:
            total += p ** sum(pattern) * (1 - p) ** (7 - sum(pattern))
    near("(7,4) block error against syndrome decoding of all 128 patterns", hamming74_block(p), total, 1e-12)
    same("data bits that r = 3 check bits protect", 2 ** 3 - 1 - 3, 4)
    same("data bits that r = 7 check bits protect", 2 ** 7 - 1 - 7, 120)
    same("fewest check bits for 64 data bits", min(r for r in range(1, 10) if 2 ** r - 1 - r >= 64), 7)


# The Mariner Reed-Muller code ------------------------------------------------------

MARINER_FIXED = "10001101110101000010010110011111"


def rm15_codeword(message: int) -> int:
    """Bit x of the codeword is a0 + a1 x1 + ... + a5 x5 mod 2, with x1 the low bit of x."""
    a0 = message & 1
    word = 0
    for x in range(32):
        bit = a0
        for i in range(5):
            if (message >> (i + 1)) & 1 and (x >> i) & 1:
                bit ^= 1
        word |= bit << x
    return word


def rm15_word_error(p: float) -> float:
    """A hard-decision decoder that corrects up to seven errors fails on eight or more."""
    return tail_at_least(32, 8, p)


def repetition5_word_error(p: float) -> float:
    bit = tail_at_least(5, 3, p)
    return 1 - (1 - bit) ** 6


def check_mariner() -> None:
    section("The Mariner Reed-Muller Code")
    words = [rm15_codeword(m) for m in range(64)]
    same("distinct codewords", len(set(words)), 64)
    dist = {}
    for w in words:
        dist[weight(w)] = dist.get(weight(w), 0) + 1
    same("weight distribution of RM(1,5)", dict(sorted(dist.items())), {0: 1, 16: 62, 32: 1})
    pair = {weight(a ^ b) for a, b in combinations(words, 2)}
    same("distances between distinct codewords", sorted(pair), [16, 32])
    same("minimum distance d", min(pair), 16)
    same("errors corrected, (d - 1) // 2", (16 - 1) // 2, 7)
    check("rate 6/32", 6 / 32, "0.1875")
    x1 = rm15_codeword(0b000010)
    x2 = rm15_codeword(0b000100)
    same("codeword for a1 = 1, first eight bits", format(x1, "032b")[::-1][:8], "01010101")
    same("codeword for a2 = 1, first eight bits", format(x2, "032b")[::-1][:8], "00110011")
    same("positions where they differ", weight(x1 ^ x2), 16)
    # Orthogonality in the plus-and-minus-one picture: the inner product of two
    # codewords of the 32 with a0 = 0 is 32 - 2 * distance, which is zero.
    ortho = all(32 - 2 * weight(rm15_codeword(2 * a) ^ rm15_codeword(2 * b)) == 0
                for a, b in combinations(range(32), 2))
    same("the 32 words with a0 = 0 are mutually orthogonal", ortho, True)
    same("length of the fixed comma-free vector", len(MARINER_FIXED), 32)
    # The fixed vector is a pseudonoise sequence of period 31 from a five-stage
    # shift register, with one bit appended.
    s = [int(c) for c in MARINER_FIXED]
    found = []
    for taps in range(1, 32):
        c = [(taps >> i) & 1 for i in range(5)]
        if c[4] == 0:
            continue
        if all(s[n] == (sum(c[i] * s[n - 1 - i] for i in range(5)) % 2) for n in range(5, 31)):
            found.append(taps)
    same("a five-stage recurrence generates its first 31 bits", bool(found), True)
    state = s[:5]
    seq = s[:5]
    c = [(found[0] >> i) & 1 for i in range(5)]
    for _ in range(62):
        seq.append(sum(c[i] * seq[-1 - i] for i in range(5)) % 2)
    period = next(t for t in range(1, 63) if seq[t:t + 31] == seq[:31])
    same("period of that recurrence", period, 31)
    same("ones among the first 31 bits", sum(s[:31]), 16)
    del state
    check("required ST_B/N_0 saved, 5.20 - 3.00", 5.20 - 3.00, "2.20")
    check("power ratio of 2.20 dB", 10 ** 0.22, "1.66")
    check("uncoded BPSK bit error at 5.20 dB", q_func(sqrt(2 * 10 ** 0.52)), "5e-3")
    check("bit rate, bits per second", 16200, "16,200")


# Voyager and concatenated codes -------------------------------------------------------

G1, G2 = 0o171, 0o133


def conv_free_distance(g1: int, g2: int, memory: int) -> int:
    """The least output weight of a path that leaves the zero state and returns to it."""
    states = 1 << memory

    def step(state: int, bit: int) -> tuple:
        reg = (bit << memory) | state
        out = weight(reg & g1) % 2 + weight(reg & g2) % 2
        return reg >> 1, out

    first, w0 = step(0, 1)
    best = {first: w0}
    heap = [(w0, first)]
    while heap:
        w, s = heapq.heappop(heap)
        if s == 0:
            return w
        if w > best.get(s, 1 << 30):
            continue
        for bit in (0, 1):
            nxt, out = step(s, bit)
            if w + out < best.get(nxt, 1 << 30):
                best[nxt] = w + out
                heapq.heappush(heap, (w + out, nxt))
    del states
    return -1


def gf_tables(poly: int) -> tuple:
    exp_t, log_t = [0] * 510, [0] * 256
    x = 1
    for i in range(255):
        exp_t[i] = x
        log_t[x] = i
        x <<= 1
        if x & 0x100:
            x ^= poly
    for i in range(255, 510):
        exp_t[i] = exp_t[i - 255]
    return exp_t, log_t


def gf_mul(a: int, b: int, exp_t: list, log_t: list) -> int:
    if a == 0 or b == 0:
        return 0
    return exp_t[log_t[a] + log_t[b]]


def ccsds_rs_generator(e: int) -> list:
    """g(x) = prod over j from 128 - E to 127 + E of (x - alpha^(11 j)), lowest degree first."""
    exp_t, log_t = gf_tables(0x187)
    g = [1]
    for j in range(128 - e, 128 + e):
        root = exp_t[(11 * j) % 255]
        nxt = [0] * (len(g) + 1)
        for i, c in enumerate(g):
            nxt[i + 1] ^= c
            nxt[i] ^= gf_mul(c, root, exp_t, log_t)
        g = nxt
    return g


VOYAGER_TABLE = [  # planet, factor in the table, expected rate, achieved rate, improvement
    ("Jupiter", 1, 115_200, 115_200, None),
    ("Saturn", 4, 29_000, 44_800, "1.5"),
    ("Uranus", 13, 9_000, 29_900, "3.3"),
    ("Neptune", 36, 3_200, 21_600, "6.8"),
]
VOYAGER_AU = {"Jupiter": 5.2, "Saturn": 10, "Uranus": 19, "Neptune": 30}


def check_voyager() -> None:
    section("Voyager and Concatenated Codes")
    same("free distance of the (171, 133) code", conv_free_distance(G1, G2, 6), 10)
    same("generators in binary", (format(G1, "07b"), format(G2, "07b")), ("1111001", "1011011"))
    check("rate of RS(255,223)", 223 / 255, "0.875")
    check("concatenated rate 223/255 x 1/2", 223 / 255 / 2, "0.437")
    check("Golay plus convolutional rate 1/2 x 1/2", 0.25, "0.25")
    same("RS check symbols", 255 - 223, 32)
    same("RS minimum distance n - k + 1", 255 - 223 + 1, 33)
    same("RS symbol errors corrected", 32 // 2, 16)
    check("RS check overhead 32/223", 32 / 223, "14.3%")
    check("RS(255,239) rate", 239 / 255, "0.937")
    exp_t, _ = gf_tables(0x187)
    order = next(i for i in range(1, 256) if exp_t[i] == 1)
    same("order of alpha in GF(256) built on 0x187", order, 255)
    g = ccsds_rs_generator(16)
    same("degree of the E = 16 generator", len(g) - 1, 32)
    same("generator coefficients read the same both ways", g == g[::-1], True)
    same("the E = 8 generator is palindromic too", ccsds_rs_generator(8) == ccsds_rs_generator(8)[::-1], True)
    check("general bound 1/16!", 1 / factorial(16), "4.78e-14")
    check("undetected word error, 1e-5 x 3e-14", 1e-5 * 3e-14, "3e-19")
    for depth, symbols, bits in ((5, 80, 640), (8, 128, 1024)):
        same(f"burst of symbols cleared at depth {depth}", 16 * depth, symbols)
        same(f"the same burst in bits at depth {depth}", 16 * depth * 8, bits)


def check_voyager_table() -> None:
    section("Voyager's Table 6-1")
    jupiter = VOYAGER_TABLE[0][2]
    exact = {"Saturn": "28,800", "Uranus": "8,862", "Neptune": "3,200"}
    for planet, factor, expected, achieved, improvement in VOYAGER_TABLE[1:]:
        check(f"{planet}: 115,200 / {factor}", jupiter / factor, exact[planet])
        near(f"{planet}: the table's rounded expected rate", expected, jupiter / factor, 0.02)
        check(f"{planet}: achieved over expected", achieved / expected, improvement)
    check("(5.2/10)^2 as one over", (10 / 5.2) ** 2, "3.70")
    check("(5.2/19)^2 as one over", (19 / 5.2) ** 2, "13.35")
    check("(5.2/30)^2 as one over", (30 / 5.2) ** 2, "33.28")
    check("table's 36 against 33.28", 36 / (30 / 5.2) ** 2 - 1, "8%")
    same("the table's 4 is (10/5)^2 exactly", (10 / 5) ** 2, 4)
    same("the table's 36 is (30/5)^2 exactly", (30 / 5) ** 2, 36)
    check("Uranus with Jupiter at 5 AU, (19/5)^2", (19 / 5) ** 2, "14.4")
    check("Neptune expected from the distances", 115_200 / (30 / 5.2) ** 2, "3,461")
    check("Neptune improvement from the distances", 21_600 / (115_200 / (30 / 5.2) ** 2), "6.2")


# The compact disc ---------------------------------------------------------------------

CHANNEL_RATE = 4_321_800
D_FRAMES = 4


def cd_max_hits(burst: int) -> int:
    """The most bytes of one C2 codeword a run of erased frames can reach, from the delay model.

    Byte i of the C2 codeword formed at frame f leaves in frame f + 4 i, so a burst that
    erases frames s to s + B - 1 erases byte i exactly when f + 4 i falls in that run.
    """
    best = 0
    for f in range(-120, 40):
        hits = sum(1 for i in range(28) if 0 <= f + D_FRAMES * i < burst)
        best = max(best, hits)
    return best


def check_cd() -> None:
    section("The Compact Disc")
    same("audio bits per second, 44,100 x 16 x 2", 44_100 * 16 * 2, 1_411_200)
    same("audio bytes per second", 44_100 * 2 * 2, 176_400)
    same("bytes in six stereo samples", 6 * 2 * 2, 24)
    same("channel bits per frame", 24 + 3 + 14 + 3 + 32 * (14 + 3), 588)
    same("frames per second, 4,321,800 / 588", CHANNEL_RATE // 588, 7_350)
    same("4,321,800 divides evenly by 588", CHANNEL_RATE % 588, 0)
    same("sections per second, 7,350 / 98", 7_350 // 98, 75)
    same("audio bytes per second from the frame rate", 24 * 7_350, 176_400)
    check("CIRC rate 24/32", 24 / 32, "0.75")
    check("share of channel bits that carry audio, 192/588", 192 / 588, "32.7%")
    same("longest C2 delay, 27 x 4 frames", 27 * D_FRAMES, 108)
    same("minimum distance of RS(28,24) and RS(32,28)", (28 - 24 + 1, 32 - 28 + 1), (5, 5))
    for b in range(0, 41):
        if cd_max_hits(b) != ceil(b / D_FRAMES):
            print(f"MISMATCH: burst {b} reaches {cd_max_hits(b)} bytes, not {ceil(b / D_FRAMES)}")
            sys.exit(1)
    same("delay model gives ceil(B/4) for every B up to 40", True, True)
    longest = max(b for b in range(0, 41) if cd_max_hits(b) <= 4)
    same("longest burst of frames that C2 can fill", longest, 16)
    check("16 frames in milliseconds", 16 / 7_350 * 1e3, "2.18")
    check("track length at 1.2 m/s, mm", 16 / 7_350 * 1.2 * 1e3, "2.6")
    check("track length at 1.4 m/s, mm", 16 / 7_350 * 1.4 * 1e3, "3.0")
    check("bound with D = 2: frames", 2 * 4, "8")
    check("bound with D = 2: milliseconds", 8 / 7_350 * 1e3, "1.09")
    same("bytes in a CD-ROM sector", 98 * 24, 2_352)


# QR codes ----------------------------------------------------------------------------

QR_DATA_CODEWORDS = {"L": 19, "M": 16, "Q": 13, "H": 9}
QR_CAPACITY = {  # DENSO WAVE's version 1 table: numeric, alphanumeric, binary, kanji
    "L": (41, 25, 17, 10), "M": (34, 20, 14, 8), "Q": (27, 16, 11, 7), "H": (17, 10, 7, 4),
}


def qr_side(version: int) -> int:
    return 21 + 4 * (version - 1)


def qr_version1_data_modules() -> int:
    n = 21
    used = set()
    for r0, c0 in ((0, 0), (0, n - 8), (n - 8, 0)):
        used |= {(r, c) for r in range(r0, r0 + 8) for c in range(c0, c0 + 8)}
    used |= {(6, c) for c in range(8, n - 8)} | {(r, 6) for r in range(8, n - 8)}
    used |= {(8, c) for c in range(0, 9) if c != 6} | {(r, 8) for r in range(0, 9) if r != 6}
    used |= {(8, c) for c in range(n - 8, n)} | {(r, 8) for r in range(n - 7, n)}
    used |= {(4 * 1 + 9, 8)}
    return n * n - len(used)


def qr_capacity(data_codewords: int) -> tuple:
    """Characters per mode in a version 1 to 9 symbol, from the mode indicator and count lengths."""
    bits = 8 * data_codewords
    numeric_bits = bits - 4 - 10
    groups, rest = divmod(numeric_bits, 10)
    numeric = 3 * groups + (2 if rest >= 7 else 1 if rest >= 4 else 0)
    alnum_bits = bits - 4 - 9
    pairs, rest = divmod(alnum_bits, 11)
    alnum = 2 * pairs + (1 if rest >= 6 else 0)
    byte = (bits - 4 - 8) // 8
    kanji = (bits - 4 - 8) // 13
    return numeric, alnum, byte, kanji


def check_qr() -> None:
    section("QR Codes")
    same("version 1 side", qr_side(1), 21)
    same("version 40 side", qr_side(40), 177)
    same("version 10 side", qr_side(10), 57)
    same("modules in version 1", 21 * 21, 441)
    same("data modules left in version 1", qr_version1_data_modules(), 208)
    same("codewords in version 1", qr_version1_data_modules() // 8, 26)
    corrected = {}
    for level, data in QR_DATA_CODEWORDS.items():
        r = 26 - data
        corrected[level] = r // 2
        same(f"level {level}: check codewords", r, {"L": 7, "M": 10, "Q": 13, "H": 17}[level])
        same(f"level {level}: capacities derived from {data} data codewords", qr_capacity(data),
             QR_CAPACITY[level])
    same("codeword errors corrected, L M Q H", tuple(corrected.values()), (3, 5, 6, 8))
    for level, stated in zip("LMQH", ("11.5%", "19.2%", "23.1%", "30.8%")):
        check(f"level {level}: share of 26 codewords corrected", corrected[level] / 26, stated)
    check("level M rate 16/26", 16 / 26, "0.615")
    same("level M byte capacity, (128 - 12) // 8", (128 - 12) // 8, 14)
    check("DENSO WAVE example: 50 of 200", 50 / 200, "25%")
    same("level H erasures filled", 26 - 9, 17)


# Memory codes after Hsiao ---------------------------------------------------------------

def is_secded(columns: list) -> bool:
    singles = set(columns)
    if len(singles) != len(columns) or 0 in singles:
        return False
    for a, b in combinations(columns, 2):
        s = a ^ b
        if s == 0 or s in singles:
            return False
    return True


def hsiao_72_64() -> list:
    units = [1 << i for i in range(8)]
    threes = [sum(1 << i for i in c) for c in combinations(range(8), 3)]
    five = 0b00011111
    fives = [((five << k) | (five >> (8 - k))) & 0xFF for k in range(8)]
    return units + threes + fives


def extended_hamming_72_64() -> list:
    """Seven Hamming rows from the lightest columns, then an all-ones overall-parity row."""
    data = []
    for w in (2, 3, 4):
        for c in combinations(range(7), w):
            if len(data) < 64:
                data.append(sum(1 << i for i in c))
    columns = data + [1 << i for i in range(7)]
    return [c | (1 << 7) for c in columns] + [1 << 7]


def hsiao_39_32() -> list:
    threes = [sum(1 << i for i in c) for c in combinations(range(7), 3)][:32]
    return [1 << i for i in range(7)] + threes


def check_hsiao() -> None:
    section("Memory Codes after Hsiao")
    h = hsiao_72_64()
    same("Hsiao (72,64): columns", len(h), 72)
    same("Hsiao (72,64): every column has odd weight", all(weight(c) % 2 for c in h), True)
    same("Hsiao (72,64): columns of weight 1, 3, 5", tuple(sum(1 for c in h if weight(c) == w) for w in (1, 3, 5)),
         (8, 56, 8))
    same("Hsiao (72,64): ones in the check matrix", sum(weight(c) for c in h), 216)
    rows = [sum((c >> r) & 1 for c in h) for r in range(8)]
    same("Hsiao (72,64): ones in each row", set(rows), {27})
    same("Hsiao (72,64): corrects one error and detects two", is_secded(h), True)
    e = extended_hamming_72_64()
    same("extended Hamming (72,64): corrects one and detects two", is_secded(e), True)
    same("extended Hamming: columns", len(e), 72)
    same("extended Hamming: data columns of weight 2, 3, 4 in the Hamming rows",
         tuple(sum(1 for c in e[:64] if weight(c & 0x7F) == w) for w in (2, 3, 4)), (21, 35, 8))
    same("extended Hamming: ones in the seven Hamming rows", sum(weight(c & 0x7F) for c in e), 186)
    same("extended Hamming: ones in all eight rows", sum(weight(c) for c in e), 258)
    same("ones Hsiao saves", 258 - 216, 42)
    check("saving as a share", 42 / 258, "16%")
    same("XOR tree depth for 27 inputs", ceil(log2(27)), 5)
    same("XOR tree depth for 72 inputs", ceil(log2(72)), 7)
    s = hsiao_39_32()
    same("Hsiao (39,32): corrects one and detects two", is_secded(s), True)
    same("Hsiao (39,32): ones", sum(weight(c) for c in s), 103)
    check("rate of the (72,64) code", 64 / 72, "0.889")
    same("double-error syndromes in the (72,64) code", comb(72, 2), 2_556)


def check_parity_against_secded() -> None:
    section("Worked example: byte parity against SEC-DED")
    p = 1e-6
    alarm, silent = byte_parity_word(p)
    check("byte parity: alarm per read", alarm, "7.2e-5")
    check("byte parity: silent corruption per read", silent, "2.88e-10")
    near("silent corruption against 8 x C(9,2) p^2", silent, 8 * comb(9, 2) * p ** 2, 1e-4)
    unc = secded_uncorrected(p)
    check("SEC-DED: two or more errors per read", unc, "2.56e-9")
    near("two or more errors against C(72,2) p^2", unc, comb(72, 2) * p ** 2, 1e-4)
    check("ratio of the two stopping rates", alarm / unc, "2.8e4")
    check("SEC-DED: three or more errors per read", tail_at_least(72, 3, p), "5.96e-14")
    same("C(9,2) and 8 x C(9,2)", (comb(9, 2), 8 * comb(9, 2)), (36, 288))


# 5G New Radio --------------------------------------------------------------------------

LIFTING_SETS = [[a * 2 ** j for j in range(10) if a * 2 ** j <= 384] for a in (2, 3, 5, 7, 9, 11, 13, 15)]
LIFTING_SIZES = sorted(z for s in LIFTING_SETS for z in s)
KCB = {1: 8448, 2: 3840}
N_OVER_Z = {1: 66, 2: 50}
K_OVER_Z = {1: 22, 2: 10}
RV_NUMERATORS = {1: (0, 17, 33, 56), 2: (0, 13, 25, 43)}


def crc_length(a: int) -> int:
    return 24 if a > 3824 else 16


def base_graph(a: int, r: float) -> int:
    return 2 if a <= 292 or (a <= 3824 and r <= 0.67) or r <= 0.25 else 1


def code_blocks(b: int, bg: int) -> int:
    return 1 if b <= KCB[bg] else ceil(b / (KCB[bg] - 24))


def kb(b: int, bg: int) -> int:
    if bg == 1:
        return 22
    return 10 if b > 640 else 9 if b > 560 else 8 if b > 192 else 6


def lifting(b: int, bg: int) -> dict:
    """One code block (C = 1): the lifting size, the code block length, fillers and buffer length."""
    zc = min(z for z in LIFTING_SIZES if kb(b, bg) * z >= b)
    k = K_OVER_Z[bg] * zc
    n = N_OVER_Z[bg] * zc
    return {"Kb": kb(b, bg), "Zc": zc, "K": k, "F": k - b, "N": n}


def k0(rv: int, bg: int, ncb: int, zc: int) -> int:
    return (RV_NUMERATORS[bg][rv] * ncb // (N_OVER_Z[bg] * zc)) * zc


SMALL_BASIS = """\
0 1 1 0 0 0 0 0 0 0 0 1
1 1 1 1 0 0 0 0 0 0 1 1
2 1 0 0 1 0 0 1 0 1 1 1
3 1 0 1 1 0 0 0 0 1 0 1
4 1 1 1 1 0 0 0 1 0 0 1
5 1 1 0 0 1 0 1 1 1 0 1
6 1 0 1 0 1 0 1 0 1 1 1
7 1 0 0 1 1 0 0 1 1 0 1
8 1 1 0 1 1 0 0 1 0 1 1
9 1 0 1 1 1 0 1 0 0 1 1
10 1 0 1 0 0 1 1 1 0 1 1
11 1 1 1 0 0 1 1 0 1 0 1
12 1 0 0 1 0 1 0 1 1 1 1
13 1 1 0 1 0 1 0 1 0 1 1
14 1 0 0 0 1 1 0 1 0 0 1
15 1 1 0 0 1 1 1 1 0 1 1
16 1 1 1 0 1 1 1 0 0 1 0
17 1 0 0 1 1 1 0 0 1 0 0
18 1 1 0 1 1 1 1 1 0 0 0
19 1 0 0 0 0 1 1 0 0 0 0
20 1 0 1 0 0 0 1 0 0 0 1
21 1 1 0 1 0 0 0 0 0 1 1
22 1 0 0 0 1 0 0 1 1 0 1
23 1 1 1 0 1 0 0 0 1 1 1
24 1 1 1 1 1 0 1 1 1 1 0
25 1 1 0 0 0 1 1 1 0 0 1
26 1 0 1 1 0 1 0 0 1 1 0
27 1 1 1 1 0 1 0 1 1 1 0
28 1 0 1 0 1 1 1 0 1 0 0
29 1 0 1 1 1 1 1 1 1 0 0
30 1 1 1 1 1 1 1 1 1 1 1
31 1 0 0 0 0 0 0 0 0 0 0
"""


def small_basis() -> list:
    """TS 38.212 Table 5.3.3.3-1: row i holds M(i,0) to M(i,10)."""
    rows = []
    for line in SMALL_BASIS.splitlines():
        cells = [int(v) for v in line.split()]
        rows.append(cells[1:])
    return rows


def small_code_distance(k: int) -> int:
    m = small_basis()
    best = 32
    for a in range(1, 1 << k):
        w = sum(sum(((a >> n) & 1) * m[i][n] for n in range(k)) % 2 for i in range(32))
        best = min(best, w)
    return best


def check_nr() -> None:
    section("5G New Radio")
    same("number of lifting sizes", len(LIFTING_SIZES), 51)
    same("distinct lifting sizes", len(set(LIFTING_SIZES)), 51)
    same("smallest and largest", (LIFTING_SIZES[0], LIFTING_SIZES[-1]), (2, 384))
    check("BG1 lowest rate 22/66", 22 / 66, "0.333")
    check("BG2 lowest rate 10/50", 10 / 50, "0.2")
    same("BG1 largest code block, 22 x 384", 22 * 384, 8_448)
    same("BG2 largest code block", KCB[2], 3_840)
    a, r = 984, 0.5
    same("CRC bits for A = 984", crc_length(a), 16)
    b = a + crc_length(a)
    same("B for A = 984", b, 1000)
    same("base graph for A = 984, R = 1/2", base_graph(a, r), 2)
    same("code blocks", code_blocks(b, 2), 1)
    bg2 = lifting(b, 2)
    same("BG2: Kb, Zc, K, F, N", (bg2["Kb"], bg2["Zc"], bg2["K"], bg2["F"], bg2["N"]), (10, 104, 1040, 40, 5200))
    same("BG2: the lifting size below 104 is 96", max(z for z in LIFTING_SIZES if z < 104), 96)
    starts = tuple(k0(rv, 2, bg2["N"], bg2["Zc"]) for rv in range(4))
    same("BG2: k0 for rv 0 to 3", starts, (0, 1352, 2600, 4472))
    same("BG2: k0 in units of Zc", tuple(s // 104 for s in starts), (0, 13, 25, 43))
    bg1 = lifting(b, 1)
    same("BG1 instead: Zc, K, N", (bg1["Zc"], bg1["K"], bg1["N"]), (48, 1056, 3168))
    same("base graph for A = 984, R = 0.8", base_graph(984, 0.8), 1)
    same("base graph for A = 200, R = 0.8", base_graph(200, 0.8), 2)
    ex = lifting(500, 2)
    same("exercise B = 500, BG2: Kb, Zc, K, F", (ex["Kb"], ex["Zc"], ex["K"], ex["F"]), (8, 64, 640, 140))
    same("exercise: lifting sizes from 63 up", [z for z in LIFTING_SIZES if 60 <= z <= 64], [60, 64])
    same("segmentation of B = 20,000 with BG1", code_blocks(20_000, 1), 3)
    same("polar mother lengths 2^9 and 2^10", (2 ** 9, 2 ** 10), (512, 1024))
    m = small_basis()
    same("(32,K) basis: M(i,0) is all ones", all(row[0] == 1 for row in m), True)
    same("(32,K) basis: M(i,1) to M(i,5) take all 32 patterns", len({tuple(row[1:6]) for row in m}), 32)
    distances = {k: small_code_distance(k) for k in range(3, 12)}
    same("(32,K) minimum distance for K = 3 to 11", tuple(distances.values()), (16, 16, 16, 16, 12, 12, 12, 12, 10))
    k6 = {}
    for a6 in range(64):
        w = sum(sum(((a6 >> n) & 1) * m[i][n] for n in range(6)) % 2 for i in range(32))
        k6[w] = k6.get(w, 0) + 1
    same("(32,6) weight distribution", dict(sorted(k6.items())), {0: 1, 16: 62, 32: 1})


# The case summary table -------------------------------------------------------------------

def check_summary() -> None:
    section("The systems at a glance")
    check("byte parity", 8 / 9, "0.889")
    check("Hamming (7,4)", 4 / 7, "0.571")
    check("Mariner (32,6)", 6 / 32, "0.1875")
    check("Voyager concatenated", 223 / 510, "0.437")
    check("compact disc CIRC", 24 / 32, "0.75")
    check("QR version 1-M", 16 / 26, "0.615")
    check("Hsiao (72,64)", 64 / 72, "0.889")


def impulse_response(g1: int, g2: int, memory: int) -> str:
    """The output pairs when a single 1 enters an empty register and zeros follow."""
    out, state = [], 0
    for bit in [1] + [0] * memory:
        reg = (bit << memory) | state
        out.append(f"{weight(reg & g1) % 2}{weight(reg & g2) % 2}")
        state = reg >> 1
    return "".join(out)


def check_prose_extras() -> None:
    section("Further numbers in the prose, examples and exercises")
    # Hamming's Code at Bell Labs and Memory Codes after Hsiao
    same("data bits that six check bits protect, 2^6 - 1 - 6", 2 ** 6 - 1 - 6, 57)
    same("six check bits are too few for 64 data bits", 2 ** 6 >= 64 + 6 + 1, False)
    # Voyager
    response = impulse_response(G1, G2, 6)
    same("impulse response of the (171, 133) code", response, "11101111000111")
    same("weight of the impulse response", response.count("1"), 10)
    check("Viterbi output bit error rate quoted by Ludwig and Taylor", 5e-3, "5e-3")
    check("Reed-Solomon output bit error rate quoted by Ludwig and Taylor", 1e-6, "1e-6")
    # The Compact Disc
    same("frames spanned by one C2 codeword, 4 x 27 + 1", 4 * 27 + 1, 109)
    # QR Codes
    same("version 1 function modules, 441 - 208", 21 * 21 - 208, 233)
    same("finder 147 + separators 45 + timing 10 + format 30 + dark module 1", 147 + 45 + 10 + 30 + 1, 233)
    # ISO/IEC 18004 holds p check codewords back from correction in version 1:
    # p = 3 at level L, 2 at M, 1 at Q and 1 at H, so t = (r - p) // 2.
    reserve = {"L": 3, "M": 2, "Q": 1, "H": 1}
    specified = {lv: (26 - QR_DATA_CODEWORDS[lv] - reserve[lv]) // 2 for lv in "LMQH"}
    same("version 1 specified correction, (r - p) // 2 for L M Q H", tuple(specified.values()), (2, 4, 6, 8))
    same("level L, (7 - 3) // 2", (7 - 3) // 2, 2)
    same("level M, (10 - 2) // 2", (10 - 2) // 2, 4)
    check("version 1-L specified correction 2/26", specified["L"] / 26, "7.7%")
    check("version 1-M specified correction 4/26", specified["M"] / 26, "15.4%")
    check("version 1-Q specified correction 6/26, below the nominal 25%", specified["Q"] / 26, "23.1%")
    # Worked example: byte parity against SEC-DED
    p = 1e-6
    _, silent = byte_parity_word(p)
    check("silent parity error over silent SEC-DED error", silent / tail_at_least(72, 3, p), "4.8e3")
    # Worked example: the Mariner code against five-fold repetition
    p = 0.05
    check("six bits sent bare at p = 0.05", 1 - 0.95 ** 6, "0.2649")
    check("five-fold repetition, one bit", tail_at_least(5, 3, p), "0.00116")
    near("five-fold repetition, one bit, by the three-term formula",
         10 * p ** 3 * (1 - p) ** 2 + 5 * p ** 4 * (1 - p) + p ** 5, tail_at_least(5, 3, p), 1e-12)
    check("five-fold repetition, six bits", repetition5_word_error(p), "0.00693")
    check("Mariner code, eight or more errors in 32", rm15_word_error(p), "0.000139")
    check("repetition over Mariner", repetition5_word_error(p) / rm15_word_error(p), "50")
    same("extra bits the Mariner code spends over repetition", 32 - 30, 2)
    # Worked example: a scratch on a compact disc
    frame_mm = 1200 / 7350
    check("frame length on the track at 1.2 m/s, mm", frame_mm, "0.163")
    check("frame lengths in a 1.5 mm speck", 1.5 / frame_mm, "9.19")
    same("frames a 9.19-frame stretch can touch, ceil + 1", ceil(1.5 / frame_mm) + 1, 11)
    same("bytes of one C2 codeword lost to 11 frames", cd_max_hits(11), 3)
    check("longest fillable dropout at 1.2 m/s, mm", 16 * frame_mm, "2.61")
    # Exercise: decoding a Hamming word
    received = [1, 1, 1, 0, 0, 1, 1]
    s = hamming_syndrome(received)
    same("syndrome of 1110011", format(s, "03b"), "001")
    corrected = received[:]
    corrected[s - 1] ^= 1
    same("corrected word", "".join(map(str, corrected)), "0110011")
    same("data bits at positions 3, 5, 6, 7", "".join(str(corrected[i - 1]) for i in (3, 5, 6, 7)), "1011")
    # Exercise: silent errors in a nine-bit group
    p = 1e-3
    exact = even_count_probability(9, p) - (1 - p) ** 9
    check("exact silent probability, p = 1e-3", exact, "3.57e-5")
    check("approximation 36 p^2", 36 * p ** 2, "3.6e-5")
    check("approximation high by", 36 * p ** 2 / exact - 1, "0.7%")
    check("(1 - p)^7", (1 - p) ** 7, "0.993")
    # Exercise: the Mariner code with erasures
    same("largest t with 2t + 5 < 16", max(t for t in range(8) if 2 * t + 5 < 16), 5)
    same("errors given up to 5 erasures", 7 - 5, 2)
    # Exercise: a lighter Reed-Solomon code
    same("RS(255,239) distance", 255 - 239 + 1, 17)
    same("RS(255,239) errors corrected", (255 - 239) // 2, 8)
    check("1/8!", 1 / factorial(8), "2.48e-5")
    check("1/8! over 1/16!", factorial(16) / factorial(8), "5e8")
    # Exercise: QR version 1 at levels L and H
    same("level L: check codewords, corrected, erasures, bytes", (26 - 19, (26 - 19) // 2, 26 - 19, (19 * 8 - 12) // 8), (7, 3, 7, 17))
    same("level H: check codewords, corrected, erasures, bytes", (26 - 9, (26 - 9) // 2, 26 - 9, (9 * 8 - 12) // 8), (17, 8, 17, 7))
    same("bytes given up and codewords gained from L to H", (17 - 7, 8 - 3), (10, 5))
    # Exercise: a 5G block of 500 bits
    check("500 / 8", 500 / 8, "62.5")
    block = lifting(500, 2)
    same("B = 500: Kb, Zc, K, F, N", (block["Kb"], block["Zc"], block["K"], block["F"], block["N"]), (8, 64, 640, 140, 3200))
    same("B = 500: start of rv2", k0(2, 2, block["N"], block["Zc"]), 1600)


def main() -> None:
    print(f"  {'quantity':<62} {'computed':<16} stated")
    check_summary()
    check_parity()
    check_hamming()
    check_mariner()
    check_voyager()
    check_voyager_table()
    check_cd()
    check_qr()
    check_hsiao()
    check_parity_against_secded()
    check_nr()
    check_prose_extras()
    print("\nAll chapter 10 numbers reproduced.")


if __name__ == "__main__":
    main()
