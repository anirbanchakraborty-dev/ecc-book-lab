"""Reproduce every number in chapter 29 of Error Correction Coding, "Case Study: ECC in Processor Caches".

Needs Python 3.10 or later and nothing beyond the standard library, and reads
nothing from the network:

    python3 python/ch29_case_study_ecc_in_processor_caches.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Whole numbers, bit patterns and hexadecimal words are compared exactly, and
every fractional value is printed in one fixed format, exponent form with
four decimals or with the one decimal the chapter states, or a decimal or a
percentage to a stated number of places, so the comparison of the printed text
is exact too.

The open design the chapter studies is the SEC-DED generator of the OpenTitan
project, read at commit 4ced702ca207ad420fbae5226ac8685a17ebef7b. The files
read, each at that commit, are

    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/util/design/secded_gen.py
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/util/design/data/secded_cfg.hjson
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_39_32_enc.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_39_32_dec.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_72_64_enc.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_72_64_dec.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_inv_39_32_enc.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_inv_39_32_dec.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_inv_72_64_enc.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/ip/prim/rtl/prim_secded_inv_72_64_dec.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/vendor/lowrisc_ibex/doc/03_reference/security.rst
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/vendor/lowrisc_ibex/doc/03_reference/icache.rst
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/vendor/lowrisc_ibex/rtl/ibex_core.sv
    https://github.com/lowRISC/opentitan/blob/4ced702ca207ad420fbae5226ac8685a17ebef7b/hw/vendor/lowrisc_ibex/rtl/ibex_icache.sv

Nothing is fetched when this program runs. The list of generated codes and the
encoder masks below were transcribed by hand from those files, and every
property the chapter claims for them is recomputed here from the transcription.

A column of a check matrix with r rows is held as an integer whose bit i is the
entry in row i, and a bit pattern is printed with row 0 first. A 64-bit mask is
held as a pair of 32-bit halves, high half first, and printed as sixteen
hexadecimal digits, so that the Octave twin can hold and print it the same way.
The book's own Hsiao matrices are imported from ch20_secded_and_hsiao_codes.py,
so the comparison uses the matrices that chapter prints.
"""

import sys
from itertools import combinations
from math import exp

from ch20_secded_and_hsiao_codes import hsiao_39_32, hsiao_72_64

COMMIT = "4ced702ca207ad420fbae5226ac8685a17ebef7b"


# ---------------------------------------------------------------- helpers

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


def weight(v: int) -> int:
    return bin(v).count("1")


def bits(v: int, r: int) -> str:
    """A column or syndrome as a string, row 0 first."""
    return "".join("1" if (v >> i) & 1 else "0" for i in range(r))


def joined(values) -> str:
    return ",".join(str(x) for x in values)


def hexw(v: int, nbits: int) -> str:
    """A word as upper-case hexadecimal, one digit per four bits, rounded up."""
    return format(v, f"0{-(-nbits // 4)}X")


def e4(x: float) -> str:
    return f"{x:.4e}"


def levels(w: int) -> int:
    """Levels of two-input gates in a balanced tree over w signals."""
    d, reach = 0, 1
    while reach < w:
        reach *= 2
        d += 1
    return d


def row_loads(cols, r: int) -> list[int]:
    return [sum((c >> i) & 1 for c in cols) for i in range(r)]


def rows_of(v: int, r: int) -> str:
    return "(" + ",".join(str(i) for i in range(r) if (v >> i) & 1) + ")"


# ---------------------------------------------------------------- OpenTitan data
# The 22 entries of util/design/data/secded_cfg.hjson at the commit above,
# as (code_type, k, m), in the order the file lists them.

CONFIG = [
    ("hsiao", 16, 6), ("hsiao", 22, 6), ("hsiao", 32, 7), ("hsiao", 57, 7),
    ("hsiao", 64, 8), ("hsiao", 128, 9),
    ("hamming", 16, 6), ("hamming", 32, 7), ("hamming", 64, 8), ("hamming", 68, 8),
    ("hamming", 128, 9),
    ("inv_hsiao", 16, 6), ("inv_hsiao", 22, 6), ("inv_hsiao", 32, 7), ("inv_hsiao", 57, 7),
    ("inv_hsiao", 64, 8), ("inv_hsiao", 128, 9),
    ("inv_hamming", 16, 6), ("inv_hamming", 32, 7), ("inv_hamming", 64, 8),
    ("inv_hamming", 68, 8), ("inv_hamming", 128, 9),
]

# The module-name infix secded_gen.py gives each code type (CODE_OPTIONS).
INFIX = {"hsiao": "", "inv_hsiao": "_inv", "hamming": "_hamming", "inv_hamming": "_inv_hamming"}

# Encoder masks of hw/ip/prim/rtl/prim_secded_39_32_enc.sv at the commit above:
# check bit 32 + i is the XOR of the data bits that mask i selects. The
# inverted encoder prim_secded_inv_39_32_enc.sv carries the same seven masks.
MASKS_39 = [0x2606BD25, 0xDEBA8050, 0x413D89AA, 0x31234ED1, 0xC2C1323B, 0x2DCC624C, 0x98505586]

# Encoder masks of hw/ip/prim/rtl/prim_secded_72_64_enc.sv at the commit above,
# as (high 32 bits, low 32 bits): check bit 64 + i is the XOR of the data bits
# that mask i selects. prim_secded_inv_72_64_enc.sv carries the same eight.
MASKS_72 = [
    (0xB9000000, 0x001FFFFF), (0x5E00000F, 0xFFE0003F), (0x67003FF0, 0x03E007C1),
    (0xCD0FC0F0, 0x3C207842), (0xB671C711, 0xC4438884), (0xB5B65926, 0x488C9108),
    (0xCBDAAA4A, 0x91152210), (0x7AED348D, 0x221A4420),
]

# The constants the inverted decoders XOR into the received word before the
# syndrome, from prim_secded_inv_39_32_dec.sv (39'h2A00000000) and
# prim_secded_inv_72_64_dec.sv (72'hAA0000000000000000), as check-bit patterns.
INV_CHECK_39 = 0x2A
INV_CHECK_72 = 0xAA


def mask_value(m) -> int:
    return (m[0] << 32) | m[1] if isinstance(m, tuple) else m


def ot_columns(masks, k: int) -> list[int]:
    """Data column j has a one in row i when bit j of mask i is set."""
    r = len(masks)
    full = [mask_value(m) for m in masks]
    data = [sum(((full[i] >> j) & 1) << i for i in range(r)) for j in range(k)]
    return data + [1 << i for i in range(r)]


def ideal_fanin(k: int, m: int) -> int:
    """The largest data fan-in of a check bit that balanced odd columns allow,
    computed as ideal_fanin() in secded_gen.py computes it."""
    fanin, needed = 0, k
    for w in range(3, m + 1, 2):
        avail = len(list(combinations(range(m), w)))
        take = min(avail, needed)
        fanin += -(-(take * w) // m)
        needed -= take
        if needed == 0:
            break
    return fanin


def distance_at_least_4(cols) -> bool:
    """No zero column, no repeated column, no column the sum of two others."""
    if any(c == 0 for c in cols) or len(set(cols)) != len(cols):
        return False
    present = set(cols)
    return all((a ^ b) not in present for a, b in combinations(cols, 2))


def weight4_codewords(cols) -> int:
    """Codewords of weight 4: pairs of column pairs with the same sum, over 3."""
    count = {}
    for a, b in combinations(cols, 2):
        x = a ^ b
        count[x] = count.get(x, 0) + 1
    return sum(c * (c - 1) // 2 for c in count.values()) // 3


def intersections(cols) -> str:
    """How many rows each pair of the given columns shares, as size:count."""
    tally = {}
    for a, b in combinations(cols, 2):
        s = weight(a & b)
        tally[s] = tally.get(s, 0) + 1
    return " ".join(f"{s}:{tally[s]}" for s in sorted(tally))


# ---------------------------------------------------------------- checks

def check_config() -> None:
    section("The codes secded_cfg.hjson asks the generator for, at commit " + COMMIT[:12])
    for kind, k, m in CONFIG:
        n = k + m
        name = f"prim_secded{INFIX[kind]}_{n}_{k}"
        print(f"  {name:<34} {kind:<12} ({n},{k})")
    same("entries in the configuration", len(CONFIG), 22)
    same("Hsiao codes", sum(1 for c in CONFIG if c[0] == "hsiao"), 6)
    same("extended Hamming codes", sum(1 for c in CONFIG if c[0] == "hamming"), 5)
    same("inverted Hsiao codes", sum(1 for c in CONFIG if c[0] == "inv_hsiao"), 6)
    same("inverted extended Hamming codes", sum(1 for c in CONFIG if c[0] == "inv_hamming"), 5)
    same("Hsiao lengths n", joined(k + m for kind, k, m in CONFIG if kind == "hsiao"), "22,28,39,64,72,137")
    same("extended Hamming lengths n", joined(k + m for kind, k, m in CONFIG if kind == "hamming"), "22,39,72,76,137")
    fewest = all(2 ** (m - 1) >= k + m and 2 ** (m - 2) < k + m - 1 for _, k, m in CONFIG)
    require("every entry uses the fewest check bits the SEC-DED bound allows", fewest)
    same("check bits of the (64,57) code", [m for kind, k, m in CONFIG if kind == "hsiao" and k == 57][0], 7)
    same("odd-weight 7-bit columns there are", sum(1 for v in range(1, 128) if weight(v) % 2), 64)
    same("so the (64,57) code uses every one of them", 57 + 7, 64)


def check_matrix(label: str, cols, r: int, k: int, stated: dict) -> None:
    section(f"The OpenTitan Hsiao {label} check matrix, from the transcribed masks")
    data = cols[:k]
    require("every column has odd weight", all(weight(c) % 2 == 1 for c in cols))
    require("the columns are distinct", len(set(cols)) == len(cols))
    require("no column is the sum of two others, so d >= 4", distance_at_least_4(cols))
    same("data column weights present", joined(sorted(set(weight(c) for c in data))), stated["wts"])
    loads = row_loads(cols, r)
    same("ones in each row", joined(loads), stated["rows"])
    same("ones in the whole matrix", sum(loads), stated["ones"])
    same("XOR gates for the syndrome, one fewer than each row", sum(w - 1 for w in loads), stated["syn"])
    same("XOR gates for the encoder, data ones less one per row", sum(w - 2 for w in loads), stated["enc"])
    same("levels of two-input gates for the longest row", levels(max(loads)), stated["levels"])
    same("ideal largest data fan-in, secded_gen.py's rule", ideal_fanin(k, r), stated["fanin"])
    same("largest data fan-in of the matrix", max(row_loads(data, r)), stated["fanin"])


def check_comparison(label: str, ot, book, r: int, k: int, stated: dict) -> None:
    section(f"OpenTitan's {label} matrix against the book's, from SECDED and Hsiao Codes")
    od, bd = ot[:k], book[:k]
    same("data columns the two matrices share", len(set(od) & set(bd)), stated["common"])
    same("data positions holding the same column", sum(1 for j in range(k) if od[j] == bd[j]), stated["agree"])
    same("columns only OpenTitan uses", " ".join(rows_of(c, r) for c in od if c not in bd), stated["ot_only"])
    same("columns only the book uses", " ".join(rows_of(c, r) for c in bd if c not in od), stated["book_only"])
    same("book row weights", joined(row_loads(book, r)), stated["book_rows"])
    same("book ones in the whole matrix", sum(row_loads(book, r)), stated["book_ones"])
    rows_differ = sorted(row_loads(ot, r)) != sorted(row_loads(book, r))
    same("row-weight lists differ after sorting", rows_differ, stated["rows_differ"])
    if "ot_w5" in stated:
        same("OpenTitan weight-5 columns, shared rows per pair",
             intersections([c for c in od if weight(c) == 5]), stated["ot_w5"])
        same("book weight-5 columns, shared rows per pair",
             intersections([c for c in bd if weight(c) == 5]), stated["book_w5"])
    same("weight-4 codewords, OpenTitan", weight4_codewords(ot), stated["ot_a4"])
    same("weight-4 codewords, book", weight4_codewords(book), stated["book_a4"])


def syndrome_class(s: int) -> str:
    if s == 0:
        return "none"
    return "single" if weight(s) % 2 else "double"


def check_inverted(label: str, cols, r: int, k: int, inv_check: int, stated: dict) -> None:
    section(f"The inverted {label} code, prim_secded_inv")
    n = k + r
    rule = sum(1 << (k + idx) for idx in range(1, r, 2))
    const = inv_check << k
    same("constant from the rule in secded_gen.py, hexadecimal", hexw(rule, n), stated["hex"])
    same("constant copied from the inverted decoder, hexadecimal", hexw(const, n), stated["hex"])
    same("check part of the constant, row 0 first", bits(inv_check, r), stated["check"])
    h_const = inv_check  # the constant is zero on the data, so its syndrome is its check part
    require("the constant is not a codeword, so the code is not linear", h_const != 0)
    ones = (1 << n) - 1
    s_ones = 0
    for c in cols:
        s_ones ^= c
    same("syndrome of the all-zero read", bits(h_const, r), stated["s_zero"])
    same("its weight", weight(h_const), stated["w_zero"])
    same("what the decoder flags for it", syndrome_class(h_const), stated["c_zero"])
    if "col_zero" in stated:
        same("data column it equals", cols.index(h_const), stated["col_zero"])
        same("data the inverted decoder returns, hexadecimal", hexw(1 << cols.index(h_const), k),
             stated["d_zero"])
    same("syndrome of the all-one word in the plain code", bits(s_ones, r), stated["s_ones_plain"])
    same("syndrome of the all-one read in the inverted code", bits(s_ones ^ h_const, r), stated["s_ones"])
    same("its weight", weight(s_ones ^ h_const), stated["w_ones"])
    same("what the decoder flags for it", syndrome_class(s_ones ^ h_const), stated["c_ones"])
    require("the all-zero word is not in the inverted code", h_const != 0)
    require("the all-one word is not in the inverted code", s_ones ^ h_const != 0)
    # the four words the generator's formal assertion is about
    data_mask = (1 << k) - 1

    def encode(d: int) -> int:
        check = 0
        for i in range(r):
            p = 0
            for j in range(k):
                if (d >> j) & 1 and (cols[j] >> i) & 1:
                    p ^= 1
            check |= p << i
        return d | ((check ^ inv_check) << k)

    words = [encode(0), encode(data_mask)]
    require("encoding all-zero or all-one data gives neither 0 nor all ones",
            all(w not in (0, ones) for w in words))
    same("check part stored with all-zero data", bits(words[0] >> k, r), stated["check"])
    same("check part stored with all-one data", bits(words[1] >> k, r), stated["enc_ones"])


def check_array() -> None:
    section("A 32 KiB data array, protected two ways")
    data_bits = 32 * 1024 * 8
    words = data_bits // 64
    same("data bytes", 32 * 1024, 32768)
    same("data bits", data_bits, 262144)
    same("64-bit words", words, 4096)
    same("extra bits for byte parity, one per byte", data_bits // 8, 32768)
    same("extra bits for (72,64) SEC-DED, eight per word", words * 8, 32768)
    same("overhead of either, percent", f"{100 * 32768 / data_bits:.1f}", "12.5")
    same("bits stored either way", data_bits + 32768, 294912)
    pairs = 72 * 71 // 2
    silent = 8 * (9 * 8 // 2)
    same("two-bit errors in a 72-bit stored word", pairs, 2556)
    same("of them inside one byte, so byte parity misses them", silent, 288)
    same("share byte parity misses, percent", f"{100 * silent / pairs:.2f}", "11.27")


def check_double_errors(cols) -> None:
    section("Every two-bit error in the OpenTitan (72,64) word")
    flagged = sum(1 for a, b in combinations(cols, 2) if syndrome_class(a ^ b) == "double")
    same("two-bit errors flagged as double", flagged, 2556)
    same("two-bit errors passing silently", 2556 - flagged, 0)


def encode72(cols, d: int) -> int:
    check = 0
    for j in range(64):
        if (d >> j) & 1:
            check ^= cols[j]
    return check


def check_partial_write(cols) -> None:
    section("A one-byte write into a (72,64) word, OpenTitan masks")
    d0 = 0x0123456789ABCDEF
    old = (d0 >> 16) & 0xFF
    new = 0x5A
    delta = old ^ new
    d1 = (d0 & ~(0xFF << 16)) | (new << 16)
    same("old data word", hexw(d0, 64), "0123456789ABCDEF")
    same("old byte 2", hexw(old, 8), "AB")
    same("new byte 2", hexw(new, 8), "5A")
    same("change in the byte", hexw(delta, 8), "F1")
    same("data bits that change", joined(16 + j for j in range(8) if (delta >> j) & 1), "16,20,21,22,23")
    same("new data word", hexw(d1, 64), "01234567895ACDEF")
    c0 = encode72(cols, d0)
    c1 = encode72(cols, d1)
    change = 0
    for j in range(8):
        if (delta >> j) & 1:
            change ^= cols[16 + j]
    same("old check bits, row 0 first", bits(c0, 8), "01101010")
    same("new check bits, re-encoded, row 0 first", bits(c1, 8), "00011111")
    same("sum of columns 16, 20, 21, 22 and 23", bits(change, 8), "01110101")
    require("old check bits plus that sum equal the re-encoded bits", c0 ^ change == c1)
    for j in (16, 20, 21, 22, 23):
        same(f"column {j}", bits(cols[j], 8), {16: "10001010", 20: "10000011", 21: "01110000",
                                                22: "01101000", 23: "01100100"}[j])
    # the same write when stored data bit 19 has flipped and the update skips the check
    stored = d0 ^ (1 << 19)
    old_bad = (stored >> 16) & 0xFF
    delta_bad = old_bad ^ new
    c_bad = c0
    for j in range(8):
        if (delta_bad >> j) & 1:
            c_bad ^= cols[16 + j]
    s = encode72(cols, d1) ^ c_bad
    same("byte 2 as read with bit 19 flipped", hexw(old_bad, 8), "A3")
    same("change computed from that byte", hexw(delta_bad, 8), "F9")
    same("check bits written by the unchecked update", bits(c_bad, 8), "10011010")
    same("syndrome of the next read", bits(s, 8), "10000101")
    same("column 19", bits(cols[19], 8), "10000101")
    same("data column the syndrome matches", cols.index(s), 19)
    same("data the decoder returns", hexw(d1 ^ (1 << cols.index(s)), 64), "012345678952CDEF")
    same("bit 19 of the data written", (d1 >> 19) & 1, 1)
    same("even-parity bit stored with byte 5A", weight(new) % 2, 0)


def at_least_two(mu: float) -> float:
    """P(two or more upsets) for a Poisson count of mean mu, summed from j = 2."""
    term = mu * mu / 2
    total = term
    for j in range(3, 31):
        term = term * mu / j
        total = total + term
    return exp(-mu) * total


def at_least_three(mu: float) -> float:
    """P(three or more upsets) for a Poisson count of mean mu, summed from j = 3."""
    term = mu * mu * mu / 6
    total = term
    for j in range(4, 31):
        term = term * mu / j
        total = total + term
    return exp(-mu) * total


def check_scrubbing() -> None:
    section("Scrubbing a 4,096-word array, illustrative upset rate")
    lam = 1e-6
    words = 4096
    year = 8760
    same("upset rate per word per hour", e4(lam), "1.0000e-06")
    same("the same rate in FIT per word", round(lam * 1e9), 1000)
    rows = []
    for t, f_stated, lost_stated in ((1, "5.0000e-13", "1.7940e-05"), (24, "2.8800e-10", "4.3056e-04"),
                                     (168, "1.4110e-08", "3.0137e-03"), (720, "2.5908e-07", "1.2911e-02")):
        mu = lam * t
        f = at_least_two(mu)
        lost = words * (year / t) * f
        rows.append((t, mu, f, lost))
        same(f"P(two or more in one word), T = {t} h", e4(f), f_stated)
        same(f"words lost a year, T = {t} h", e4(lost), lost_stated)
        require(f"mu^2/2 - mu^3/3 <= P <= mu^2/2, T = {t} h", mu * mu / 2 - mu ** 3 / 3 <= f <= mu * mu / 2)
        require(f"mu^2/2 overstates P by at most 2mu/3 of itself, T = {t} h",
                mu * mu / 2 - f <= (2 * mu / 3) * (mu * mu / 2))
    same("lambda T for a daily scrub", e4(rows[1][1]), "2.4000e-05")
    same("square law's error bound 2 lambda T / 3, daily scrub", f"{2 * rows[1][1] / 3:.1e}", "1.6e-05")
    same("linear estimate W Y lambda^2 T / 2, T = 24 h", e4(words * year * lam * lam * 24 / 2), "4.3057e-04")
    for t, share_stated in ((24, "8.0e-06"), (720, "2.4e-04"), (year, "2.9e-03")):
        mu = lam * t
        same(f"share of lost words with three or more upsets, T = {t} h",
             f"{at_least_three(mu) / at_least_two(mu):.1e}", share_stated)
    mu_year = lam * year
    lost_year = words * at_least_two(mu_year)
    same("lambda T with no scrub for a year", e4(mu_year), "8.7600e-03")
    same("words holding two or more upsets after that year", e4(lost_year), "1.5624e-01")
    same("the fraction of the 4,096 words that is", e4(at_least_two(mu_year)), "3.8145e-05")
    same("ratio of a year unscrubbed to a daily scrub", f"{lost_year / rows[1][3]:.0f}", "363")
    same("factor 1 - 2 lambda T / 3 with no scrub for a year", f"{1 - 2 * mu_year / 3:.3f}", "0.994")


def check_interleaving() -> None:
    section("Bursts in a 288-cell row holding four 72-bit words")
    cells = 288

    def worst(owner, start: int, length: int) -> int:
        hits = {}
        for c in range(start, start + length):
            w = owner(c)
            hits[w] = hits.get(w, 0) + 1
        return max(hits.values())

    def interleaved(c: int) -> int:
        return c % 4

    def side_by_side(c: int) -> int:
        return c // 72

    for lengths, stated_count, stated_worst in (((1, 2, 3, 4), 1146, 1), ((5, 6, 7, 8), 1130, 2)):
        count, most = 0, 0
        for length in lengths:
            for start in range(cells - length + 1):
                count += 1
                most = max(most, worst(interleaved, start, length))
        span = f"{lengths[0]} to {lengths[-1]}"
        same(f"bursts of {span} cells, interleaved", count, stated_count)
        same(f"most flips one word takes from such a burst", most, stated_worst)
    same("word and bit held by cell 9", f"{9 % 4},{9 // 4}", "1,2")
    bad = 0
    for length in (2, 3, 4):
        for start in range(cells - length + 1):
            if worst(side_by_side, start, length) >= 2:
                bad += 1
    same("bursts of 2 to 4 cells, words side by side, uncorrectable", bad, 855)
    same("bursts of 2 cells across a word boundary", sum(
        1 for s in range(cells - 1) if side_by_side(s) != side_by_side(s + 1)), 3)


def main() -> None:
    o39 = ot_columns(MASKS_39, 32)
    o72 = ot_columns(MASKS_72, 64)
    section("The transcribed encoder masks, hexadecimal")
    same("(39,32) masks", " ".join(hexw(m, 32) for m in MASKS_39),
         "2606BD25 DEBA8050 413D89AA 31234ED1 C2C1323B 2DCC624C 98505586")
    for i, m in enumerate(MASKS_72):
        print(f"  (72,64) mask {i}: {m[0]:08X}{m[1]:08X}")
    check_config()
    check_matrix("(39,32)", o39, 7, 32, {
        "wts": "3", "rows": "15,15,15,15,15,15,13", "ones": 103, "syn": 96, "enc": 89,
        "levels": 4, "fanin": 14})
    check_matrix("(72,64)", o72, 8, 64, {
        "wts": "3,5", "rows": "27,27,27,27,27,27,27,27", "ones": 216, "syn": 208, "enc": 200,
        "levels": 5, "fanin": 26})
    check_comparison("(39,32)", o39, hsiao_39_32(), 7, 32, {
        "common": 29, "agree": 1, "ot_only": "(0,3,4) (0,1,2) (1,5,6)",
        "book_only": "(0,1,6) (2,5,6) (3,4,6)", "book_rows": "14,14,15,15,15,15,15",
        "book_ones": 103, "rows_differ": True, "ot_a4": 1366, "book_a4": 1363})
    check_comparison("(72,64)", o72, hsiao_72_64(), 8, 64, {
        "common": 57, "agree": 56,
        "ot_only": "(0,2,3,5,6) (1,2,4,6,7) (0,1,3,6,7) (0,1,4,5,7) (0,2,4,5,7) (1,2,3,6,7) (0,3,4,5,6)",
        "book_only": "(0,1,2,3,4) (2,3,4,5,6) (3,4,5,6,7) (0,4,5,6,7) (0,1,5,6,7) (0,1,2,6,7) (0,1,2,3,7)",
        "book_rows": "27,27,27,27,27,27,27,27", "book_ones": 216, "rows_differ": False,
        "ot_w5": "2:8 3:16 4:4", "book_w5": "2:12 3:8 4:8", "ot_a4": 8396, "book_a4": 8392})
    check_inverted("(39,32)", o39, 7, 32, INV_CHECK_39, {
        "hex": "2A00000000", "check": "0101010", "s_zero": "0101010", "w_zero": 3,
        "c_zero": "single", "col_zero": 6, "d_zero": "00000040", "s_ones_plain": "1111111",
        "s_ones": "1010101", "w_ones": 4, "c_ones": "double", "enc_ones": "0101010"})
    check_inverted("(72,64)", o72, 8, 64, INV_CHECK_72, {
        "hex": "AA0000000000000000", "check": "01010101", "s_zero": "01010101", "w_zero": 4,
        "c_zero": "double", "s_ones_plain": "11111111", "s_ones": "10101010", "w_ones": 4,
        "c_ones": "double", "enc_ones": "01010101"})
    check_array()
    check_double_errors(o72)
    check_partial_write(o72)
    check_scrubbing()
    check_interleaving()
    print("\nAll chapter 29 numbers reproduced.")


if __name__ == "__main__":
    main()
