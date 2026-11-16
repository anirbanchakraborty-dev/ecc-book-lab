"""Software model of the (72,64) SEC-DED encoder and decoder of the project chapter.

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
This module is the single source of every check matrix the project builds: the
generator gen_secded72_pkg.py writes the SystemVerilog package from the column
lists below, the cocotb tests in test_secded72.py compare the hardware with the
functions below, and check_secded72.py checks the column lists against the ones
"SECDED and Hsiao Codes" prints, column for column. It uses the Python standard
library only.

Bit conventions match hdl/secded72/:

    column     an integer whose bit i is row i of the check matrix, the form
               python/ch20_secded_and_hsiao_codes.py uses
    codeword   bit j is codeword position j, in the order of the column list
    data word  bit i is data bit i, which sits at codeword position
               data_positions[i]
    syndrome   an integer whose bit i is syndrome row i

Four codes are built. HSIAO is the balanced Hsiao (72,64) code, hsiao_72_64() of
"SECDED and Hsiao Codes", and HAMMING the extended Hamming (72,64) code,
extended_hamming(72, 7) there. Both are the chapter's subject. HSIAO_DICT, the same
weight-3 and weight-5 columns taken in dictionary order, and HSIAO_39_32, the
(39,32) Hsiao code, exist for the chapter's exercises.

A decoder output is a triple (data, corrected, uncorrectable). Two forms of the
uncorrectable flag are modeled:

    nomatch   the syndrome is nonzero and equals no column, the default, which
              raises the flag on every received word the decoder cannot correct
    parity    the parity form: for a Hsiao code the syndrome is nonzero with an
              even number of ones, and for the extended Hamming code the seven
              Hamming bits are nonzero while the overall parity bit is 0
"""

from itertools import combinations
import random

HSIAO = "hsiao"
HAMMING = "hamming"
HSIAO_DICT = "hsiao_dictionary"
HSIAO_39_32 = "hsiao_39_32"
CODE_NAMES = (HSIAO, HAMMING, HSIAO_DICT, HSIAO_39_32)
FLAGS = ("nomatch", "parity")

RANDOM_WORDS = 1000          # seeded random data words in the bench sweep
RANDOM_SEED = 30             # the seed those words are drawn with


def col(rows) -> int:
    """The column with ones in the given rows."""
    v = 0
    for i in rows:
        v |= 1 << i
    return v


# ---------------------------------------------------------------- check matrices

def hsiao_72_64_columns() -> list[int]:
    """All 56 weight-3 columns of 8 rows, then 8 rotations of a weight-5 column, then the identity."""
    data = [col(c) for c in combinations(range(8), 3)]
    data += [col((i + j) % 8 for j in range(5)) for i in range(8)]
    return data + [1 << i for i in range(8)]


def hsiao_dictionary_columns() -> list[int]:
    """The 56 weight-3 columns, then the first 8 weight-5 columns in dictionary order."""
    data = [col(c) for c in combinations(range(8), 3)]
    data += [col(c) for c in list(combinations(range(8), 5))[:8]]
    return data + [1 << i for i in range(8)]


def extended_hamming_72_64_columns() -> list[int]:
    """Position j holds j in rows 0 to 6 and a 1 in the parity row 7."""
    return [j | (1 << 7) for j in range(72)]


def hsiao_39_32_columns() -> list[int]:
    """The 35 weight-3 columns of 7 rows less three, then the identity."""
    dropped = {col((0, 1, 2)), col((0, 3, 4)), col((1, 5, 6))}
    data = [col(c) for c in combinations(range(7), 3) if col(c) not in dropped]
    return data + [1 << i for i in range(7)]


class Code:
    """One SEC-DED code: its columns, the positions of its data and check bits, and its encoder."""

    def __init__(self, name: str, columns: list[int], r: int, check_positions: list[int]):
        self.name = name
        self.columns = columns
        self.n = len(columns)
        self.r = r
        self.check_positions = check_positions
        self.data_positions = [j for j in range(self.n) if j not in check_positions]
        self.k = len(self.data_positions)
        self.parity_row = r - 1 if name == HAMMING else None
        # row_masks[i] has bit j set when position j is in row i
        self.row_masks = [sum(1 << j for j, c in enumerate(columns) if (c >> i) & 1) for i in range(r)]
        self.position_of = {c: j for j, c in enumerate(columns)}
        self.check_masks = self._solve_check_masks()

    def _solve_check_masks(self) -> list[int]:
        """For each check position, the data bits whose XOR sets it so that every syndrome bit is 0.

        The check bits c satisfy H_C c = H_D d over GF(2), where H_C and H_D are the
        columns at the check and the data positions. Gauss-Jordan elimination on H_C
        gives its inverse, and row t of H_C^-1 H_D is the mask for check t.
        """
        r = self.r
        # augmented rows: (check-part bits, data-part mask) for each syndrome row
        rows = []
        for i in range(r):
            left = sum(1 << t for t, p in enumerate(self.check_positions) if (self.columns[p] >> i) & 1)
            right = sum(1 << b for b, p in enumerate(self.data_positions) if (self.columns[p] >> i) & 1)
            rows.append([left, right])
        for t in range(r):
            pivot = next(i for i in range(t, r) if (rows[i][0] >> t) & 1)
            rows[t], rows[pivot] = rows[pivot], rows[t]
            for i in range(r):
                if i != t and (rows[i][0] >> t) & 1:
                    rows[i][0] ^= rows[t][0]
                    rows[i][1] ^= rows[t][1]
        return [rows[t][1] for t in range(r)]

    # -------------------------------------------------------- encoder and decoder

    def encode(self, d: int) -> int:
        c = 0
        for b, p in enumerate(self.data_positions):
            c |= ((d >> b) & 1) << p
        for t, p in enumerate(self.check_positions):
            c |= ((d & self.check_masks[t]).bit_count() & 1) << p
        return c

    def syndrome(self, word: int) -> int:
        s = 0
        for i, mask in enumerate(self.row_masks):
            s |= ((word & mask).bit_count() & 1) << i
        return s

    def data_bits(self, word: int) -> int:
        d = 0
        for b, p in enumerate(self.data_positions):
            d |= ((word >> p) & 1) << b
        return d

    def uncorrectable(self, s: int, matched: bool, flag: str) -> bool:
        if flag == "nomatch":
            return s != 0 and not matched
        if self.parity_row is None:
            return s != 0 and s.bit_count() % 2 == 0
        hamming_bits = s & ((1 << self.parity_row) - 1)
        return hamming_bits != 0 and not (s >> self.parity_row) & 1

    def decode(self, received: int, flag: str = "nomatch") -> tuple:
        """Return (data, corrected, uncorrectable) for one received word."""
        s = self.syndrome(received)
        j = self.position_of.get(s)
        matched = j is not None
        if matched:
            received ^= 1 << j
        return (self.data_bits(received), int(matched), int(self.uncorrectable(s, matched, flag)))


def build(name: str) -> Code:
    if name == HSIAO:
        return Code(name, hsiao_72_64_columns(), 8, list(range(64, 72)))
    if name == HSIAO_DICT:
        return Code(name, hsiao_dictionary_columns(), 8, list(range(64, 72)))
    if name == HAMMING:
        return Code(name, extended_hamming_72_64_columns(), 8, [0] + [1 << i for i in range(7)])
    if name == HSIAO_39_32:
        return Code(name, hsiao_39_32_columns(), 7, list(range(32, 39)))
    raise ValueError(f"unknown code {name!r}")


# ---------------------------------------------------------------- the bench sweep

def envelope_patterns(n: int) -> list[int]:
    """Every error pattern of weight 0, 1 and 2 on n bits: 1 + n + n(n-1)/2 of them."""
    patterns = [0]
    patterns += [1 << i for i in range(n)]
    patterns += [(1 << i) | (1 << j) for i, j in combinations(range(n), 2)]
    return patterns


def bench_words(k: int) -> list[int]:
    """The all-zero word, the all-ones word, the k walking-one words and the seeded random words."""
    rng = random.Random(RANDOM_SEED)
    words = [0, (1 << k) - 1] + [1 << i for i in range(k)]
    words += [rng.getrandbits(k) for _ in range(RANDOM_WORDS)]
    return words


def spec(code: Code, d: int, e: int) -> tuple:
    """The guarantee inside the envelope, wt(e) <= 2."""
    w = e.bit_count()
    if w == 0:
        return (d, 0, 0)
    if w == 1:
        return (d, 1, 0)
    return (d ^ code.data_bits(e), 0, 1)
