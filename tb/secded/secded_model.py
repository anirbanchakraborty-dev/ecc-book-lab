"""Software model of the Hsiao (13,8) SEC-DED encoder and decoder, mutants included.

The chapter is "Verifying ECC Hardware" in Error Correction Coding. This module is
the reference the cocotb tests in test_secded.py compare the hardware with, and the
oracle formal/secded/run_formal.py checks each solver counterexample against. It
uses the Python standard library only, and it is written from the code's check
matrix, not from the SystemVerilog, so a slip in one does not hide in the other.

Bit conventions match hdl/secded/:

    codeword   bits 0 to 7 are the data bits d0 to d7, bits 8 to 12 the check
               bits c0 to c4
    syndrome   an integer whose bit i is syndrome row i

A decoder output is a triple (data, corrected, uncorrectable). The specification
inside the design envelope, wt(e) <= 2, is

    wt(e) = 0   (d, 0, 0)
    wt(e) = 1   (d, 1, 0)
    wt(e) = 2   (d ^ (e & 0xFF), 0, 1)

and a case is wrong when the decoder's triple differs from that. A wrong case is
silent when the data is wrong and uncorrectable is low, so nothing warns the
system that consumes the word.
"""

from itertools import combinations

K = 8                       # data bits
N = 13                      # codeword bits
R = N - K                   # check bits, the syndrome width
MUTANTS = (0, 1, 2, 3, 4)   # 0 is the correct decoder


def column(rows: str) -> int:
    """A check-matrix column written row 0 first, as the book prints it, to an integer."""
    return sum(1 << i for i, bit in enumerate(rows) if bit == "1")


# The Hsiao (13,8) check matrix of "SECDED and Hsiao Codes", one column per codeword bit.
DATA_COLUMNS = [column(c) for c in ("11010", "11001", "10110", "10101",
                                    "01110", "01101", "01011", "00111")]
CHECK_COLUMNS = [1 << i for i in range(R)]
COLUMNS = DATA_COLUMNS + CHECK_COLUMNS


def weight(x: int) -> int:
    """The number of 1 bits in x."""
    return bin(x).count("1")


def encode(d: int) -> int:
    """The 13-bit codeword for data word d: check bit c_i makes syndrome row i vanish."""
    check = 0
    for j in range(K):
        if d >> j & 1:
            check ^= DATA_COLUMNS[j]
    return d | check << K


def syndrome(r: int, mutant: int = 0) -> int:
    """H times the received word. Mutant 1 drops received bit 6 from row 1."""
    s = 0
    for j in range(N):
        if r >> j & 1:
            col = COLUMNS[j]
            if mutant == 1 and j == 6:
                col &= ~0b00010
            s ^= col
    return s


def decode(r: int, mutant: int = 0) -> tuple:
    """(data, corrected, uncorrectable) for received word r, with mutant 0 to 4 planted."""
    s = syndrome(r, mutant)
    hit = [s == col for col in COLUMNS]
    flip = [hit[j] for j in range(K)]
    if mutant == 3:                       # d5's column steers the flip to d4
        flip[4] = hit[4] or hit[5]
        flip[5] = False
    if mutant == 4:                       # d0's flip is gated by received bit 7
        flip[0] = hit[0] and not (r >> 7 & 1)
    data = (r & 0xFF) ^ sum(1 << j for j in range(K) if flip[j])
    corrected = int(any(hit))
    uncorrectable = int(s != 0 and not corrected)
    if mutant == 2:                       # raised only for a syndrome with two 1s
        uncorrectable = int(uncorrectable and weight(s) == 2)
    return data, corrected, uncorrectable


def spec(d: int, e: int) -> tuple:
    """The output the specification requires for data d and error e with wt(e) <= 2."""
    w = weight(e)
    if w == 0:
        return d, 0, 0
    if w == 1:
        return d, 1, 0
    if w == 2:
        return d ^ (e & 0xFF), 0, 1
    raise ValueError(f"error pattern {e:#x} has weight {w}, outside the envelope")


def envelope_patterns() -> list:
    """The 92 error patterns of weight at most 2: zero, the 13 singles, the 78 doubles."""
    return ([0] + [1 << i for i in range(N)]
            + [1 << i | 1 << j for i, j in combinations(range(N), 2)])


def check_case(d: int, e: int, mutant: int = 0) -> tuple:
    """(wrong, silent) for one case: the decoder's verdict against the specification."""
    got = decode(encode(d) ^ e, mutant)
    wrong = got != spec(d, e)
    silent = got[0] != d and not got[2]
    return wrong, silent


def count_wrong(mutant: int) -> tuple:
    """(wrong, silent) over all 2^8 data words times the 92 envelope patterns."""
    wrong = silent = 0
    patterns = envelope_patterns()
    for d in range(2**K):
        for e in patterns:
            w, s = check_case(d, e, mutant)
            wrong += w
            silent += s
    return wrong, silent
