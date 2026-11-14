"""Reproduce every number in chapter 1 of Error Correction Coding, "Why Codes Exist".

Needs Python 3.9 or later and nothing beyond the standard library:

    python3 python/ch01_why_codes_exist.py

Each check prints the value the chapter quotes and raises AssertionError if the
computation disagrees with it, so a silent run past the last line means every
number in the chapter was reproduced.
"""

from fractions import Fraction
from math import ceil, comb, exp, expm1, isclose, log, log1p, log2


# One flipped bit ------------------------------------------------------------

def flip(word: int, bit: int) -> int:
    """Return word with bit number `bit` inverted (bit 0 is the least significant)."""
    return word ^ (1 << bit)


def as_signed8(u: int) -> int:
    """Read an 8-bit pattern as a two's-complement integer."""
    return u - 256 if u >= 128 else u


# How often errors happen -----------------------------------------------------

def p_at_least_one(p: float, n: float) -> float:
    """Probability that n independent bits, each wrong with probability p, hold one error or more.

    Computed as -expm1(n * log1p(-p)) so that p = 1e-12 does not round to zero.
    """
    return -expm1(n * log1p(-p))


# Repetition codes ------------------------------------------------------------

def rep_encode(bits: str, n: int = 3) -> str:
    return " ".join(b * n for b in bits)


def rep_decode_block(block: str) -> str:
    """Majority vote over one block of a repetition code of odd length."""
    return "1" if block.count("1") > len(block) // 2 else "0"


def rep_error(p, n: int):
    """Probability that majority decoding of R_n returns the wrong bit on a BSC(p).

    The vote is wrong when more than half of the n copies are flipped. Works on
    floats and on Fractions, so identities can be checked exactly.
    """
    assert n % 2 == 1
    return sum(comb(n, k) * p**k * (1 - p) ** (n - k) for k in range((n + 1) // 2, n + 1))


def rep_bound(p: float, n: int) -> float:
    """The bound (4p(1-p))^(n/2), proved in the chapter for p < 1/2."""
    return (4 * p * (1 - p)) ** (n / 2)


def majority(a: int, b: int, c: int) -> int:
    return (a & b) | (b & c) | (c & a)


# Information ---------------------------------------------------------------

def binary_entropy(p: float) -> float:
    if p in (0, 1):
        return 0.0
    return -p * log2(p) - (1 - p) * log2(1 - p)


# Single parity ---------------------------------------------------------------

def parity_undetected_direct(p, n: int):
    """Sum over every nonzero even error weight: those are the patterns parity cannot see."""
    return sum(comb(n, w) * p**w * (1 - p) ** (n - w) for w in range(2, n + 1, 2))


def parity_undetected_closed(p, n: int):
    """The closed form [1 + (1-2p)^n] / 2 - (1-p)^n derived in the chapter."""
    return (1 + (1 - 2 * p) ** n) / 2 - (1 - p) ** n


def xor_tree(n: int) -> tuple:
    """Build a balanced tree of two-input XOR gates over n inputs; return (gates, depth)."""
    level = [0] * n          # depth at which each signal becomes available
    gates = 0
    while len(level) > 1:
        nxt = []
        for i in range(0, len(level) - 1, 2):
            nxt.append(max(level[i], level[i + 1]) + 1)
            gates += 1
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    return gates, level[0]


def xor_bits(*words: str) -> str:
    width = len(words[0])
    acc = 0
    for w in words:
        acc ^= int(w, 2)
    return format(acc, f"0{width}b")


# Checks ----------------------------------------------------------------------

def show(label: str, value) -> None:
    print(f"  {label:<58} {value}")


def main() -> None:
    print("One flipped bit")
    x = 0b01100100
    assert x == 100
    hi, lo = flip(x, 7), flip(x, 0)
    assert format(hi, "08b") == "11100100" and hi == 228 and as_signed8(hi) == -28
    assert format(lo, "08b") == "01100101" and lo == 101
    show("01100100 is", x)
    show("flip bit 7: 11100100 unsigned, signed", (hi, as_signed8(hi)))
    show("flip bit 0: 01100101", lo)

    print("Where flips come from")
    p_dimm = 0.08
    p16 = 1 - (1 - p_dimm) ** 16
    assert isclose(p16, 0.7366, abs_tol=5e-5)
    show("P(one of 16 DIMMs errs in a year | 8% each, independent)", f"{p16:.4f}")

    print("How often")
    rate, ber = 8e9, 1e-12
    per_s = rate * ber
    assert isclose(per_s, 0.008) and isclose(per_s * 3600, 28.8)
    assert isclose(per_s * 86400, 691.2)
    show("errors per second at 8 Gb/s, BER 1e-12", per_s)
    show("errors per hour", round(per_s * 3600, 6))
    show("errors per day", round(per_s * 86400, 6))
    n_tb = 8e12                                    # one terabyte, 10^12 bytes
    exact = p_at_least_one(ber, n_tb)
    np_ = n_tb * ber
    assert isclose(np_, 8.0)
    assert isclose(exact, 1 - exp(-8), rel_tol=1e-9)
    assert isclose(exact, 0.99966, abs_tol=5e-6)
    show("terabyte: Np", np_)
    show("terabyte: P(at least one error)", f"{exact:.6f}")
    show("terabyte: 1 - e^-8", f"{1 - exp(-8):.6f}")
    n_pg = 4096 * 8
    small = p_at_least_one(ber, n_pg)
    assert n_pg == 32768
    assert isclose(small, n_pg * ber, rel_tol=1e-7)
    assert 1 - exp(-n_pg * ber) <= small <= n_pg * ber
    show("4 KiB page: Np", f"{n_pg * ber:.4e}")
    show("4 KiB page: P(at least one error)", f"{small:.6e}")
    # The lemma 1 - x <= e^-x and Bernoulli's inequality, spot-checked.
    for t in (0.0, 1e-9, 0.01, 0.5, 1.0, 3.0):
        assert 1 - t <= exp(-t)
    for p in (Fraction(1, 10), Fraction(1, 1000)):
        for n in range(0, 60):
            assert (1 - p) ** n >= 1 - n * p

    print("The binary symmetric channel")
    # Three bits at p = 0.1: the probability of each error weight, which the
    # chapter's table lists, and the check that the four add to one.
    p = Fraction(1, 10)
    one_pattern = {w: p ** w * (1 - p) ** (3 - w) for w in range(4)}
    by_weight = {w: comb(3, w) * one_pattern[w] for w in range(4)}
    assert [float(one_pattern[w]) for w in range(4)] == [0.729, 0.081, 0.009, 0.001]
    assert [float(by_weight[w]) for w in range(4)] == [0.729, 0.243, 0.027, 0.001]
    assert sum(by_weight.values()) == 1
    show("one pattern of weight 0..3, n = 3, p = 0.1", [float(one_pattern[w]) for w in range(4)])
    show("all patterns of weight 0..3", [float(by_weight[w]) for w in range(4)])

    print("The repetition code R3")
    msg = "1011"
    sent = rep_encode(msg)
    assert sent == "111 000 111 111"
    recv = "101 000 100 110"
    decoded = "".join(rep_decode_block(b) for b in recv.split())
    assert decoded == "1001"
    show("encode 1011", sent)
    show("decode 101 000 100 110", decoded)
    table = {format(w, "03b"): rep_decode_block(format(w, "03b")) for w in range(8)}
    assert [w for w, d in table.items() if d == "0"] == ["000", "001", "010", "100"]
    show("words decoded to 0", [w for w, d in table.items() if d == "0"])
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                assert majority(a, b, c) == (1 if a + b + c >= 2 else 0)
    show("maj(a,b,c) = ab + bc + ca matches the vote on all 8 rows", True)

    print("How much R3 buys")
    for p in (Fraction(1, 10), Fraction(1, 100), Fraction(1, 1000), Fraction(1, 3), Fraction(2, 5)):
        assert rep_error(p, 3) == 3 * p**2 - 2 * p**3
    quoted = {0.1: 0.028, 0.01: 2.98e-4, 0.001: 2.998e-6}
    for p, want in quoted.items():
        got = rep_error(p, 3)
        assert isclose(got, want, rel_tol=1e-9)
        show(f"P_e(R3) at p = {p}", f"{got:.4g}  (gain {p / got:.1f}x)")
    for p in (Fraction(1, 10), Fraction(49, 100), Fraction(1, 2), Fraction(51, 100), Fraction(9, 10)):
        better = rep_error(p, 3) < p
        assert better == (p < Fraction(1, 2))
    show("R3 beats uncoded exactly when p < 1/2", True)

    print("The price: rate")
    r5 = rep_error(0.01, 5)
    assert isclose(r5, 9.8506e-6, rel_tol=1e-4)
    show("P_e(R5) at p = 0.01", f"{r5:.5g}")
    for p in (Fraction(1, 100), Fraction(1, 10), Fraction(1, 3)):
        assert rep_error(p, 5) == 10 * p**3 - 15 * p**4 + 6 * p**5
    show("P_e(R5) = 10p^3 - 15p^4 + 6p^5, checked exactly", True)
    for p in (0.001, 0.01, 0.1, 0.3, 0.49):
        for n in range(1, 41, 2):
            assert rep_error(p, n) <= rep_bound(p, n) + 1e-15
    show("P_e(R_n) <= (4p(1-p))^(n/2) for odd n < 41, five values of p", True)
    show("(4p(1-p))^(1/2) at p = 0.1", round((4 * 0.1 * 0.9) ** 0.5, 12))
    h = binary_entropy(0.1)
    assert isclose(h, 0.4690, abs_tol=5e-5) and isclose(1 - h, 0.5310, abs_tol=5e-5)
    show("H(0.1)", f"{h:.4f}")
    show("capacity 1 - H(0.1)", f"{1 - h:.4f}")
    for n in (1, 3, 5, 7, 9, 11, 15, 25):
        show(f"R{n}: rate {1 / n:.4f}, P_e at p = 0.1", f"{rep_error(0.1, n):.3e}")

    print("Single parity")
    for p in (Fraction(1, 100), Fraction(1, 10), Fraction(1, 3)):
        for n in range(2, 17):
            assert parity_undetected_direct(p, n) == parity_undetected_closed(p, n)
    pud = parity_undetected_closed(0.01, 8)
    raw = 1 - 0.99**8
    assert isclose(pud, 0.0026368, abs_tol=5e-8)
    assert isclose(raw, 0.077255, abs_tol=5e-7)
    two = comb(8, 2) * 0.01**2 * 0.99**6
    assert isclose(two, 0.0026361, abs_tol=5e-8)
    show("P(undetected) for n = 8, p = 0.01", f"{pud:.7f}")
    show("of which weight-2 patterns", f"{two:.7f}")
    show("P(any error) in the same 8 bits", f"{raw:.6f}")
    show("fraction of erroneous words that slip through", f"{pud / raw:.4f}")
    for n in range(1, 129):
        g, d = xor_tree(n)
        assert g == n - 1 and d == (ceil(log2(n)) if n > 1 else 0)
    show("XOR tree, 8 inputs: (gates, depth)", xor_tree(8))
    show("XOR tree, 64 inputs: (gates, depth)", xor_tree(64))
    d1, d2, d3 = "1011", "0110", "1100"
    par = xor_bits(d1, d2, d3)
    assert par == "0001"
    assert xor_bits(d1, d3, par) == d2
    show("parity block P = D1 xor D2 xor D3", par)
    show("rebuilt D2 = D1 xor D3 xor P", xor_bits(d1, d3, par))
    word = "1011001"
    pbit = str(word.count("1") % 2)
    assert word + pbit == "10110010"
    show("even-parity codeword for 1011001", word + pbit)

    print("Detection, retransmission and correction")
    ok = (1 - 1e-4) ** 1000
    assert isclose(ok, 0.904833, abs_tol=5e-7)
    assert isclose(1 / ok, 1.105176, abs_tol=5e-6)
    show("P(1000-bit frame arrives clean at p = 1e-4)", f"{ok:.6f}")
    show("expected transmissions per frame", f"{1 / ok:.6f}")

    print("Worked examples")
    # The length at which an error becomes as likely as not, at BER 1e-12.
    ber = 1e-12
    n_half = log(2) / ber
    assert isclose(n_half, 6.931e11, rel_tol=1e-4)
    assert isclose(n_half / 8 / 1e9, 86.6, abs_tol=0.05)
    assert isclose(p_at_least_one(ber, n_half), 0.5, rel_tol=1e-9)
    show("bits until an error is as likely as not, BER 1e-12", f"{n_half:.4e}")
    show("the same length in gigabytes (10^9 bytes)", f"{n_half / 8 / 1e9:.1f}")
    # R3 on a channel with p = 0.05.
    p = Fraction(1, 20)
    pe = rep_error(p, 3)
    assert pe == Fraction(29, 4000) and float(pe) == 0.00725
    assert isclose(float(p / pe), 6.9, abs_tol=0.05)
    assert isclose(rep_bound(0.05, 3), 0.0828, abs_tol=5e-5)
    show("P_e(R3) at p = 0.05, and the gain", (float(pe), round(float(p / pe), 1)))
    show("the bound (4p(1-p))^(3/2) at p = 0.05", round(rep_bound(0.05, 3), 4))
    # A 64-bit word with one parity bit, so n = 65, at p = 1e-6.
    p = Fraction(1, 10**6)
    any65 = 1 - (1 - p) ** 65
    ud65 = parity_undetected_closed(p, 65)
    assert isclose(float(any65), 6.4998e-5, rel_tol=1e-4)
    assert isclose(float(ud65), 2.0799e-9, rel_tol=1e-4)
    assert isclose(float(comb(65, 2) * p**2), 2.08e-9, rel_tol=1e-12)
    assert isclose(float(ud65 / any65), 3.2e-5, rel_tol=1e-3)
    show("n = 65, p = 1e-6: P(any error)", f"{float(any65):.4e}")
    show("n = 65, p = 1e-6: P(undetected)", f"{float(ud65):.4e}")
    show("n = 65, p = 1e-6: C(65,2) p^2", f"{float(comb(65, 2) * p**2):.4e}")
    # A second RAID rebuild, this time losing the third disk.
    d1, d2, d3 = "1110", "0101", "1000"
    par = xor_bits(d1, d2, d3)
    assert par == "0011" and xor_bits(d1, d2, par) == d3
    show("P for 1110, 0101, 1000; rebuilt D3", (par, xor_bits(d1, d2, par)))
    # Decoding a received R3 sequence.
    recv = "110 011 001"
    assert "".join(rep_decode_block(b) for b in recv.split()) == "110"
    show("decode 110 011 001", "110")
    print("All chapter 1 numbers reproduced.")


if __name__ == "__main__":
    main()
