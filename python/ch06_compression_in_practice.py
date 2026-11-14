"""Reproduce every number in chapter 6 of Error Correction Coding, "Compression in Practice".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch06_compression_in_practice.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. Codewords, parses,
tables and decoded strings must match the chapter exactly. Every encoder here
has a decoder, and each pair is run as a round trip. A run that reaches the
final line has reproduced every number in the chapter.

The module also holds the arithmetic-coding, LZ78 and quantizer functions that
ch06_charts.py plots, so the charts and this check share one source.
"""

import random
import sys
from fractions import Fraction
from math import ceil, comb, erf, exp, lgamma, log, log2, log10, pi, sqrt


# Helpers ------------------------------------------------------------------------

def half_unit(stated: str) -> float:
    text = stated.replace(",", "")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", ""))
    print(f"  {label:<62} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<62} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<62} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH in '{label}'")
        sys.exit(1)


def entropy(probs) -> float:
    return sum(float(p) * log2(1 / float(p)) for p in probs if p > 0)


def binary_entropy(p: float) -> float:
    return 0.0 if p in (0, 1) else -p * log2(p) - (1 - p) * log2(1 - p)


# The running source: a 0.7, b 0.2, c 0.1 ---------------------------------------

LETTERS = "abc"
PROBS = {"a": Fraction(7, 10), "b": Fraction(2, 10), "c": Fraction(1, 10)}
FREQ = {"a": 7, "b": 2, "c": 1}          # integer counts for range coding and rANS
TOTAL = 10


def cumulative(probs: dict, order: str) -> dict:
    """Start of each letter's slice of [0, 1), letters taken in the given order."""
    out, run = {}, type(next(iter(probs.values())))(0)
    for s in order:
        out[s] = run
        run += probs[s]
    return out


# Arithmetic coding ---------------------------------------------------------------

def arith_intervals(message: str, probs: dict, order: str) -> list:
    """The interval (low, width) after each letter, starting from (0, 1), in exact fractions."""
    cum = cumulative(probs, order)
    low, width = Fraction(0), Fraction(1)
    steps = [(low, width)]
    for s in message:
        low, width = low + width * cum[s], width * probs[s]
        steps.append((low, width))
    return steps


def arith_codeword(low: Fraction, width: Fraction) -> str:
    """Shannon-Fano-Elias tag: l = ceil(log2(1/width)) + 1 bits of ceil(low * 2^l) / 2^l."""
    l0 = 0
    while Fraction(1, 2 ** l0) > width:     # smallest l0 with 2^-l0 <= width, which is ceil(log2(1/width))
        l0 += 1
    l = l0 + 1
    v = ceil(low * 2 ** l)
    return format(v, f"0{l}b")


def arith_encode(message: str, probs: dict, order: str) -> str:
    low, width = arith_intervals(message, probs, order)[-1]
    return arith_codeword(low, width)


def arith_decode(bits: str, n: int, probs: dict, order: str):
    """Decode n letters from a tag; returns the letters and the scaled value before each step."""
    cum = cumulative(probs, order)
    t = Fraction(int(bits, 2), 2 ** len(bits))
    out, values = [], []
    for _ in range(n):
        values.append(t)
        for s in order:
            if cum[s] <= t < cum[s] + probs[s]:
                out.append(s)
                t = (t - cum[s]) / probs[s]
                break
    return "".join(out), values


def arith_expected_rate(p: float, n: int) -> float:
    """Expected bits per letter of the arithmetic code on n letters of a binary source with P(1) = p."""
    q, total = 1 - p, 0.0
    for k in range(n + 1):
        log_prob = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + k * log(p) + (n - k) * log(q)
        info = -(k * log2(p) + (n - k) * log2(q))
        total += exp(log_prob) * (ceil(info - 1e-12) + 1)
    return total / n


def settle_bits(low: Fraction, high: Fraction):
    """Emit the leading bits shared by every number in [low, high), rescaling as they go."""
    out = ""
    while True:
        if high <= Fraction(1, 2):
            out += "0"
            low, high = 2 * low, 2 * high
        elif low >= Fraction(1, 2):
            out += "1"
            low, high = 2 * low - 1, 2 * high - 1
        else:
            return out, low, high


def laplace_probability(message: str, alphabet: str) -> Fraction:
    """Probability of a message under the adaptive model whose counts start at 1 for every letter."""
    counts = {s: 1 for s in alphabet}
    prob = Fraction(1)
    for s in message:
        prob *= Fraction(counts[s], sum(counts.values()))
        counts[s] += 1
    return prob


def laplace_closed_form(message: str, alphabet: str) -> Fraction:
    from math import factorial
    m, n = len(alphabet), len(message)
    num = factorial(m - 1)
    for s in alphabet:
        num *= factorial(message.count(s))
    return Fraction(num, factorial(n + m - 1))


# Range coding with integers ------------------------------------------------------

def range_step(low: int, rng: int, s: str) -> tuple:
    cum = cumulative({k: FREQ[k] for k in LETTERS}, LETTERS)
    return low + rng * cum[s] // TOTAL, rng * FREQ[s] // TOTAL


# rANS ------------------------------------------------------------------------------

CUM = cumulative({k: FREQ[k] for k in LETTERS}, LETTERS)


def rans_push(x: int, s: str, M: int = TOTAL) -> int:
    f, c = FREQ[s], CUM[s]
    return (x // f) * M + c + (x % f)


def rans_pop(x: int, M: int = TOTAL) -> tuple:
    slot = x % M
    for s in LETTERS:
        if CUM[s] <= slot < CUM[s] + FREQ[s]:
            return s, FREQ[s] * (x // M) + slot - CUM[s]
    raise ValueError("slot outside every letter")


def rans_encode_plain(message: str, x0: int) -> list:
    """States while encoding the message backward with no renormalization: [x0, x1, ...]."""
    states = [x0]
    for s in reversed(message):
        states.append(rans_push(states[-1], s))
    return states


def rans_decode_plain(x: int, n: int) -> tuple:
    out, states = [], [x]
    for _ in range(n):
        s, x = rans_pop(x)
        out.append(s)
        states.append(x)
    return "".join(out), states


def rans_encode_renorm(message: str, L: int = 10, b: int = 10) -> tuple:
    """rANS keeping x in [L, bL): before encoding s, emit digits while x >= (L // M) * b * f_s."""
    x, digits, log_ = L, [], []
    for s in reversed(message):
        emitted = []
        while x >= (L // TOTAL) * b * FREQ[s]:
            emitted.append(x % b)
            digits.append(x % b)
            x //= b
        new = rans_push(x, s)
        log_.append((s, x, emitted, new))
        x = new
    return x, digits, log_


def rans_decode_renorm(x: int, digits: list, n: int, L: int = 10, b: int = 10) -> tuple:
    stack, out = list(digits), []
    for _ in range(n):
        s, x = rans_pop(x)
        out.append(s)
        while x < L:
            x = x * b + stack.pop()
    return "".join(out), x


def kth_state(s: str, k: int, limit: int = 1000) -> int:
    """The k-th number, counting from 0, whose slot x mod M belongs to letter s."""
    seen = -1
    for x in range(limit):
        if CUM[s] <= x % TOTAL < CUM[s] + FREQ[s]:
            seen += 1
            if seen == k:
                return x
    raise ValueError


# LZ77 ------------------------------------------------------------------------------

def lz77_encode(text: str) -> list:
    """Greedy triples (offset, length, next letter) with an unbounded window; ties go to the nearest match."""
    i, out = 0, []
    while i < len(text):
        best_len, best_off = 0, 0
        for off in range(1, i + 1):
            length = 0
            while i + length < len(text) - 1 and text[i + length - off] == text[i + length]:
                length += 1
            if length > best_len:
                best_len, best_off = length, off
        out.append((best_off, best_len, text[i + best_len]))
        i += best_len + 1
    return out


def lz77_decode(triples: list) -> str:
    out = []
    for off, length, ch in triples:
        for _ in range(length):
            out.append(out[-off])
        out.append(ch)
    return "".join(out)


# LZ78 and LZW ----------------------------------------------------------------------

def lz78_encode(text: str) -> tuple:
    phrases, out, w = {"": 0}, [], ""
    for ch in text:
        if w + ch in phrases:
            w += ch
        else:
            out.append((phrases[w], ch))
            phrases[w + ch] = len(phrases)
            w = ""
    if w:
        out.append((phrases[w], ""))
    return out, phrases


def lz78_decode(pairs: list) -> str:
    table, out = [""], []
    for idx, ch in pairs:
        phrase = table[idx] + ch
        table.append(phrase)
        out.append(phrase)
    return "".join(out)


def lz78_bits(bits: list) -> float:
    """Bits per letter of LZ78 on a binary sequence: phrase i costs ceil(log2 i) + 1 bits."""
    trie, node, nxt, phrases, total = {}, 0, 1, 0, 0
    for bit in bits:
        if (node, bit) in trie:
            node = trie[(node, bit)]
        else:
            phrases += 1
            total += (ceil(log2(phrases)) if phrases > 1 else 0) + 1
            trie[(node, bit)] = nxt
            nxt += 1
            node = 0
    if node:
        phrases += 1
        total += ceil(log2(phrases))
    return total / len(bits)


def binary_source(p: float, n: int, seed: int = 2026) -> list:
    rng = random.Random(seed)
    return [1 if rng.random() < p else 0 for _ in range(n)]


def lzw_encode(text: str, alphabet: str) -> tuple:
    table = {ch: i for i, ch in enumerate(alphabet)}
    out, steps, w = [], [], ""
    for ch in text:
        if w + ch in table:
            w += ch
        else:
            out.append(table[w])
            table[w + ch] = len(table)
            steps.append((w, ch, table[w], w + ch, table[w + ch]))
            w = ch
    out.append(table[w])
    steps.append((w, "", table[w], "", None))
    return out, steps


def lzw_decode(codes: list, alphabet: str) -> tuple:
    table = {i: ch for i, ch in enumerate(alphabet)}
    prev = table[codes[0]]
    out, special = [prev], []
    for code in codes[1:]:
        if code in table:
            entry = table[code]
        else:
            entry = prev + prev[0]
            special.append(code)
        table[len(table)] = prev + entry[0]
        out.append(entry)
        prev = entry
    return "".join(out), special, table


# Run-length coding -----------------------------------------------------------------

def run_lengths(bits: str) -> list:
    runs, cur, count = [], "0", 0
    for b in bits:
        if b == cur:
            count += 1
        else:
            runs.append(count)
            cur, count = b, 1
    runs.append(count)
    return runs


def run_decode(runs: list) -> str:
    return "".join(("0" if i % 2 == 0 else "1") * r for i, r in enumerate(runs))


# Burrows-Wheeler and move-to-front -------------------------------------------------

def bwt(text: str) -> tuple:
    rows = sorted(text[i:] + text[:i] for i in range(len(text)))
    return "".join(r[-1] for r in rows), rows


def lf_map(last: str) -> list:
    first = sorted(last)
    start = {}
    for i, ch in enumerate(first):
        start.setdefault(ch, i)
    seen, out = {}, []
    for ch in last:
        out.append(start[ch] + seen.get(ch, 0))
        seen[ch] = seen.get(ch, 0) + 1
    return out


def inverse_bwt(last: str, end: str = "$") -> tuple:
    lf = lf_map(last)
    row, letters, path = 0, [], [0]
    while last[row] != end:
        letters.append(last[row])
        row = lf[row]
        path.append(row)
    return "".join(reversed(letters)) + end, path


def mtf_encode(text: str, alphabet: str) -> tuple:
    lst, out, lists = list(alphabet), [], []
    for ch in text:
        lists.append("".join(lst))
        i = lst.index(ch)
        out.append(i)
        lst.insert(0, lst.pop(i))
    return out, lists


def mtf_decode(indices: list, alphabet: str) -> str:
    lst, out = list(alphabet), []
    for i in indices:
        ch = lst.pop(i)
        out.append(ch)
        lst.insert(0, ch)
    return "".join(out)


def empirical_entropy(seq) -> float:
    n = len(seq)
    return entropy([Fraction(list(seq).count(s), n) for s in set(seq)])


# Quantization, transforms and rate-distortion ---------------------------------------

def quantize(x: float, step: float) -> tuple:
    k = round(x / step)
    return k, k * step


def phi(x: float) -> float:
    return exp(-x * x / 2) / sqrt(2 * pi) if abs(x) < 40 else 0.0


def Phi(x: float) -> float:
    return 0.5 * (1 + erf(x / sqrt(2)))


def cell_mse(a: float, b: float, y: float) -> float:
    """Integral of (x - y)^2 phi(x) over [a, b] for the standard Gaussian density phi."""
    afa = a * phi(a) if abs(a) < 40 else 0.0
    bfb = b * phi(b) if abs(b) < 40 else 0.0
    return (1 + y * y) * (Phi(b) - Phi(a)) - (bfb - afa) + 2 * y * (phi(b) - phi(a))


def fixed_rate_mse(levels: int, step: float) -> float:
    total = 0.0
    for i in range(levels):
        a = -50.0 if i == 0 else (i - levels / 2) * step
        b = 50.0 if i == levels - 1 else (i + 1 - levels / 2) * step
        total += cell_mse(a, b, (i - (levels - 1) / 2) * step)
    return total


def golden_min(f, lo: float, hi: float, iterations: int = 200) -> float:
    g = (sqrt(5) - 1) / 2
    a, b = lo, hi
    for _ in range(iterations):
        c, d = b - g * (b - a), a + g * (b - a)
        if f(c) < f(d):
            b = d
        else:
            a = c
    return (a + b) / 2


def best_uniform(levels: int) -> tuple:
    step = golden_min(lambda d: fixed_rate_mse(levels, d), 0.01, 3.0)
    return step, fixed_rate_mse(levels, step)


def entropy_coded_uniform(step: float) -> tuple:
    """Rate (index entropy) and mean squared error of a uniform quantizer on a unit Gaussian."""
    rate, mse, K = 0.0, 0.0, int(12 / step) + 2
    for k in range(-K, K + 1):
        a, b = (k - 0.5) * step, (k + 0.5) * step
        p = Phi(b) - Phi(a)
        if p > 0:
            rate -= p * log2(p)
        mse += cell_mse(a, b, k * step)
    return rate, mse


def db(x: float) -> float:
    return 10 * log10(x)


# Base-delta-immediate compression ----------------------------------------------------

def fits_signed(v: int, nbytes: int) -> bool:
    return -(1 << (8 * nbytes - 1)) <= v < (1 << (8 * nbytes - 1))


def bdi_compress(words: list, word_bytes: int, delta_bytes: int, immediate: bool):
    """One explicit base (the first word not small enough for the zero base) and signed deltas.
    Returns (base, deltas, mask, size in bytes) or None when some word does not fit."""
    base = None
    for w in words:
        if not (immediate and fits_signed(w, delta_bytes)):
            base = w
            break
    if base is None:
        base = words[0]
    deltas, mask = [], []
    for w in words:
        if immediate and fits_signed(w, delta_bytes):
            deltas.append(w)
            mask.append(0)
        elif fits_signed(w - base, delta_bytes):
            deltas.append(w - base)
            mask.append(1)
        else:
            return None
    size = word_bytes + delta_bytes * len(words) + (ceil(len(words) / 8) if immediate else 0)
    return base, deltas, mask, size


def bdi_decompress(base: int, deltas: list, mask: list) -> list:
    return [base + d if m else d for d, m in zip(deltas, mask)]


# The chapter's numbers, section by section ----------------------------------------------

def bits_of(prob: Fraction) -> float:
    """log2(1/prob) for an exact fraction too small for a float."""
    return log2(prob.denominator) - log2(prob.numerator)


def check_intervals() -> None:
    print("Coding with Intervals")
    check("entropy of the source a 0.7, b 0.2, c 0.1", entropy(PROBS.values()), "1.1568")
    huff = {"a": "0", "b": "10", "c": "11"}
    check("average length of Huffman's code", sum(float(PROBS[s]) * len(huff[s]) for s in LETTERS), "1.3")
    msg = "aabaaa"
    prob = PROBS["a"] ** 5 * PROBS["b"]
    same("probability of aabaaa", prob, Fraction("0.033614"))
    check("information in aabaaa", log2(1 / prob), "4.8948")
    same("bits Huffman's code spends on aabaaa", sum(len(huff[s]) for s in msg), 7)
    check("information in one a", log2(10 / 7), "0.5146")
    steps = arith_intervals(msg, PROBS, LETTERS)
    table = [("0", "0.7", "0.7"), ("0", "0.49", "0.49"), ("0.343", "0.441", "0.098"),
             ("0.343", "0.4116", "0.0686"), ("0.343", "0.39102", "0.04802"), ("0.343", "0.376614", "0.033614")]
    for k, (lo, hi, w) in enumerate(table, start=1):
        low, width = steps[k]
        same(f"interval and width after step {k}", (low, low + width, width),
             (Fraction(lo), Fraction(hi), Fraction(w)))
    same("b keeps 70 to 90 percent of [0, 0.49)", (Fraction("0.49") * Fraction("0.7"), Fraction("0.49") * Fraction("0.9")),
         (Fraction("0.343"), Fraction("0.441")))

    low, width = steps[-1]
    tag = arith_codeword(low, width)
    same("tag of aabaaa", tag, "010110")
    same("tag length l", len(tag), 6)
    check("u times 2^6", float(low) * 64, "21.952")
    same("tag integer v", int(tag, 2), 22)
    v = int(tag, 2)
    same("dyadic interval of the tag", (Fraction(v, 64), Fraction(v + 1, 64)),
         (Fraction("0.34375"), Fraction("0.359375")))
    require("the dyadic interval lies inside [0.343, 0.376614)", low <= Fraction(v, 64) and Fraction(v + 1, 64) <= low + width)
    check("tag bits beyond the information", 6 - log2(1 / prob), "1.1052")
    decoded, values = arith_decode(tag, 6, PROBS, LETTERS)
    same("decoding the tag", decoded, msg)
    for k, stated in enumerate(["0.34375", "0.491071", "0.701531", "0.007653", "0.010933", "0.015618"]):
        check(f"decoder value z at step {k + 1}", values[k], stated)

    from itertools import product
    for n in range(1, 6):
        tags, total, inside, decodes = [], Fraction(0), True, True
        for letters in product(LETTERS, repeat=n):
            m = "".join(letters)
            lo, w = arith_intervals(m, PROBS, LETTERS)[-1]
            t = arith_codeword(lo, w)
            vv, ll = int(t, 2), len(t)
            inside &= lo <= Fraction(vv, 2 ** ll) and Fraction(vv + 1, 2 ** ll) <= lo + w
            decodes &= arith_decode(t, n, PROBS, LETTERS)[0] == m
            tags.append(t)
            total += w
        tags.sort()
        prefix_free = all(not tags[i + 1].startswith(tags[i]) for i in range(len(tags) - 1))
        require(f"n = {n}: intervals tile [0, 1), tags inside, prefix-free, decodable",
                total == 1 and inside and prefix_free and decodes)


def check_rates(block_rate) -> None:
    print("Why Arithmetic Coding Reaches the Entropy")
    hb = binary_entropy(0.1)
    check("h(0.1)", hb, "0.4690")
    rows = [(1, "2.3", "2.4690"), (2, "1.3", "1.4690"), (4, "0.8009", "0.9690"),
            (8, "0.6751", "0.7190"), (100, "0.4841", "0.4890"), (1000, "0.4705", "0.4710")]
    for n, rate, bound in rows:
        r = arith_expected_rate(0.1, n)
        check(f"arithmetic coding, n = {n}", r, rate)
        check(f"bound h(0.1) + 2/n, n = {n}", hb + 2 / n, bound)
        require(f"h(0.1) <= rate < bound, n = {n}", hb <= r < hb + 2 / n)
    for n, stated in [(1, "1"), (2, "0.6450"), (4, "0.4926"), (8, "0.4758")]:
        check(f"Huffman's code on blocks of {n}", block_rate(0.1, n), stated)
    same("tag length of the letter 0 alone", len(arith_codeword(Fraction(0), Fraction(9, 10))), 2)
    same("tag length of the letter 1 alone", len(arith_codeword(Fraction(9, 10), Fraction(1, 10))), 5)
    check("0.9 x 2 + 0.1 x 5", 0.9 * 2 + 0.1 * 5, "2.3")
    r6, r7 = arith_expected_rate(0.1, 6), arith_expected_rate(0.1, 7)
    check("arithmetic coding, n = 6", r6, "0.7114")
    check("arithmetic coding, n = 7", r7, "0.7286")
    require("the rate rises from n = 6 to n = 7", r7 > r6)
    check("distance from the entropy at n = 1000", arith_expected_rate(0.1, 1000) - hb, "0.0015")
    h01 = binary_entropy(0.01)
    check("h(0.01)", h01, "0.0808")
    check("1 bit per letter over h(0.01)", 1 / h01, "12.4")
    check("arithmetic coding at p = 0.01, n = 1000", arith_expected_rate(0.01, 1000), "0.0823")
    check("bound h(0.01) + 2/1000", h01 + 2 / 1000, "0.0828")


def check_range_coding() -> None:
    print("Range Coding and Adaptive Models")
    bits, lo, hi = settle_bits(Fraction(0), Fraction("0.49"))
    same("bits settled after aa", (bits, lo, hi), ("0", Fraction(0), Fraction("0.98")))
    w = hi - lo
    lo, hi = lo + w * PROBS["a"], lo + w * (PROBS["a"] + PROBS["b"])
    same("interval after the b", (lo, hi), (Fraction("0.686"), Fraction("0.882")))
    bits, lo, hi = settle_bits(lo, hi)
    same("bits settled after the b", (bits, lo, hi), ("1", Fraction("0.372"), Fraction("0.764")))
    for k in range(3):
        hi = lo + (hi - lo) * PROBS["a"]
        bits, lo, hi = settle_bits(lo, hi)
        same(f"bits settled by a number {k + 4}", bits, "")
    same("interval after aabaaa", (lo, hi), (Fraction("0.372"), Fraction("0.506456")))
    require("it straddles 1/2 inside [1/4, 3/4)", Fraction(1, 4) <= lo < Fraction(1, 2) < hi <= Fraction(3, 4))
    require("the settled bits 01 begin the tag 010110", "010110".startswith("01"))

    msg = "aabaaa"
    low, rng, ranges = 0, 1000, [1000]
    for s in msg:
        low, rng = range_step(low, rng, s)
        ranges.append(rng)
    same("ranges for aabaaa", ranges, [1000, 700, 490, 98, 68, 47, 32])
    same("final interval in thousandths", (low, low + rng), (343, 375))
    check("cost of truncation", log2(33.614 / 32), "0.0710")
    low, rng = 0, 1000
    for s in "aaba":
        low, rng = range_step(low, rng, s)
    same("interval after aaba", (low, low + rng), (343, 411))
    require("[343, 411) holds numbers beginning 3 and 4", str(low)[0] == "3" and str(low + rng - 1)[0] == "4")
    low, rng = range_step(low, rng, "a")
    same("interval after aabaa", (low, low + rng), (343, 390))
    require("every number in [343, 390) begins with 3", str(low)[0] == str(low + rng - 1)[0] == "3")
    no_overlap = all(r * CUM[s] // TOTAL + r * FREQ[s] // TOTAL <= r * (CUM[s] + FREQ[s]) // TOTAL
                     for r in range(1, 5000) for s in LETTERS)
    require("truncated slices never overlap, R = 1 to 4999", no_overlap)

    same("Laplace probability of aab", laplace_probability("aab", LETTERS), Fraction(1, 30))
    pl = laplace_probability(msg, LETTERS)
    same("Laplace probability of aabaaa", pl, Fraction(1, 168))
    same("the closed form agrees", laplace_closed_form(msg, LETTERS), pl)
    check("bits under Laplace's rule", log2(168), "7.3923")
    same("tag length under Laplace's rule", len(arith_codeword(Fraction(0), pl)), 9)
    from itertools import product
    holds = True
    for n in range(1, 7):
        for letters in product(LETTERS, repeat=n):
            m = "".join(letters)
            p_true = Fraction(1)
            for s in m:
                p_true *= PROBS[s]
            holds &= laplace_probability(m, LETTERS) * comb(n + 2, 2) >= p_true
            holds &= laplace_closed_form(m, LETTERS) == laplace_probability(m, LETTERS)
    require("price of learning holds for every message of 1 to 6 letters", holds)
    check("price bound at n = 100", log2(comb(102, 2)), "12.331")
    check("price bound at n = 1000", log2(comb(1002, 2)), "18.936")


def check_rans() -> None:
    print("Asymmetric Numeral Systems")
    msg = "aabaaa"
    check("log2 10", log2(10), "3.32")
    same("a-integer number 10", kth_state("a", 10), 13)
    same("pushing a onto 10", rans_push(10, "a"), 13)
    states = rans_encode_plain(msg, 10)
    same("states while encoding aabaaa backward", states, [10, 13, 16, 22, 117, 165, 234])
    qr = [divmod(states[i], FREQ[s]) for i, s in enumerate(reversed(msg))]
    same("(q, r) at each push", qr, [(1, 3), (1, 6), (2, 2), (11, 0), (16, 5), (23, 4)])
    for i, stated in enumerate(["1.3", "1.2308", "1.375", "5.3182", "1.4103", "1.4182"]):
        check(f"growth factor of push {i + 1}", states[i + 1] / states[i], stated)
    check("ideal factor for a", 10 / 7, "1.4286")
    check("ideal factor for b", 10 / 2, "5")
    same("b-integer number 22", kth_state("b", 22), 117)
    check("total growth log2(234/10)", log2(23.4), "4.5484")
    same("bits to write the state 234", (234).bit_length(), 8)
    same("slot of 234", 234 % TOTAL, 4)
    same("popping 234", rans_pop(234), ("a", 165))
    same("7 x 23 + 4", 7 * 23 + 4, 165)
    same("plain decoding from 234", rans_decode_plain(234, 6), (msg, [234, 165, 117, 22, 16, 13, 10]))
    inverse = all(rans_pop(rans_push(x, s)) == (s, x) for s in LETTERS for x in range(0, 3000))
    require("D(C(s, x)) = x for x = 0 to 2999", inverse)
    growth = all(abs(rans_push(x, s) * FREQ[s] / (TOTAL * x) - 1) < FREQ[s] / x
                 for s in LETTERS for x in range(1, 3000))
    require("|delta| < f_s / x for x = 1 to 2999", growth)

    x, digits, steps = rans_encode_renorm(msg)
    same("renormalized states after each push", [step[3] for step in steps], [13, 16, 22, 17, 23, 32])
    same("digits emitted before each push", [step[2] for step in steps], [[], [], [], [2], [], []])
    same("state sent and digits sent", (x, digits), (32, [2]))
    same("C(b, 2)", rans_push(2, "b"), 17)
    same("renormalized decoding", rans_decode_renorm(x, digits, 6), (msg, 10))
    in_range = True
    for s in LETTERS:
        for x0 in range(10, 100):
            y, emitted = x0, 0
            while y >= 10 * FREQ[s]:
                y //= 10
                emitted += 1
            in_range &= emitted <= 1 and 10 <= rans_push(y, s) < 100 and (rans_pop(rans_push(y, s))[1] < 10) == (emitted == 1)
    require("renormalization keeps every state in [10, 100)", in_range)
    rng = random.Random(6)
    for trial in range(20):
        m = "".join(rng.choices(LETTERS, weights=[7, 2, 1], k=300))
        x, digits, _ = rans_encode_renorm(m)
        ok = rans_decode_renorm(x, digits, len(m)) == (m, 10)
        ok &= rans_decode_plain(rans_encode_plain(m, 10)[-1], len(m))[0] == m
        if not ok:
            require(f"rANS round trip {trial}", False)
    require("rANS round trips, plain and renormalized, 20 messages", True)

    same("states while encoding ca", rans_encode_plain("ca", 10), [10, 13, 139])
    same("10 x 13 + 9 + 0", 10 * 13 + 9 + 0, 139)
    check("growth for ca", log2(13.9), "3.797")
    check("information in ca", log2(1 / 0.07), "3.8365")
    same("popping 139", rans_pop(139), ("c", 13))
    same("popping 13", rans_pop(13), ("a", 10))


def check_dictionary() -> None:
    print("Dictionary Coding with LZ77")
    text = "abcbcbcabcd"
    triples = lz77_encode(text)
    same("LZ77 parse of abcbcbcabcd", triples, [(0, 0, "a"), (0, 0, "b"), (0, 0, "c"), (2, 4, "a"), (3, 2, "d")])
    same("LZ77 decoding", lz77_decode(triples), text)
    same("offsets back to bc from position 9", [o for o in range(1, 9) if text[8 - o:10 - o] == "bc"], [3, 5, 7])
    require("the copy at position 4 is longer than its offset", triples[3][1] > triples[3][0])
    same("LZ77 parse of aaaaaaab", lz77_encode("aaaaaaab"), [(0, 0, "a"), (1, 6, "b")])
    same("its decoding", lz77_decode([(0, 0, "a"), (1, 6, "b")]), "aaaaaaab")

    print("LZ78 and LZW")
    pairs, _ = lz78_encode("aababcabcd")
    same("LZ78 pairs for aababcabcd", pairs, [(0, "a"), (1, "b"), (2, "c"), (3, "d")])
    same("LZ78 decoding", lz78_decode(pairs), "aababcabcd")
    same("B(4) on the phrases 1, 0, 11, 00", lz78_bits([1, 0, 1, 1, 0, 0]) * 6, 1 + 2 + 3 + 3)
    hb = binary_entropy(0.1)
    lengths = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]
    bits = binary_source(0.1, lengths[-1])
    rates = {n: lz78_bits(bits[:n]) for n in lengths}
    check("LZ78 after 100 letters", rates[100], "0.7600")
    check("LZ78 after 10,000 letters", rates[10000], "0.6237")
    check("LZ78 after 1,000,000 letters", rates[1000000], "0.5586")
    check("percent above the entropy at 1,000,000", 100 * (rates[1000000] / hb - 1), "19")
    scaled = [(rates[n] - hb) * log2(n) for n in lengths]
    check("smallest excess times log2 n", min(scaled), "1.68")
    check("largest excess times log2 n", max(scaled), "2.34")

    codes, steps = lzw_encode("abababa", "ab")
    same("LZW numbers for abababa", codes, [0, 1, 2, 4])
    same("LZW matches and new entries", [(st[0], st[3], st[4]) for st in steps],
         [("a", "ab", 2), ("b", "ba", 3), ("ab", "aba", 4), ("aba", "", None)])
    same("LZW decoding of 0 1 2 4", lzw_decode(codes, "ab")[:2], ("abababa", [4]))
    codes, _ = lzw_encode("aaaaaaa", "ab")
    same("LZW numbers for aaaaaaa", codes, [0, 2, 3, 0])
    same("its decoding and special cases", lzw_decode(codes, "ab")[:2], ("aaaaaaa", [2, 3]))
    codes, _ = lzw_encode("a" * 1000, "ab")
    same("numbers for a run of 1000 a's", len(codes), 45)
    same("44 x 45 / 2", 44 * 45 // 2, 990)
    same("the 45th phrase", codes[44], 10)
    decoded, special, table = lzw_decode(codes, "ab")
    same("decoding the run", decoded, "a" * 1000)
    same("special cases are the 2nd to the 44th numbers", special, codes[1:44])
    same("entries in the dictionary at the end", len(table), 46)
    require("every number fits in 6 bits", max(codes) < 2 ** 6 and len(table) <= 2 ** 6)
    same("bits for the run", 45 * 6, 270)
    same("bits to send the length 1000", (1000).bit_length(), 10)
    decoded, special, table = lzw_decode([1, 0, 2, 4, 1], "ab")
    same("decoding 1 0 2 4 1", (decoded, special), ("babababb", [4]))
    same("entries 2, 3 and 4", (table[2], table[3], table[4]), ("ba", "ab", "bab"))
    same("encoding babababb", lzw_encode("babababb", "ab")[0], [1, 0, 2, 4, 1])

    rng = random.Random(7)
    for trial in range(30):
        t = "".join(rng.choices("abc", weights=[6, 3, 1], k=400))
        ok = lz77_decode(lz77_encode(t)) == t and lz78_decode(lz78_encode(t)[0]) == t
        ok &= lzw_decode(lzw_encode(t, "abc")[0], "abc")[0] == t
        if not ok:
            require(f"dictionary round trip {trial}", False)
    require("LZ77, LZ78 and LZW round trips, 30 texts", True)


def check_runs_bwt_mtf() -> None:
    print("Run-Length Coding")
    row = "0" * 12 + "1" * 3 + "0" * 9 + "1" * 2
    same("pixels in the row", len(row), 26)
    same("runs of the row", run_lengths(row), [12, 3, 9, 2])
    require("each run fits a 4-bit field", all(r < 16 for r in run_lengths(row)))
    same("bits for four 4-bit runs", 4 * len(run_lengths(row)), 16)
    same("run decoding", run_decode([12, 3, 9, 2]), row)
    rng = random.Random(8)
    ok = all(run_decode(run_lengths(t)) == t
             for t in ("".join(rng.choice("01") for _ in range(200)) for _ in range(30)))
    require("run-length round trips, 30 rows", ok)
    p, hb = 0.1, binary_entropy(0.1)
    gap_entropy = -sum((1 - p) ** (r - 1) * p * ((r - 1) * log2(1 - p) + log2(p)) for r in range(1, 3000))
    check("entropy of a gap, summed term by term", gap_entropy, "4.6900")
    check("h(p) / p", hb / p, "4.6900")
    check("mean gap E[R]", sum(r * (1 - p) ** (r - 1) * p for r in range(1, 3000)), "10")
    check("h(p) + p", hb + p, "0.5690")

    print("The Burrows-Wheeler Transform")
    last, rows = bwt("couscous$")
    same("sorted rotations", rows, ["$couscous", "cous$cous", "couscous$", "ous$cousc", "ouscous$c",
                                    "s$couscou", "scous$cou", "us$cousco", "uscous$co"])
    same("transform of couscous", last, "ss$ccuuoo")
    same("first column", "".join(r[0] for r in rows), "$ccoossuu")
    same("last-to-first map", lf_map(last), [5, 6, 0, 1, 2, 7, 8, 3, 4])
    text, path = inverse_bwt(last)
    same("rows visited from row 0", path, [0, 5, 7, 3, 1, 6, 8, 4, 2])
    same("letters collected", "".join(last[r] for r in path[:-1]), "suocsuoc")
    same("inverse transform", text, "couscous$")
    first = sorted("ii$bbrree")
    same("first column of ii$bbrree", "".join(first), "$bbeeiirr")
    same("C for $, b, e, i, r", [first.index(ch) for ch in "$beir"], [0, 1, 3, 5, 7])
    same("its last-to-first map", lf_map("ii$bbrree"), [5, 6, 0, 1, 2, 7, 8, 3, 4])
    text, path = inverse_bwt("ii$bbrree")
    same("letters collected", "".join("ii$bbrree"[r] for r in path[:-1]), "irebireb")
    same("inverse of ii$bbrree", text, "beriberi$")
    same("forward transform of beriberi", bwt("beriberi$")[0], "ii$bbrree")
    rng = random.Random(9)
    ok = all(inverse_bwt(bwt(t + "$")[0])[0] == t + "$"
             for t in ("".join(rng.choices("abcd", k=60)) for _ in range(30)))
    require("Burrows-Wheeler round trips, 30 blocks", ok)

    print("Move-to-Front")
    require("the output has the letter counts of couscous", sorted(last) == sorted("couscous$"))
    idx, lists = mtf_encode(last, "$cosu")
    same("move-to-front output", idx, [3, 0, 1, 2, 0, 4, 0, 4, 0])
    same("lists before each step", lists, ["$cosu", "s$cou", "s$cou", "$scou", "c$sou", "c$sou",
                                           "uc$so", "uc$so", "ouc$s"])
    same("move-to-front decoding", mtf_decode(idx, "$cosu"), last)
    check("empirical entropy of the transform output", empirical_entropy(last), "2.2810")
    check("empirical entropy of the move-to-front output", empirical_entropy(idx), "2.0588")
    same("decoding 1 0 2 0 0 2 from abc", mtf_decode([1, 0, 2, 0, 0, 2], "abc"), "bbccca")
    same("encoding bbccca", mtf_encode("bbccca", "abc")[0], [1, 0, 2, 0, 0, 2])
    ok = all(mtf_decode(mtf_encode(t, "abcd")[0], "abcd") == t
             for t in ("".join(rng.choices("abcd", k=60)) for _ in range(30)))
    require("move-to-front round trips, 30 texts", ok)


def check_lossy() -> None:
    print("Quantization")
    xs = [0.37, -1.12, 2.05, 0.61, -0.28]
    q = [quantize(x, 0.5) for x in xs]
    same("indices", [k for k, _ in q], [1, -2, 4, 1, -1])
    same("reconstructions", [r for _, r in q], [0.5, -1.0, 2.0, 0.5, -0.5])
    check("mean squared error", sum((x - r) ** 2 for x, (_, r) in zip(xs, q)) / len(xs), "0.01886")
    check("Delta^2 / 12 at Delta = 0.5", 0.5 ** 2 / 12, "0.020833")
    check("decibels per bit", db(4), "6.0206")
    check("signal-to-noise ratio at 8 bits", db(2 ** 16), "48.16")
    check("signal-to-noise ratio at 16 bits", db(2 ** 32), "96.33")

    print("Transform Coding")
    x1, x2 = 100, 104
    s, d = (x1 + x2) / sqrt(2), (x1 - x2) / sqrt(2)
    check("s", s, "144.2498")
    check("d", d, "-2.8284")
    same("energy x1^2 + x2^2", x1 * x1 + x2 * x2, 20816)
    check("energy s^2 + d^2", s * s + d * d, "20816")
    check("fraction of the energy in s", s * s / (x1 * x1 + x2 * x2), "0.99962")
    sq, dq = quantize(s, 4)[1], quantize(d, 8)[1]
    same("s and d quantized", (sq, dq), (144, 0))
    y1, y2 = (sq + dq) / sqrt(2), (sq - dq) / sqrt(2)
    check("reconstructed pixels", y1, "101.8234")
    require("both pixels come back equal", y1 == y2)
    check("squared error over the pixels", (x1 - y1) ** 2 + (x2 - y2) ** 2, "8.0624")
    check("squared error over the coefficients", (s - sq) ** 2 + (d - dq) ** 2, "8.0624")

    print("The Rate-Distortion Trade")
    y = sqrt(2 / pi)
    check("best one-bit level E|X|", y, "0.7979")
    check("its distortion 1 - 2/pi", 1 - 2 / pi, "0.36338")
    check("the same, integrated over the density", fixed_rate_mse(2, 2 * y), "0.36338")
    check("in decibels", db(1 - 2 / pi), "-4.40")
    check("the bound at one bit", db(2 ** -2), "-6.02")
    table = [(1, "1.5958", "0.36338", "-4.40", "-6.02", "1.62"), (2, "0.9957", "0.11885", "-9.25", "-12.04", "2.79"),
             (3, "0.5860", "0.03744", "-14.27", "-18.06", "3.80"), (4, "0.3352", "0.01154", "-19.38", "-24.08", "4.71"),
             (5, "0.1881", "0.00350", "-24.57", "-30.10", "5.54")]
    for rate, step_s, d_s, db_s, bound_s, gap_s in table:
        step, dist = best_uniform(2 ** rate)
        check(f"best step, {2 ** rate} levels", step, step_s)
        check(f"distortion, {2 ** rate} levels", dist, d_s)
        check(f"distortion in dB, {2 ** rate} levels", db(dist), db_s)
        check(f"bound in dB at {rate} bits", db(2 ** (-2 * rate)), bound_s)
        check(f"gap in dB at {rate} bits", db(dist) - db(2 ** (-2 * rate)), gap_s)
    check("differential entropy of a unit Gaussian", 0.5 * log2(2 * pi * exp(1)), "2.0471")
    check("limiting gap in bits", 0.5 * log2(pi * exp(1) / 6), "0.2546")
    check("limiting distortion ratio pi e / 6", pi * exp(1) / 6, "1.4233")
    check("limiting gap in dB", db(pi * exp(1) / 6), "1.5329")
    check("limiting gap in dB, as the caption rounds it", db(pi * exp(1) / 6), "1.53")
    for step, stated in [(1, "1.88"), (0.5, "1.62"), (0.25, "1.56")]:
        rate, mse = entropy_coded_uniform(step)
        check(f"entropy-coded gap in dB at step {step}", db(mse) - db(2 ** (-2 * rate)), stated)
    check("index entropy at step 0.25", entropy_coded_uniform(0.25)[0], "4.0508")
    check("R(0.25) for variance 4", 0.5 * log2(4 / 0.25), "2")
    check("D(3) for variance 4", 4 * 2 ** -6, "0.0625")
    check("a quarter of D(2)", (4 * 2 ** -6) / (4 * 2 ** -4), "0.25")


def check_memory() -> None:
    print("Compression Inside Memory Systems")
    ptrs = [0x00403A10, 0x00403A18, 0x00403A20, 0x00403A28, 0x00403A40, 0x00403A48, 0x00403A60, 0x00403A70]
    require("only the last byte varies", all(w >> 8 == ptrs[0] >> 8 for w in ptrs))
    base, deltas, mask, size = bdi_compress(ptrs, 4, 1, False)
    same("base", format(base, "08X"), "00403A10")
    same("deltas", deltas, [0, 8, 16, 24, 48, 56, 80, 96])
    same("deltas in hexadecimal", [format(x, "02X") for x in deltas], ["00", "08", "10", "18", "30", "38", "50", "60"])
    same("compressed bytes", size, 12)
    check("compression ratio", 32 / size, "2.67")
    same("decompression", bdi_decompress(base, deltas, mask), ptrs)
    flipped = base ^ (1 << 4)
    same("base after bit 4 flips", format(flipped, "08X"), "00403A00")
    same("error in each word", [a - b for a, b in zip(bdi_decompress(flipped, deltas, mask), ptrs)], [-16] * 8)
    same("stored bits", 8 * size, 96)
    same("bits in the line", 8 * 32, 256)

    mixed = [0, 0x00403A10, 3, 0x00403A28, 1, 0x00403A18, 0xC, 0x00403A40]
    require("no 1-byte or 2-byte delta reaches every word from one base",
            bdi_compress(mixed, 4, 1, False) is None and bdi_compress(mixed, 4, 2, False) is None)
    base, deltas, mask, size = bdi_compress(mixed, 4, 1, True)
    same("deltas with zero as a second base", deltas, [0, 0, 3, 24, 1, 8, 12, 48])
    same("mask", "".join(map(str, mask)), "01010101")
    same("compressed bytes", size, 4 + 8 + 1)
    check("compression ratio", 32 / size, "2.46")
    same("decompression", bdi_decompress(base, deltas, mask), mixed)

    print("Worked Example: Base-Delta Compression of Timestamps")
    start = 0x17F3A2C1000
    stamps = [start + x for x in (0, 0x120, 0x250, 0x3A8, 0x400, 0x5F0, 0x710, 0x7F8)]
    require("1-byte deltas fail", bdi_compress(stamps, 8, 1, False) is None)
    base, deltas, mask, size = bdi_compress(stamps, 8, 2, False)
    same("deltas", deltas, [0, 288, 592, 936, 1024, 1520, 1808, 2040])
    require("the second delta exceeds 127", deltas[1] > 127)
    same("compressed bytes", size, 24)
    check("compression ratio", 64 / size, "2.67")
    same("bytes with 4-byte deltas", bdi_compress(stamps, 8, 4, False)[3], 40)
    same("decompression", bdi_decompress(base, deltas, mask), stamps)
    require("bit 12 is set in the base", (base >> 12) & 1 == 1)
    flipped = base ^ (1 << 12)
    same("base after bit 12 flips", format(flipped, "X"), "17F3A2C0000")
    same("error in each timestamp", [a - b for a, b in zip(bdi_decompress(flipped, deltas, mask), stamps)], [-4096] * 8)
    same("stored bits", 8 * size, 192)
    same("bits in the line", 8 * 64, 512)


def check_examples() -> None:
    print("Worked Example: Arithmetic Coding of a Sparse Binary Message")
    bp = {"0": Fraction(9, 10), "1": Fraction(1, 10)}
    msg = "0000100000"
    steps = arith_intervals(msg, bp, "01")
    same("interval after four 0s", (steps[4][0], steps[4][0] + steps[4][1]), (Fraction(0), Fraction("0.6561")))
    same("interval after the 1", (steps[5][0], steps[5][0] + steps[5][1]), (Fraction("0.59049"), Fraction("0.6561")))
    low, width = steps[10]
    check("final upper end", low + width, "0.629232")
    check("final width", width, "0.038742")
    check("information in the message", log2(1 / width), "4.690")
    check("10 h(0.1)", 10 * binary_entropy(0.1), "4.690")
    tag = arith_codeword(low, width)
    same("tag", tag, "100110")
    check("u times 2^6", float(low) * 64, "37.79")
    v = int(tag, 2)
    same("tag integer v", v, 38)
    same("dyadic interval", (Fraction(v, 64), Fraction(v + 1, 64)), (Fraction("0.59375"), Fraction("0.609375")))
    require("inside the message interval", low <= Fraction(v, 64) and Fraction(v + 1, 64) <= low + width)
    same("decoding", arith_decode(tag, 10, bp, "01")[0], msg)

    print("Worked Example: Learning the Model as It Goes")
    for n, true_s, price_s, per_s, bound_s in [(100, "115.68", "6.096", "0.061", "12.331"),
                                               (1000, "1156.78", "9.396", "0.0094", "18.936")]:
        ka, kb, kc = 7 * n // 10, 2 * n // 10, n // 10
        true_bits = ka * log2(10 / 7) + kb * log2(10 / 2) + kc * log2(10)
        check(f"bits from the true model, n = {n}", true_bits, true_s)
        check(f"n H(X), n = {n}", n * entropy(PROBS.values()), true_s)
        message = "a" * ka + "b" * kb + "c" * kc
        laplace_bits = bits_of(laplace_closed_form(message, LETTERS))
        price = laplace_bits - true_bits
        check(f"price of learning, n = {n}", price, price_s)
        check(f"price per letter, n = {n}", price / n, per_s)
        bound = log2(comb(n + 2, 2))
        check(f"guarantee, n = {n}", bound, bound_s)
        require(f"price is about half the guarantee, n = {n}", 0.45 < price / bound < 0.55)
    p100 = bits_of(laplace_closed_form("a" * 70 + "b" * 20 + "c" * 10, LETTERS)) - 100 * entropy(PROBS.values())
    p1000 = bits_of(laplace_closed_form("a" * 700 + "b" * 200 + "c" * 100, LETTERS)) - 1000 * entropy(PROBS.values())
    check("growth of the price from 100 to 1000 letters", p1000 - p100, "3.30")

    print("Exercise: Arithmetic coding of ba")
    low, width = arith_intervals("ba", PROBS, LETTERS)[-1]
    same("interval of ba", (low, low + width, width), (Fraction("0.7"), Fraction("0.84"), Fraction("0.14")))
    check("information in ba", log2(1 / 0.14), "2.8365")
    tag = arith_codeword(low, width)
    same("tag", tag, "1100")
    check("u times 2^4", float(low) * 16, "11.2")
    v = int(tag, 2)
    same("dyadic interval", (Fraction(v, 16), Fraction(v + 1, 16)), (Fraction("0.75"), Fraction("0.8125")))
    require("inside the interval of ba", low <= Fraction(v, 16) and Fraction(v + 1, 16) <= low + width)
    same("Huffman bits on ba", len("10" + "0"), 3)

    print("Looking Ahead")
    same("110110 decodes to", arith_decode("110110", 6, PROBS, LETTERS)[0], "bbaaaa")
    same("011110 decodes to", arith_decode("011110", 6, PROBS, LETTERS)[0], "aacaba")
    every = all(len(arith_decode(format(t, "06b"), 6, PROBS, LETTERS)[0]) == 6 for t in range(64))
    require("every 6-bit string decodes to some message", every)
    hb = binary_entropy(0.1)
    check("capacity 1 - h(0.1)", 1 - hb, "0.5310")
    require("0.4690 < 0.5310", hb < 1 - hb)


def main() -> None:
    from ch05_source_coding_and_prefix_codes import block_rate
    check_intervals()
    check_rates(block_rate)
    check_range_coding()
    check_rans()
    check_dictionary()
    check_runs_bwt_mtf()
    check_lossy()
    check_memory()
    check_examples()
    print("All chapter 6 numbers reproduced.")


if __name__ == "__main__":
    main()
