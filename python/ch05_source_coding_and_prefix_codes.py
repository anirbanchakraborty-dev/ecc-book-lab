"""Reproduce every number in chapter 5 of Error Correction Coding, "Source Coding and Prefix Codes".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch05_source_coding_and_prefix_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. Codewords, parses
and decoded strings must match the chapter exactly. A run that reaches the
final line has reproduced every number in the chapter.

The module also holds the Huffman, Shannon-code and block-coding functions that
ch05_charts.py plots, so the charts and this check share one source.
"""

import sys
from fractions import Fraction
from itertools import product
from math import log2


# Helpers ------------------------------------------------------------------------

def half_unit(stated: str) -> float:
    text = stated.replace(",", "")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", ""))
    print(f"  {label:<60} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<60} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<60} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH in '{label}'")
        sys.exit(1)


# Sources, entropy and average length --------------------------------------------

def entropy(probs) -> float:
    """Entropy in bits of a distribution given as a list of probabilities."""
    return sum(float(p) * log2(1 / float(p)) for p in probs if p > 0)


def average_length(probs, lengths):
    return sum(p * l for p, l in zip(probs, lengths))


def length_variance(probs, lengths):
    mean = average_length(probs, lengths)
    return sum(p * l * l for p, l in zip(probs, lengths)) - mean * mean


def relative_entropy(p_dist, q_dist) -> float:
    return sum(float(p) * log2(float(p) / float(q)) for p, q in zip(p_dist, q_dist) if p > 0)


def kraft_sum(lengths) -> Fraction:
    return sum(Fraction(1, 2 ** l) for l in lengths)


# Codes as lists of codewords ----------------------------------------------------

def is_prefix_free(words) -> bool:
    return not any(a != b and b.startswith(a) for a in words for b in words)


def parses(bits: str, words: dict) -> list:
    """Every way to split a bit string into codewords, as strings of letters."""
    if bits == "":
        return [""]
    found = []
    for letter, word in words.items():
        if bits.startswith(word):
            found += [letter + rest for rest in parses(bits[len(word):], words)]
    return sorted(found)


def decode_prefix(bits: str, words: dict):
    """Decode left to right with a prefix code, returning the letters and any leftover bits."""
    inverse = {w: c for c, w in words.items()}
    letters, current = [], ""
    for b in bits:
        current += b
        if current in inverse:
            letters.append(inverse[current])
            current = ""
    return "".join(letters), current


def encode(message: str, words: dict) -> str:
    return "".join(words[c] for c in message)


def sardinas_patterson(words):
    """Return (uniquely decodable, list of dangling-suffix sets in order)."""
    code = set(words)

    def dangling(a_set, b_set):
        out = set()
        for a in a_set:
            for b in b_set:
                if a != b and b.startswith(a):
                    out.add(b[len(a):])
        return out

    current = dangling(code, code)
    history, seen = [], []
    while current:
        history.append(sorted(current))
        if current & code:
            return False, history
        if current in seen:
            return True, history
        seen.append(current)
        current = dangling(code, current) | dangling(current, code)
    return True, history


def binary_digits(x: Fraction, n: int) -> str:
    """The first n bits after the binary point of a number in [0, 1)."""
    out = []
    for _ in range(n):
        x *= 2
        bit = int(x >= 1)
        out.append(str(bit))
        x -= bit
    return "".join(out)


def kraft_code(lengths) -> list:
    """The prefix code built in the proof of Kraft's inequality, in the order the lengths are given."""
    order = sorted(range(len(lengths)), key=lambda i: lengths[i])
    words = [None] * len(lengths)
    point = Fraction(0)
    for i in order:
        words[i] = binary_digits(point, lengths[i])
        point += Fraction(1, 2 ** lengths[i])
    return words


# Shannon's code, Fano's splitting and Huffman's algorithm -----------------------

def shannon_lengths(probs) -> list:
    """The smallest l with 2^-l at most p, which is the ceiling of log2(1/p), computed exactly."""
    lengths = []
    for p in probs:
        l = 0
        while Fraction(1, 2 ** l) > p:
            l += 1
        lengths.append(l)
    return lengths


def shannon_code(probs) -> list:
    """Shannon's construction: probabilities in decreasing order, codeword i the first l_i bits of F_i."""
    lengths = shannon_lengths(probs)
    words, cumulative = [], Fraction(0)
    for p, l in zip(probs, lengths):
        words.append(binary_digits(cumulative, l))
        cumulative += p
    return words


def fano_code(probs):
    """Fano's splitting of a list in decreasing order. Returns the codewords and whether a split was tied."""
    words = [""] * len(probs)
    tied = [False]

    def split(indices):
        if len(indices) == 1:
            return
        total = sum(probs[i] for i in indices)
        best, best_gap, upper = None, None, Fraction(0)
        gaps = []
        for k in range(1, len(indices)):
            upper += probs[indices[k - 1]]
            gap = abs(2 * upper - total)
            gaps.append(gap)
            if best_gap is None or gap < best_gap:
                best, best_gap = k, gap
        if gaps.count(best_gap) > 1:
            tied[0] = True
        for i in indices[:best]:
            words[i] += "0"
        for i in indices[best:]:
            words[i] += "1"
        split(indices[:best])
        split(indices[best:])

    split(list(range(len(probs))))
    return words, tied[0]


def huffman(probs, merged_high: bool = True):
    """Huffman's algorithm by the sorted-list method.

    The list is kept in decreasing order of probability and the last two entries are
    merged at each step. A merged entry is placed above every entry of equal
    probability when merged_high is true and below them otherwise, which is the only
    choice the algorithm leaves open. Returns the codeword lengths in the order of
    probs and the list of merged probabilities in the order they were formed.
    """
    nodes = [(p, (i,)) for i, p in sorted(enumerate(probs), key=lambda t: -t[1])]
    lengths = [0] * len(probs)
    merges = []
    while len(nodes) > 1:
        (p1, m1), (p2, m2) = nodes[-2], nodes[-1]
        nodes = nodes[:-2]
        weight, members = p1 + p2, m1 + m2
        for i in members:
            lengths[i] += 1
        merges.append(weight)
        if merged_high:
            spot = next((k for k, (p, _) in enumerate(nodes) if p <= weight), len(nodes))
        else:
            spot = next((k for k, (p, _) in enumerate(nodes) if p < weight), len(nodes))
        nodes.insert(spot, (weight, members))
    return lengths, merges


def huffman_length(probs) -> float:
    """Average Huffman codeword length, as the sum of the merged probabilities."""
    return float(sum(huffman([float(p) for p in probs])[1]))


# Blocks of a binary source ------------------------------------------------------

def block_probs(p_one, n: int) -> list:
    """Probabilities of the 2^n blocks of n independent letters, each 1 with probability p_one."""
    return [p_one ** sum(b) * (1 - p_one) ** (n - sum(b)) for b in product((0, 1), repeat=n)]


def block_rate(p_one, n: int) -> float:
    """Huffman bits per letter when blocks of n letters are coded."""
    return huffman_length(block_probs(p_one, n)) / n


def binary_entropy(p: float) -> float:
    return entropy([p, 1 - p])


# McMillan -----------------------------------------------------------------------

def mcmillan_first_failure(lengths) -> int:
    """Smallest n with K^n above n(l_max - 1) + 1, the bound every uniquely decodable code obeys."""
    k, top = kraft_sum(lengths), max(lengths)
    n = 1
    while k ** n <= n * (top - 1) + 1:
        n += 1
    return n


# The chapter --------------------------------------------------------------------

FOUR = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
CODE_FIXED = {"a": "00", "b": "01", "c": "10", "d": "11"}
CODE_PREFIX = {"a": "0", "b": "10", "c": "110", "d": "111"}
CODE_SUFFIX = {"a": "0", "b": "01", "c": "011", "d": "111"}
CODE_AMBIG = {"a": "0", "b": "01", "c": "100", "d": "101"}
CODE_SHORT = {"a": "0", "b": "1", "c": "01", "d": "10"}

SIX = [Fraction(n, 100) for n in (45, 20, 10, 10, 10, 5)]
LETTERS6 = "abcdef"
TREE_A = {"a": "0", "b": "10", "c": "1100", "d": "1101", "e": "1110", "f": "1111"}
TREE_B = {"a": "0", "b": "110", "c": "100", "d": "101", "e": "1110", "f": "1111"}


def lengths_of(code: dict) -> list:
    return [len(w) for w in code.values()]


def siblings(u: str, v: str) -> bool:
    return len(u) == len(v) and u[:-1] == v[:-1] and u != v


def main() -> None:
    print("Codes for a Source")
    check("four-letter source entropy H", entropy(FOUR), "1.75")
    check("fixed-length code 00 01 10 11, average length", average_length(FOUR, lengths_of(CODE_FIXED)), "2")
    check("code 0 10 110 111, average length", average_length(FOUR, lengths_of(CODE_PREFIX)), "1.75")
    check("fixed-length code, bits above entropy", average_length(FOUR, [2] * 4) - entropy(FOUR), "0.25")
    check("code 0 1 01 10, average length", average_length(FOUR, lengths_of(CODE_SHORT)), "1.25")
    check("code 0 1 01 10, average below entropy by", entropy(FOUR) - float(average_length(FOUR, lengths_of(CODE_SHORT))), "0.5")
    same("code 0 1 01 10, readings of 01", parses("01", CODE_SHORT), ["ab", "c"])
    same("code 0 1 01 10, readings of 0110", parses("0110", CODE_SHORT), ["abba", "abd", "cba", "cd"])

    print("\nWhen a Code Can Be Decoded")
    for name, code in (("0 10 110 111", CODE_PREFIX), ("0 01 011 111", CODE_SUFFIX), ("0 01 100 101", CODE_AMBIG)):
        check(f"code {name}, Kraft sum", kraft_sum(lengths_of(code)), "1")
        check(f"code {name}, average length", average_length(FOUR, lengths_of(code)), "1.75")
    same("code 0 01 100 101, readings of 0101", parses("0101", CODE_AMBIG), ["ad", "bb"])
    ud, hist = sardinas_patterson(CODE_AMBIG.values())
    same("code 0 01 100 101, dangling suffixes", hist, [["1"], ["00", "01"]])
    require("code 0 01 100 101 is not uniquely decodable", not ud)
    ud, hist = sardinas_patterson(CODE_SUFFIX.values())
    require("code 0 01 011 111 is uniquely decodable", ud)
    require("code 0 01 011 111 is not a prefix code", not is_prefix_free(CODE_SUFFIX.values()))
    require("code 0 10 110 111 is a prefix code", is_prefix_free(CODE_PREFIX.values()))
    same("code 0 01 011 111, readings of 0 then six 1s", parses("0" + "1" * 6, CODE_SUFFIX), ["add"])
    same("code 0 01 011 111, readings of 0 then seven 1s", parses("0" + "1" * 7, CODE_SUFFIX), ["bdd"])
    same("code 0 01 011 111, readings of 0 then eight 1s", parses("0" + "1" * 8, CODE_SUFFIX), ["cdd"])

    print("\nPrefix Codes and Code Trees")
    same("decode 1100100111 with 0 10 110 111", decode_prefix("1100100111", CODE_PREFIX), ("cabad", ""))
    same("encode cabad", encode("cabad", CODE_PREFIX), "1100100111")

    print("\nThe Kraft Inequality")
    seven = [2, 2, 3, 3, 3, 4, 4]
    check("lengths 2 2 3 3 3 4 4, Kraft sum", kraft_sum(seven), "1")
    same("depth-4 leaves below each codeword", [2 ** (4 - l) for l in seven], [4, 4, 2, 2, 2, 1, 1])
    same("total depth-4 leaves", sum(2 ** (4 - l) for l in seven), 16)
    same("Kraft construction for 2 2 3 3 3 4 4", kraft_code(seven), ["00", "01", "100", "101", "110", "1110", "1111"])
    same("starting points in sixteenths", [int(sum(Fraction(1, 2 ** l) for l in seven[:i]) * 16) for i in range(7)], [0, 4, 8, 10, 12, 14, 15])
    four_lengths = [1, 3, 3, 4]
    check("lengths 1 3 3 4, Kraft sum", kraft_sum(four_lengths), "0.8125")
    same("Kraft construction for 1 3 3 4", kraft_code(four_lengths), ["0", "100", "101", "1100"])
    same("unused depth-4 leaves for 1 3 3 4", 16 - sum(2 ** (4 - l) for l in four_lengths), 3)
    check("lengths 1 2 2 3, Kraft sum", kraft_sum([1, 2, 2, 3]), "1.125")

    print("\nMcMillan's Inequality")
    counts = {}
    for x, y in product(CODE_PREFIX.values(), repeat=2):
        counts[len(x) + len(y)] = counts.get(len(x) + len(y), 0) + 1
    same("code 0 10 110 111, pairs by encoded length 2..6", [counts[m] for m in range(2, 7)], [1, 2, 5, 4, 4])
    check("sum of A_m 2^-m for pairs", sum(Fraction(a, 2 ** m) for m, a in counts.items()), "1")
    k = kraft_sum([1, 1, 2])
    check("lengths 1 1 2, Kraft sum", k, "1.25")
    check("1.25^11", k ** 11, "11.64")
    check("1.25^12", k ** 12, "14.55")
    same("first n with 1.25^n above n + 1", mcmillan_first_failure([1, 1, 2]), 12)

    print("\nThe Entropy Lower Bound")
    check("four-letter source, D(P||U)", relative_entropy(FOUR, [Fraction(1, 4)] * 4), "0.25")
    h6 = entropy(SIX)
    check("six-letter source entropy", h6, "2.1955")
    sh_len = shannon_lengths(SIX)
    k_sh = kraft_sum(sh_len)
    check("Shannon code, Kraft sum", k_sh, "0.59375")
    same("Shannon code, Kraft sum as a fraction", k_sh, Fraction(19, 32))
    check("Shannon code, log2(1/K)", log2(1 / float(k_sh)), "0.7521")
    q_sh = [Fraction(1, 2 ** l) / k_sh for l in sh_len]
    check("Shannon code, D(P||Q)", relative_entropy(SIX, q_sh), "0.0025")
    check("Shannon code, L - H", float(average_length(SIX, sh_len)) - h6, "0.7545")
    check("Shannon code, D + log2(1/K)", relative_entropy(SIX, q_sh) + log2(1 / float(k_sh)), "0.7545")

    print("\nShannon's Code and Fano's Splitting")
    for p, stated in zip((0.45, 0.20, 0.10, 0.05), ("1.152", "2.322", "3.322", "4.322")):
        check(f"surprise of a letter with p = {p}", log2(1 / p), stated)
    same("Shannon lengths", sh_len, [2, 3, 4, 4, 4, 5])
    check("Shannon code, average length", average_length(SIX, sh_len), "2.95")
    check("H + 1", h6 + 1, "3.1955")
    cumulative = [sum(SIX[:i]) for i in range(6)]
    same("cumulative probabilities F in hundredths", [int(c * 100) for c in cumulative], [0, 45, 65, 75, 85, 95])
    same("binary expansion of 0.45, eight bits", binary_digits(Fraction(45, 100), 8), "01110011")
    same("binary expansion of 0.65, eight bits", binary_digits(Fraction(65, 100), 8), "10100110")
    same("binary expansion of 0.85, eight bits", binary_digits(Fraction(85, 100), 8), "11011001")
    same("binary expansion of 0.95, eight bits", binary_digits(Fraction(95, 100), 8), "11110011")
    sh_words = shannon_code(SIX)
    same("Shannon codewords", sh_words, ["00", "011", "1010", "1100", "1101", "11110"])
    require("Shannon codewords form a prefix code", is_prefix_free(sh_words))
    fano_words, fano_tied = fano_code(SIX)
    same("Fano codewords", fano_words, ["0", "100", "101", "110", "1110", "1111"])
    require("no Fano split of the six-letter source is tied", not fano_tied)
    check("Fano code, average length", average_length(SIX, [len(w) for w in fano_words]), "2.25")
    check("Fano first split, upper group a", SIX[0], "0.45")
    check("Fano second split, b and c", SIX[1] + SIX[2], "0.30")
    check("Fano second split, d e f", SIX[3] + SIX[4] + SIX[5], "0.25")
    check("Fano third split, e and f", SIX[4] + SIX[5], "0.15")

    print("\nHuffman's Algorithm")
    len_b, merges_b = huffman(SIX, merged_high=True)
    len_a, merges_a = huffman(SIX, merged_high=False)
    same("merges, merged entry placed high", [float(m) for m in merges_b], [0.15, 0.2, 0.35, 0.55, 1.0])
    same("merges, merged entry placed low", [float(m) for m in merges_a], [0.15, 0.2, 0.35, 0.55, 1.0])
    same("lengths, merged entry placed low (tree A)", len_a, [1, 2, 4, 4, 4, 4])
    same("lengths, merged entry placed high (tree B)", len_b, [1, 3, 3, 3, 4, 4])
    check("tree A average length", average_length(SIX, len_a), "2.25")
    check("tree B average length", average_length(SIX, len_b), "2.25")
    check("sum of merged probabilities", sum(merges_b), "2.25")
    check("tree A, 0.45 + 2(0.20) + 4(0.35)", Fraction(45, 100) + 2 * Fraction(20, 100) + 4 * Fraction(35, 100), "2.25")
    check("tree B, 0.45 + 3(0.40) + 4(0.15)", Fraction(45, 100) + 3 * Fraction(40, 100) + 4 * Fraction(15, 100), "2.25")
    check("tree A variance of lengths", length_variance(SIX, len_a), "1.7875")
    check("tree B variance of lengths", length_variance(SIX, len_b), "1.3875")
    same("tree A codewords have tree A lengths", [len(TREE_A[c]) for c in LETTERS6], len_a)
    same("tree B codewords have tree B lengths", [len(TREE_B[c]) for c in LETTERS6], len_b)
    require("tree A and tree B codewords are prefix codes", is_prefix_free(TREE_A.values()) and is_prefix_free(TREE_B.values()))
    require("tree A: c,d siblings, e,f siblings, cd and ef siblings, b and cdef siblings",
            siblings(TREE_A["c"], TREE_A["d"]) and siblings(TREE_A["e"], TREE_A["f"])
            and siblings(TREE_A["c"][:-1], TREE_A["e"][:-1]) and siblings(TREE_A["b"], TREE_A["c"][:-2]))
    require("tree B: c,d siblings, e,f siblings, b and ef siblings, cd and bef siblings",
            siblings(TREE_B["c"], TREE_B["d"]) and siblings(TREE_B["e"], TREE_B["f"])
            and siblings(TREE_B["b"], TREE_B["e"][:-1]) and siblings(TREE_B["c"][:-1], TREE_B["b"][:-1]))
    check("Huffman redundancy L - H", 2.25 - h6, "0.0545")
    check("Huffman efficiency H / L", h6 / 2.25, "0.976")
    same("Huffman on the four-letter source", huffman(FOUR)[0], [1, 2, 3, 3])

    print("\nShannon's Source Coding Theorem")
    hb = binary_entropy(0.1)
    check("h(0.1)", hb, "0.4690")
    pair = block_probs(Fraction(1, 10), 2)
    same("pair probabilities", sorted([float(p) for p in pair], reverse=True), [0.81, 0.09, 0.09, 0.01])
    pair_sorted = sorted(pair, reverse=True)
    same("pair Huffman lengths", huffman(pair_sorted)[0], [1, 2, 3, 3])
    same("pair merges", [float(m) for m in huffman(pair_sorted)[1]], [0.1, 0.19, 1.0])
    l2 = sum(huffman(pair_sorted)[1])
    check("L_2", l2, "1.29")
    check("L_2 / 2", l2 / 2, "0.645")
    triple = sorted(block_probs(Fraction(1, 10), 3), reverse=True)
    l3 = sum(huffman(triple)[1])
    same("triple merges", [float(m) for m in huffman(triple)[1]], [0.01, 0.018, 0.028, 0.109, 0.162, 0.271, 1.0])
    check("L_3", l3, "1.598")
    check("L_3 / 3", l3 / 3, "0.5327")
    quad = sorted(block_probs(Fraction(1, 10), 4), reverse=True)
    l4 = sum(huffman(quad)[1])
    check("L_4", l4, "1.9702")
    check("L_4 / 4", l4 / 4, "0.4926")
    for n, stated in ((1, "1.4690"), (2, "0.9690"), (3, "0.8023"), (4, "0.7190")):
        check(f"H + 1/{n}", hb + 1 / n, stated)
    for n, stated in ((5, "0.4802"), (6, "0.4702"), (7, "0.4743"), (8, "0.4758")):
        check(f"L_{n} / {n}", block_rate(0.1, n), stated)
    require("L_7 / 7 exceeds L_6 / 6, so the rate is not monotone in n", block_rate(0.1, 7) > block_rate(0.1, 6))
    for n, block_total, per_roll in ((1, "2.6667", "2.6667"), (2, "5.2222", "2.6111"), (3, "7.8148", "2.6049")):
        die_block = [Fraction(1, 6 ** n)] * (6 ** n)
        lengths_die, merges_die = huffman(die_block)
        check(f"fair die, blocks of {n}, L_{n}", sum(merges_die), block_total)
        check(f"fair die, blocks of {n}, per roll", sum(merges_die) / n, per_roll)
        same(f"fair die, blocks of {n}, lengths are the two nearest whole numbers to log2 6^{n}",
             sorted(set(lengths_die)), [int(n * log2(6)), int(n * log2(6)) + 1])
    same("fair die, pairs: codewords of length 6", sum(1 for l in huffman([Fraction(1, 36)] * 36)[0] if l == 6), 8)
    same("fair die, triples: codewords of length 8", sum(1 for l in huffman([Fraction(1, 216)] * 216)[0] if l == 8), 176)
    check("fair die entropy log2 6", log2(6), "2.5850")

    print("\nCompression and Its Fragility")
    message = "fadedcab"
    bits = encode(message, TREE_B)
    same("fadedcab encoded with tree B", bits, "1111010111101011000110")
    same("number of bits", len(bits), 22)
    flipped = bits[:2] + ("1" if bits[2] == "0" else "0") + bits[3:]
    same("third bit flipped", flipped, "1101010111101011000110")
    same("decoded after the flip", decode_prefix(flipped, TREE_B), ("bdafadcab", ""))
    same("letters decoded", len(decode_prefix(flipped, TREE_B)[0]), 9)

    print("\nLooking Ahead")
    check("h(0.01)", binary_entropy(0.01), "0.0808")
    check("1 bit against h(0.01)", 1 / binary_entropy(0.01), "12.4")
    check("p = 0.01, pairs, bits per letter", block_rate(0.01, 2), "0.5150")
    check("p = 0.01, triples, bits per letter", block_rate(0.01, 3), "0.3533")

    print("\nWorked Examples")
    five = [Fraction(n, 100) for n in (42, 15, 15, 15, 13)]
    check("five-letter source entropy", entropy(five), "2.1399")
    len5, merges5 = huffman(five)
    same("five-letter Huffman lengths", len5, [1, 3, 3, 3, 3])
    same("five-letter merges", [float(m) for m in merges5], [0.28, 0.3, 0.58, 1.0])
    check("five-letter Huffman average", average_length(five, len5), "2.16")
    f5, tied5 = fano_code(five)
    same("five-letter Fano codewords", f5, ["00", "01", "10", "110", "111"])
    require("no Fano split of the five-letter source is tied", not tied5)
    check("five-letter Fano average", average_length(five, [len(w) for w in f5]), "2.28")
    check("Fano first split gap with a alone, 0.58 - 0.42", Fraction(58, 100) - Fraction(42, 100), "0.16")
    check("Fano first split gap with a and b, 0.57 - 0.43", Fraction(57, 100) - Fraction(43, 100), "0.14")
    check("five-letter, Fano minus Huffman", Fraction(228, 100) - Fraction(216, 100), "0.12")

    lean = [Fraction(n, 100) for n in (40, 30, 25, 5)]
    check("four-letter lean source entropy", entropy(lean), "1.7660")
    len_lean, merges_lean = huffman(lean)
    same("lean source Huffman lengths", len_lean, [1, 2, 3, 3])
    same("lean source merges", [float(m) for m in merges_lean], [0.3, 0.6, 1.0])
    check("lean source Huffman average", average_length(lean, len_lean), "1.9")
    same("lean source Shannon lengths", shannon_lengths(lean), [2, 2, 2, 5])
    check("lean source Shannon average", average_length(lean, shannon_lengths(lean)), "2.15")
    check("surprise of 0.25", log2(4), "2")
    check("surprise of 0.05", log2(20), "4.322")

    steep = [Fraction(n, 100) for n in (40, 25, 15, 10, 5, 5)]
    len_steep, merges_steep = huffman(steep)
    same("steep source Huffman lengths", len_steep, [1, 2, 3, 4, 5, 5])
    same("steep source merges", [float(m) for m in merges_steep], [0.1, 0.2, 0.35, 0.6, 1.0])
    check("steep source Huffman average", average_length(steep, len_steep), "2.25")
    check("steep source entropy", entropy(steep), "2.2037")
    check("steep source redundancy", 2.25 - entropy(steep), "0.0463")
    check("surprise of 0.05 rounded up", shannon_lengths([Fraction(5, 100)])[0], "5")

    uniform6 = [Fraction(1, 6)] * 6
    check("tree B used on a uniform source, average length", average_length(uniform6, len_b), "3")
    check("log2 6", log2(6), "2.585")
    check("excess over log2 6", 3 - log2(6), "0.415")
    check("D(U || 2^-l) for tree B", relative_entropy(uniform6, [Fraction(1, 2 ** l) for l in len_b]), "0.415")
    len_u, _ = huffman(uniform6)
    same("Huffman lengths for the uniform source", sorted(len_u), [2, 2, 3, 3, 3, 3])
    check("Huffman average for the uniform source", average_length(uniform6, len_u), "2.667")
    check("tree B on the six-letter source, cross-check", average_length(SIX, len_b), "2.25")

    print("\nExercises")
    for words, ud_stated, prefix_stated in ((["1", "01", "001", "000"], True, True),
                                            (["1", "10", "00"], True, False),
                                            (["0", "10", "01"], False, False)):
        ud, _ = sardinas_patterson(words)
        same(f"code {' '.join(words)}: uniquely decodable, prefix", (ud, is_prefix_free(words)), (ud_stated, prefix_stated))
    same("code 1 10 00, readings of 1000", parses("1000", {"x": "1", "y": "10", "z": "00"}), ["yz"])
    same("code 0 10 01, readings of 010", parses("010", {"x": "0", "y": "10", "z": "01"}), ["xy", "zx"])
    check("code 1 01 001 000, Kraft sum", kraft_sum([1, 2, 3, 3]), "1")
    ex_lengths = [1, 2, 4, 4, 4]
    check("lengths 1 2 4 4 4, Kraft sum", kraft_sum(ex_lengths), "0.9375")
    same("Kraft construction for 1 2 4 4 4", kraft_code(ex_lengths), ["0", "10", "1100", "1101", "1110"])
    k = kraft_sum([1, 2, 2, 3])
    n_fail = mcmillan_first_failure([1, 2, 2, 3])
    same("first n with 1.125^n above 2n + 1", n_fail, 37)
    check("1.125^36", k ** 36, "69.42")
    check("1.125^37", k ** 37, "78.10")
    q4 = [Fraction(n, 10) for n in (4, 3, 2, 1)]
    check("source 0.4 0.3 0.2 0.1, entropy", entropy(q4), "1.8464")
    same("source 0.4 0.3 0.2 0.1, Huffman lengths", huffman(q4)[0], [1, 2, 3, 3])
    check("source 0.4 0.3 0.2 0.1, Huffman average", average_length(q4, huffman(q4)[0]), "1.9")
    same("source 0.4 0.3 0.2 0.1, Shannon lengths", shannon_lengths(q4), [2, 2, 3, 4])
    check("source 0.4 0.3 0.2 0.1, Shannon average", average_length(q4, shannon_lengths(q4)), "2.4")
    check("h(0.2)", binary_entropy(0.2), "0.7219")
    pair8 = sorted(block_probs(Fraction(2, 10), 2), reverse=True)
    same("pairs of the 0.8 coin", [float(p) for p in pair8], [0.64, 0.16, 0.16, 0.04])
    same("pairs of the 0.8 coin, merges", [float(m) for m in huffman(pair8)[1]], [0.2, 0.36, 1.0])
    check("pairs of the 0.8 coin, L_2", sum(huffman(pair8)[1]), "1.56")
    check("pairs of the 0.8 coin, per letter", sum(huffman(pair8)[1]) / 2, "0.78")
    uniform5 = [Fraction(1, 5)] * 5
    len_u5, merges_u5 = huffman(uniform5)
    same("five equal letters, Huffman lengths", sorted(len_u5), [2, 2, 2, 3, 3])
    check("five equal letters, Huffman average", average_length(uniform5, len_u5), "2.4")
    check("log2 5", log2(5), "2.3219")
    check("five equal letters, redundancy", 2.4 - log2(5), "0.0781")
    sent = encode("decafd", TREE_B)
    same("decafd encoded with tree B", sent, "101111010001111101")
    received = "001111010001111101"
    same("received stream differs from sent in the first bit only", [i for i in range(18) if sent[i] != received[i]], [0])
    same("received tree-B stream decoded", decode_prefix(received, TREE_B), ("aafacafd", ""))
    same("seventeenth bit flipped instead, decoded with leftover", decode_prefix(sent[:16] + "1" + sent[17:], TREE_B), ("decaf", "111"))
    same("seventeenth bit of the sent stream, before the flip", sent[16], "0")
    same("first bit flipped: the last four letters are the message's", decode_prefix(received, TREE_B)[0][-4:], "decafd"[-4:])

    further_checks()

    print()
    print("All chapter 5 numbers reproduced.")


def codeword_ends(bits: str, words: dict) -> list:
    """Bit positions, counted from 1, at which each decoded codeword ends."""
    inverse = {w: c for c, w in words.items()}
    ends, current = [], ""
    for position, b in enumerate(bits, start=1):
        current += b
        if current in inverse:
            ends.append(position)
            current = ""
    return ends


def split_gaps(probs) -> list:
    """The gap between the two groups for every place a sorted list can be split."""
    total = sum(probs)
    return [abs(2 * sum(probs[:k]) - total) for k in range(1, len(probs))]


def further_checks() -> None:
    """Numbers stated in the prose and the figure texts that the sections above do not reach."""
    print()
    print("Further checks on the text")

    same("code 0 1 01 10, Kraft sum", kraft_sum(lengths_of(CODE_SHORT)), Fraction(3, 2))
    len_four, merges_four = huffman(FOUR)
    same("Huffman merges on the four-letter source", merges_four, [Fraction(1, 4), Fraction(1, 2), Fraction(1)])
    check("four-letter merges add to", sum(merges_four), "1.75")
    same("Huffman lengths on the four-letter source are code A's", len_four, lengths_of(CODE_PREFIX))

    same("interval named by the string 10", (Fraction(0b10, 4), Fraction(0b10 + 1, 4)), (Fraction(1, 2), Fraction(3, 4)))
    same("lengths 1 3 3 4, Kraft sum as a fraction", kraft_sum([1, 3, 3, 4]), Fraction(13, 16))
    same("lengths 1 3 3 4, unused fraction of the line", 1 - kraft_sum([1, 3, 3, 4]), Fraction(3, 16))
    shortened = ["0", "100", "101", "11"]
    same("shortened code 0 100 101 11, Kraft sum", kraft_sum([len(w) for w in shortened]), Fraction(1))
    require("shortened code 0 100 101 11 is a prefix code", is_prefix_free(shortened))

    pair_lengths = [len(CODE_PREFIX[x]) + len(CODE_PREFIX[y]) for x in "abcd" for y in "abcd"]
    same("code A pairs, A_m for m = 2 to 6", [pair_lengths.count(m) for m in range(2, 7)], [1, 2, 5, 4, 4])

    shannon_k = kraft_sum(shannon_lengths(SIX))
    same("Shannon code, Kraft sum as 1/4 + 1/8 + 3/16 + 1/32",
         Fraction(1, 4) + Fraction(1, 8) + Fraction(3, 16) + Fraction(1, 32), shannon_k)
    same("Shannon code, unused fraction of the tree", 1 - shannon_k, Fraction(13, 32))
    check("Shannon code, unused fraction in percent", 100 * float(1 - shannon_k), "41")
    same("binary expansion of 0.75, four bits", binary_digits(Fraction(3, 4), 4), "1100")
    fano_words, _ = fano_code(SIX)
    check("Shannon average minus Fano average",
          average_length(SIX, shannon_lengths(SIX)) - average_length(SIX, [len(w) for w in fano_words]), "0.70")
    first_gaps = split_gaps(SIX)
    check("Fano first split, gap with a alone", first_gaps[0], "0.10")
    check("Fano first split, gap with a and b", first_gaps[1], "0.30")
    rest_gaps = split_gaps(SIX[1:])
    same("Fano split of b to f, gaps after b, c, d, e", [float(g) for g in rest_gaps], [0.15, 0.05, 0.25, 0.45])
    def_gaps = split_gaps(SIX[3:])
    same("Fano split of d e f, gaps after d and after e", [float(g) for g in def_gaps], [0.05, 0.15])
    check("tree A, 0.45 + 2(0.20) + 4(0.35)", average_length(SIX, [len(TREE_A[c]) for c in LETTERS6]), "2.25")

    five = [Fraction(n, 100) for n in (42, 15, 15, 15, 13)]
    cde_gaps = split_gaps(five[2:])
    check("five-letter Fano split of c d e, gap with c alone", cde_gaps[0], "0.13")
    check("five-letter Fano split of c d e, gap splitting off e alone", cde_gaps[1], "0.17")
    lean = [Fraction(n, 100) for n in (40, 30, 25, 5)]
    same("lean source Shannon lengths, Kraft sum", kraft_sum(shannon_lengths(lean)), Fraction(25, 32))
    uniform6 = [Fraction(1, 6)] * 6
    check("tree B on the uniform source minus the uniform Huffman code",
          average_length(uniform6, lengths_of(TREE_B)) - average_length(uniform6, huffman(uniform6)[0]), "0.3333")

    h01 = binary_entropy(0.1)
    for n, stated in zip(range(5, 9), ("0.6690", "0.6357", "0.6119", "0.5940")):
        check(f"H + 1/{n}", h01 + 1 / n, stated)
    require("every rate from n = 1 to 8 lies in [h(0.1), h(0.1) + 1/n)",
            all(h01 <= block_rate(0.1, n) < h01 + 1 / n for n in range(1, 9)))
    for count, expected in ((5, (3, 2)), (6, (2, 4)), (36, (28, 8)), (216, (40, 176))):
        k = count.bit_length() - 1
        lengths_n, _ = huffman([Fraction(1, count)] * count)
        same(f"{count} equal letters, Huffman leaves at depths {k} and {k + 1}",
             (lengths_n.count(k), lengths_n.count(k + 1)), expected)
        same(f"{count} equal letters, 2^(k+1) - N and 2(N - 2^k)", (2 ** (k + 1) - count, 2 * (count - 2 ** k)), expected)

    sent = encode("fadedcab", TREE_B)
    received = sent[:2] + ("0" if sent[2] == "1" else "1") + sent[3:]
    sent_ends, received_ends = codeword_ends(sent, TREE_B), codeword_ends(received, TREE_B)
    common = sorted(set(sent_ends) & set(received_ends))
    same("fadedcab, first codeword end shared by sent and received", common[0], 12)
    require("fadedcab, every end from bit 12 on is shared", [e for e in sent_ends if e >= 12] == [e for e in received_ends if e >= 12])
    decoded = decode_prefix(received, TREE_B)[0]
    require("fadedcab, the first five decoded letters are all wrong", all(a != b for a, b in zip(decoded[:5], "fadedcab"[:5])))
    same("fadedcab, the last four decoded letters, one place late", (decoded[5:], "fadedcab"[4:]), ("dcab", "dcab"))

    same("code 1 10 00 written backward", [w[::-1] for w in ("1", "10", "00")], ["1", "01", "00"])
    require("code 1 10 00 written backward is a prefix code", is_prefix_free(["1", "01", "00"]))
    same("starting points for lengths 1 2 4 4 4, in sixteenths",
         [int(sum(Fraction(1, 2 ** l) for l in [1, 2, 4, 4, 4][:i]) * 16) for i in range(5)], [0, 8, 12, 13, 14])
    same("lengths 1 2 4 4 4 with 1110 shortened to 111, Kraft sum", kraft_sum([1, 2, 4, 4, 3]), Fraction(1))
    require("codewords 0 10 1100 1101 111 form a prefix code", is_prefix_free(["0", "10", "1100", "1101", "111"]))
    dyadic = [Fraction(1, 4)] * 3 + [Fraction(1, 8)] * 2
    dyadic_lengths = shannon_lengths(dyadic)
    same("dyadic source, lengths log2(1/p)", dyadic_lengths, [2, 2, 2, 3, 3])
    same("dyadic source, Kraft construction", kraft_code(dyadic_lengths), ["00", "01", "10", "110", "111"])
    check("dyadic source, average length", average_length(dyadic, dyadic_lengths), "2.25")
    check("dyadic source, entropy", entropy(dyadic), "2.25")


if __name__ == "__main__":
    main()
