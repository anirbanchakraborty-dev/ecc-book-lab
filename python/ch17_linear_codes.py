"""Reproduce every number in chapter 17 of Error Correction Coding, "Linear Codes".

    python3 python/ch17_linear_codes.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in this chapter is a whole number, a binary word, a binary matrix
or a rate printed to six places, so every comparison is exact. The general
claims are also checked exhaustively on small cases: closure on every subset
of the eight words of length 3, G H^T = 0 and C = null space of H for all
4,096 blocks P of four rows and three columns, the column view of H x^T on all
128 words of length 7, the half lemma for every length up to 6, and the
automorphisms of the (7,4) Hamming code on all 5,040 orderings of its seven
positions. A run that reaches the final line has reproduced every number in
the chapter.

A word is a list of bits, a matrix is a list of rows, and every sum is taken
modulo 2, so addition is XOR and multiplication is AND. Positions are numbered
from 1 on the left, as in the chapter.
"""

import sys
from itertools import combinations, permutations


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


# Words and matrices over GF(2) ----------------------------------------------------

def parse(text: str) -> list:
    """Rows written as bit strings separated by spaces, such as '10110 01101'."""
    return [[int(b) for b in w] for w in text.split()]


def word(v: list) -> str:
    return "".join(str(b) for b in v)


def show(M: list) -> str:
    return " ".join(word(r) for r in M)


def join_numbers(v: list) -> str:
    return " ".join(str(t) for t in v)


def add(a: list, b: list) -> list:
    return [(x + y) % 2 for x, y in zip(a, b)]


def weight(v: list) -> int:
    return sum(v)


def dot(a: list, b: list) -> int:
    return sum(x * y for x, y in zip(a, b)) % 2


def transpose(A: list) -> list:
    return [list(col) for col in zip(*A)]


def matmul(A: list, B: list) -> list:
    Bt = transpose(B)
    return [[dot(r, c) for c in Bt] for r in A]


def encode(m: list, G: list) -> list:
    """The codeword m G, the sum of the rows of G where m holds a 1."""
    c = [0] * len(G[0])
    for bit, row in zip(m, G):
        if bit:
            c = add(c, row)
    return c


def all_words(n: int) -> list:
    """Every word of n bits, in counting order with the first bit most significant."""
    return [[(i >> (n - 1 - j)) & 1 for j in range(n)] for i in range(2 ** n)]


def span(G: list) -> list:
    """Every codeword m G, sorted, as bit strings."""
    return sorted({word(encode(m, G)) for m in all_words(len(G))})


def null_space(H: list) -> list:
    """Every word x with H x^T = 0, in counting order, as bit strings."""
    return [word(x) for x in all_words(len(H[0])) if all(dot(r, x) == 0 for r in H)]


def syndrome(H: list, x: list) -> str:
    return word([dot(r, x) for r in H])


def rank(A: list) -> int:
    rows = [r[:] for r in A]
    r = 0
    for c in range(len(rows[0])):
        pivot = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = add(rows[i], rows[r])
        r += 1
    return r


def identity(k: int) -> list:
    return [[1 if j == i else 0 for j in range(k)] for i in range(k)]


def systematic(P: list) -> list:
    """G = [I | P]."""
    return [e + p for e, p in zip(identity(len(P)), P)]


def check_matrix(P: list) -> list:
    """H = [P^T | I]."""
    return [q + e for q, e in zip(transpose(P), identity(len(P[0])))]


def weights(code: list) -> list:
    """The weight distribution A_0, ..., A_n of a list of bit strings."""
    n = len(code[0])
    A = [0] * (n + 1)
    for c in code:
        A[c.count("1")] += 1
    return A


def min_weight(code: list) -> int:
    return min(c.count("1") for c in code if "1" in c)


def min_distance(code: list) -> int:
    return min(sum(a != b for a, b in zip(x, y)) for x, y in combinations(code, 2))


def zero_column_sets(H: list, w: int) -> int:
    """How many sets of w columns of H add to the zero column."""
    cols = transpose(H)
    count = 0
    for chosen in combinations(range(len(cols)), w):
        total = [0] * len(H)
        for j in chosen:
            total = add(total, cols[j])
        if not any(total):
            count += 1
    return count


def least_zero_column_set(H: list) -> int:
    for w in range(1, len(H[0]) + 1):
        if zero_column_sets(H, w):
            return w
    return 0


def support(v: list) -> str:
    return join_numbers([j + 1 for j, b in enumerate(v) if b])


def rate6(k: int, n: int) -> str:
    return f"{k / n:.6f}"


def ceil_log2(w: int) -> int:
    d = 0
    while (1 << d) < w:
        d += 1
    return d


def tree_gates(P: list) -> list:
    """Two-input XOR gates for each parity bit: one fewer than its column weight."""
    return [max(weight(col) - 1, 0) for col in transpose(P)]


def permute(x: str, sigma: list) -> str:
    """New position j takes old position sigma(j)."""
    return "".join(x[s - 1] for s in sigma)


# The codes of the chapter ------------------------------------------------------------

C5 = ["00000", "01101", "10110", "11011"]
G5 = parse("10110 01101")
P5 = parse("110 101")
H5 = parse("11100 10010 01001")

P = parse("110 011 101 111")
G = systematic(P)
H = check_matrix(P)
HAMMING = span(G)

P16 = parse("110 011 100 111")
G16 = systematic(P16)
H16 = check_matrix(P16)
CODE16 = span(G16)

SIGMA = [5, 6, 1, 7, 3, 2, 4]


def hamming_weights() -> list:
    return weights(HAMMING)


def code16_weights() -> list:
    return weights(CODE16)


def main() -> None:
    print("Chapter 17, Linear Codes: every number in the chapter, recomputed")

    section("A code closed under XOR")
    same("C5 contains the zero word", "00000" in C5, True)
    sums = [word(add(parse(a)[0], parse(b)[0])) for a in C5 for b in C5]
    same("ordered pairs of C5 words checked", len(sums), 16)
    same("every XOR of two C5 words lies in C5", all(s in C5 for s in sums), True)
    even3 = [w for w in (word(v) for v in all_words(3)) if w.count("1") % 2 == 0]
    odd3 = [w for w in (word(v) for v in all_words(3)) if w.count("1") % 2 == 1]
    same("even-weight words of length 3", " ".join(even3), "000 011 101 110")
    same("even-weight words closed under XOR",
         all(word(add(parse(a)[0], parse(b)[0])) in even3 for a in even3 for b in even3), True)
    same("odd-weight words of length 3", " ".join(odd3), "001 010 100 111")
    same("001 xor 010", word(add([0, 0, 1], [0, 1, 0])), "011")
    same("odd-weight words closed under XOR", "011" in odd3, False)
    closed = []
    subsets = 0
    for mask in range(1, 256):
        S = [i for i in range(8) if (mask >> i) & 1]
        subsets += 1
        if all((a ^ b) in S for a in S for b in S):
            closed.append(S)
    same("nonempty subsets of the 8 words of length 3", subsets, 255)
    same("of those, closed under XOR", len(closed), 16)
    same("every closed subset contains 000", all(0 in S for S in closed), True)
    same("sizes of the closed subsets", join_numbers(sorted(len(S) for S in closed)),
         "1 2 2 2 2 2 2 2 4 4 4 4 4 4 4 8")

    section("The generator matrix")
    same("G5", show(G5), "10110 01101")
    for m, c in (("00", "00000"), ("01", "01101"), ("10", "10110"), ("11", "11011")):
        same(f"message {m} encodes to", word(encode(parse(m)[0], G5)), c)
    same("the four codewords of G5 are C5", span(G5) == C5, True)
    same("c = (m1, m2, m1+m2, m1, m2) for every message",
         all(encode(m, G5) == [m[0], m[1], (m[0] + m[1]) % 2, m[0], m[1]] for m in all_words(2)), True)
    same("rank of G5", rank(G5), 2)
    same("codewords 2^k for k = 2", len(span(G5)), 4)
    same("rate of C5", rate6(2, 5), "0.400000")
    same("transpose of G5", show(transpose(G5)), "10 01 11 10 01")

    section("The parity-check matrix")
    same("H5", show(H5), "11100 10010 01001")
    same("H5 c^T = 000 for every codeword of C5",
         all(syndrome(H5, parse(c)[0]) == "000" for c in C5), True)
    same("G5 H5^T", show(matmul(G5, transpose(H5))), "000 000")
    same("rank of H5", rank(H5), 3)
    same("words x of length 5 with H5 x^T = 0", len(null_space(H5)), 4)
    same("null space of H5 equals C5", null_space(H5) == C5, True)
    same("syndrome of 10111", syndrome(H5, parse("10111")[0]), "001")
    same("syndrome of 10110", syndrome(H5, parse("10110")[0]), "000")

    section("From the generator to the checks")
    same("P5", show(P5), "110 101")
    same("G5 = [I2 | P5]", systematic(P5) == G5, True)
    same("P5^T", show(transpose(P5)), "11 10 01")
    same("[P5^T | I3] = H5", check_matrix(P5) == H5, True)
    same("P of the (7,4) Hamming code", show(P), "110 011 101 111")
    same("G = [I4 | P]", show(G), "1000110 0100011 0010101 0001111")
    same("H = [P^T | I3]", show(H), "1011100 1101010 0111001")
    same("G H^T", show(matmul(G, transpose(H))), "000 000 000 000")
    ok_zero = ok_null = True
    for bits in range(4096):
        Q = [[(bits >> (3 * i + j)) & 1 for j in range(3)] for i in range(4)]
        GQ, HQ = systematic(Q), check_matrix(Q)
        ok_zero &= not any(any(r) for r in matmul(GQ, transpose(HQ)))
        ok_null &= span(GQ) == null_space(HQ)
    require("G H^T = 0 for all 4096 blocks P of 4 rows, 3 columns", ok_zero)
    require("row space of G = null space of H for all 4096 blocks P", ok_null)
    same("rank of H", rank(H), 3)
    same("null space of H has 2^4 words", len(null_space(H)), 16)
    same("null space of H equals the row space of G", null_space(H) == HAMMING, True)

    section("Systematic encoding with the (7,4) code")
    for j, stated in enumerate(("1 3 4", "1 2 4", "2 3 4")):
        same(f"message bits feeding p{j + 1}", support(transpose(P)[j]), stated)
    same("message 1011 encodes to", word(encode(parse("1011")[0], G)), "1011100")
    same("message 0110 encodes to", word(encode(parse("0110")[0], G)), "0110110")
    same("message 1110 encodes to", word(encode(parse("1110")[0], G)), "1110000")
    same("first four bits of every codeword are its message",
         all(word(encode(m, G))[:4] == word(m) for m in all_words(4)), True)
    listed = [word(encode(m, G)) for m in all_words(4)]
    same("the 16 codewords in message order", " ".join(listed),
         "0000000 0001111 0010101 0011010 0100011 0101100 0110110 0111001 "
         "1000110 1001001 1010011 1011100 1100101 1101010 1110000 1111111")
    same("their weights", join_numbers([c.count("1") for c in listed]),
         "0 4 3 3 3 3 4 4 3 3 4 4 4 4 3 7")

    section("The encoder as XOR trees")
    same("column weights of P", join_numbers([weight(c) for c in transpose(P)]), "3 3 3")
    same("two-input gates per parity bit", join_numbers(tree_gates(P)), "2 2 2")
    same("gates in the encoder, one tree per parity bit", sum(tree_gates(P)), 6)
    same("gate depth of a three-input tree", ceil_log2(3), 2)
    same("gates in the C5 encoder", sum(tree_gates(P5)), 1)
    shared_ok = True
    for m1, m2, m3, m4 in all_words(4):
        t = (m1 + m4) % 2
        u = (m2 + m3) % 2
        shared = [(t + m3) % 2, (t + m2) % 2, (u + m4) % 2]
        shared_ok &= shared == encode([m1, m2, m3, m4], G)[4:]
    require("five shared gates give p1, p2, p3 for all 16 messages", shared_ok)

    section("Distance is the smallest weight")
    same("weight distribution of the (7,4) code", join_numbers(weights(HAMMING)), "1 0 0 7 7 0 0 1")
    same("nonzero weights checked", len(HAMMING) - 1, 15)
    same("pairs of codewords", len(list(combinations(HAMMING, 2))), 120)
    same("minimum nonzero weight", min_weight(HAMMING), 3)
    same("minimum distance over all 120 pairs", min_distance(HAMMING), 3)
    same("distances from each codeword repeat the weights",
         all(sorted(sum(a != b for a, b in zip(c, x)) for x in HAMMING)
             == sorted(x.count("1") for x in HAMMING) for c in HAMMING), True)
    same("minimum weight of C5", min_weight(C5), 3)
    same("minimum distance of C5", min_distance(C5), 3)
    same("G of the code from Matrices and Elimination", show(G16), "1000110 0100011 0010100 0001111")
    same("weight distribution of that code", join_numbers(weights(CODE16)), "1 0 1 6 5 2 1 0")
    same("its minimum weight", min_weight(CODE16), 2)
    same("0010100 is one of its codewords", "0010100" in CODE16, True)
    same("message 1011 encodes to", word(encode(parse("1011")[0], G16)), "1011101")

    section("Reading the distance from the columns of H")
    same("columns of H", show(transpose(H)), "110 011 101 111 100 010 001")
    cols = [word(c) for c in transpose(H)]
    same("every column nonzero", all("1" in c for c in cols), True)
    same("all seven columns different", len(set(cols)), 7)
    same("nonzero words of length 3", 2 ** 3 - 1, 7)
    same("columns 1 + 5 + 6", word(add(add(transpose(H)[0], transpose(H)[4]), transpose(H)[5])), "000")
    same("1000110 is a codeword", "1000110" in HAMMING, True)
    column_view = True
    for x in all_words(7):
        total = [0, 0, 0]
        for j in range(7):
            if x[j]:
                total = add(total, transpose(H)[j])
        column_view &= word(total) == syndrome(H, x)
    require("H x^T = sum of the columns at the 1s of x, all 128 x", column_view)
    same("sets of w columns adding to zero, w = 1 to 7",
         join_numbers([zero_column_sets(H, w) for w in range(1, 8)]), "0 0 7 7 0 0 1")
    same("codewords of weight w, w = 1 to 7", join_numbers(weights(HAMMING)[1:]), "0 0 7 7 0 0 1")
    same("fewest columns of H adding to zero", least_zero_column_set(H), 3)
    same("columns of H5", show(transpose(H5)), "110 101 100 010 001")
    same("columns 1 + 3 + 4 of H5", word(add(add(transpose(H5)[0], transpose(H5)[2]), transpose(H5)[3])), "000")
    same("fewest columns of H5 adding to zero", least_zero_column_set(H5), 3)
    same("columns of H for the Matrices and Elimination code", show(transpose(H16)),
         "110 011 100 111 100 010 001")
    same("its columns 3 and 5 are equal", transpose(H16)[2] == transpose(H16)[4], True)
    same("fewest of its columns adding to zero", least_zero_column_set(H16), 2)

    section("The Venn diagram of the (7,4) code")
    for i, stated in enumerate(("1 3 4 5", "1 2 4 6", "2 3 4 7")):
        same(f"positions inside circle {i + 1}", support(H[i]), stated)
    same("circles holding positions 1 to 7",
         " ".join("".join(str(i + 1) for i in range(3) if H[i][j]) for j in range(7)),
         "12 23 13 123 1 2 3")
    same("circle parities of 1011100", syndrome(H, parse("1011100")[0]), "000")
    same("every codeword has three even circles",
         all(syndrome(H, parse(c)[0]) == "000" for c in HAMMING), True)
    same("words of length 7 with three even circles", len(null_space(H)), 16)
    same("every position lies in some circle", all(any(H[i][j] for i in range(3)) for j in range(7)), True)

    section("The dual code")
    dual5 = span(H5)
    same("C5 dual, the row space of H5", " ".join(dual5),
         "00000 00111 01001 01110 10010 10101 11011 11100")
    same("C5 dual = words orthogonal to all of C5",
         dual5 == [word(y) for y in all_words(5) if all(dot(y, parse(c)[0]) == 0 for c in C5)], True)
    same("dimension of C5 + dimension of its dual", f"{rank(G5)} + {rank(H5)} = {rank(G5) + rank(H5)}", "2 + 3 = 5")
    same("weight distribution of C5 dual", join_numbers(weights(dual5)), "1 0 2 4 1 0")
    simplex = span(H)
    same("dual of the (7,4) code, the row space of H", " ".join(simplex),
         "0000000 0001111 0110110 0111001 1010011 1011100 1100101 1101010")
    same("its weight distribution", join_numbers(weights(simplex)), "1 0 0 0 7 0 0 0")
    same("H H^T", show(matmul(H, transpose(H))), "000 000 000")
    same("H G^T", show(matmul(H, transpose(G))), "0000 0000 0000")
    same("the dual lies inside the (7,4) code", all(y in HAMMING for y in simplex), True)
    same("null space of G equals the dual", null_space(G) == simplex, True)
    same("dimension of the code + dimension of its dual", f"{rank(G)} + {rank(H)} = {rank(G) + rank(H)}", "4 + 3 = 7")
    same("rate of the dual", rate6(3, 7), "0.428571")
    half = True
    for m in range(1, 7):
        for y in all_words(m)[1:]:
            half &= sum(dot(y, h) for h in all_words(m)) == 2 ** (m - 1)
    same("h in GF(2)^3 with y . h = 1, each nonzero y",
         join_numbers([sum(dot(y, h) for h in all_words(3)) for y in all_words(3)[1:]]), "4 4 4 4 4 4 4")
    require("exactly half of GF(2)^m meets y . h = 1, m = 1 to 6", half)

    section("Equivalent codes")
    same("the permutation sigma", join_numbers(SIGMA), "5 6 1 7 3 2 4")
    Hp = transpose([transpose(H)[s - 1] for s in SIGMA])
    same("H' with columns sigma(1), ..., sigma(7) of H", show(Hp), "1010101 0110011 0001111")
    same("column j of H' read as binary, row 1 lowest",
         join_numbers([c[0] + 2 * c[1] + 4 * c[2] for c in transpose(Hp)]), "1 2 3 4 5 6 7")
    same("1011100 permuted", permute("1011100", SIGMA), "1010101")
    same("1010101 is in the (7,4) code", "1010101" in HAMMING, False)
    same("message 1010 encodes to", word(encode(parse("1010")[0], G)), "1010011")
    Cp = sorted(permute(c, SIGMA) for c in HAMMING)
    same("permuted code is the null space of H'", Cp == null_space(Hp), True)
    same("weight distribution of the permuted code", join_numbers(weights(Cp)), "1 0 0 7 7 0 0 1")
    same("minimum distance of the permuted code", min_distance(Cp), 3)
    fixed = 0
    orderings = 0
    images = set()
    for p in permutations(range(1, 8)):
        orderings += 1
        image = tuple(sorted(permute(c, p) for c in HAMMING))
        images.add(image)
        if list(image) == HAMMING:
            fixed += 1
    same("orderings of 7 positions", orderings, 5040)
    same("orderings that map the (7,4) code to itself", fixed, 168)
    same("different codes equivalent to it", len(images), 30)
    same("5040 / 168", orderings // fixed, len(images))

    section("How many codewords, and at what rate")
    table = (
        ("repetition", "111", "3 1 3", 2, "0.333333"),
        ("even weight, length 3", "101 011", "3 2 2", 4, "0.666667"),
        ("C5", "10110 01101", "5 2 3", 4, "0.400000"),
        ("(6,3) code", "100011 010101 001110", "6 3 3", 8, "0.500000"),
        ("(7,4) Hamming", show(G), "7 4 3", 16, "0.571429"),
        ("its dual", show(H), "7 3 4", 8, "0.428571"),
        ("Matrices and Elimination code", show(G16), "7 4 2", 16, "0.571429"),
    )
    for name, gen, nkd, size, rate in table:
        Gt = parse(gen)
        code = span(Gt)
        n, k = len(Gt[0]), rank(Gt)
        same(f"{name}: n k d", f"{n} {k} {min_weight(code)}", nkd)
        same(f"{name}: codewords", len(code), size)
        same(f"{name}: rate", rate6(k, n), rate)
    same("fraction of 7-bit words in the (7,4) code", f"1/{2 ** 7 // len(HAMMING)}", "1/8")
    same("(15,11) code: codewords", 2 ** 11, 2048)
    same("(15,11) code: rate", rate6(11, 15), "0.733333")
    same("(15,11) code: fraction of words that are codewords", f"1/{2 ** 15 // 2 ** 11}", "1/16")

    section("Worked examples")
    P6 = parse("011 101 110")
    G6, H6 = systematic(P6), check_matrix(P6)
    code6 = span(G6)
    same("(6,3): G", show(G6), "100011 010101 001110")
    same("(6,3): H", show(H6), "011100 101010 110001")
    same("(6,3): columns of H", show(transpose(H6)), "011 101 110 100 010 001")
    same("(6,3): columns distinct and nonzero",
         len({word(c) for c in transpose(H6)}) == 6 and all(any(c) for c in transpose(H6)), True)
    same("(6,3): columns 1 + 5 + 6", word(add(add(transpose(H6)[0], transpose(H6)[4]), transpose(H6)[5])), "000")
    same("(6,3): 100011 is a codeword", "100011" in code6, True)
    same("(6,3): weight distribution", join_numbers(weights(code6)), "1 0 0 4 3 0 0")
    same("(6,3): the eight codewords", " ".join(code6),
         "000000 001110 010101 011011 100011 101101 110110 111000")
    same("(6,3): message 101 encodes to", word(encode(parse("101")[0], G6)), "101101")
    same("(6,3): gates in the encoder", sum(tree_gates(P6)), 3)
    G4 = systematic(parse("1 1 1"))
    H4 = check_matrix(parse("1 1 1"))
    same("even weight, length 4: G", show(G4), "1001 0101 0011")
    same("even weight, length 4: H", show(H4), "1111")
    same("even weight, length 4: codewords", " ".join(span(G4)),
         "0000 0011 0101 0110 1001 1010 1100 1111")
    same("even weight, length 4: minimum distance", min_weight(span(G4)), 2)
    same("even weight, length 4: dual", " ".join(span(H4)), "0000 1111")
    same("even weight, length 4: rate", rate6(3, 4), "0.750000")
    Gn = parse("11011 01101")
    same("non-systematic G spans C5", span(Gn) == C5, True)
    same("row 1 + row 2 of it", word(add(Gn[0], Gn[1])), "10110")
    same("0110: codeword", word(encode(parse("0110")[0], G)), "0110110")
    same("0110: weight of the codeword", word(encode(parse("0110")[0], G)).count("1"), 4)
    same("0110: H c^T", syndrome(H, parse("0110110")[0]), "000")
    m1, m2, m3, m4 = parse("0110")[0]
    t, u = (m1 + m4) % 2, (m2 + m3) % 2
    same("0110: t, p1, p2, u, p3 from the shared gates",
         join_numbers([t, (t + m3) % 2, (t + m2) % 2, u, (u + m4) % 2]), "0 1 1 0 0")

    section("Exercises")
    S = ["0000", "1100", "0011", "1111", "1010"]
    same("1100 xor 1010", word(add(parse("1100")[0], parse("1010")[0])), "0110")
    same("0110 is in the five-word set", "0110" in S, False)
    P53 = parse("11 10 01")
    G53, H53 = systematic(P53), check_matrix(P53)
    same("(5,3): G", show(G53), "10011 01010 00101")
    same("(5,3): H", show(H53), "11010 10101")
    same("(5,3): columns of H", show(transpose(H53)), "11 10 01 10 01")
    same("(5,3): columns 2 and 4 equal", transpose(H53)[1] == transpose(H53)[3], True)
    same("(5,3): 01010 is a codeword", "01010" in span(G53), True)
    same("(5,3): minimum distance", min_weight(span(G53)), 2)
    same("1110: H c^T", syndrome(H, parse("1110000")[0]), "000")
    same("dot products of the rows of H, 11 12 13 22 23 33",
         join_numbers([sum(a * b for a, b in zip(H[i], H[j])) for i, j in
                       ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))]), "4 2 2 4 2 4")
    even = [c for c in HAMMING if c.count("1") % 2 == 0]
    same("even-weight codewords of the (7,4) code", len(even), 8)
    same("they are closed under XOR",
         all(word(add(parse(a)[0], parse(b)[0])) in even for a in even for b in even), True)
    same("they are exactly the dual code", sorted(even) == span(H), True)

    print("\nAll chapter 17 numbers reproduced.")


if __name__ == "__main__":
    main()
