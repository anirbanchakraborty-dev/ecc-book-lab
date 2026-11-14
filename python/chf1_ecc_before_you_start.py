"""Reproduce every number in front-matter chapter f1 of Error Correction Coding, "Before You Start".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/chf1_ecc_before_you_start.py

The chapter is an orientation and a notation primer, so its numbers are of two
kinds. The first describes the book itself: the chapters in each Part, the
totals at each level, and the layout of Part I. Those counts are embedded below
as a table copied from the book's table of contents, and the script checks that
they add up to the 187 chapters the book plans. The second kind is the worked
arithmetic of the primer: sums and products, factorials and binomial
coefficients, powers and base-2 logarithms, floors and ceilings, XOR on bits and
bit vectors, small sets, a few probabilities, and the cases the three proofs
illustrate. Each of those is recomputed from scratch, and each theorem the
chapter proves is also checked by brute force over a range of inputs.

Every value computed here is printed beside the value the chapter states. The
script stops with exit status 1 at the first pair that disagrees by more than
half a unit in the last digit the chapter prints, and at the first count, list
or bit pattern that differs at all.
"""

import sys
from itertools import permutations, product
from math import ceil, comb, exp, factorial, floor, log, log2, prod

# The book's table of contents: label, title, level, chapters.
PARTS = [
    ("Front matter", "Before You Start", "front", 2),
    ("Part I", "Errors, Information and Channels", "undergraduate", 12),
    ("Part II", "Linear Block Codes", "undergraduate", 19),
    ("Part III", "Finite Fields and Their Arithmetic", "undergraduate", 15),
    ("Part IV", "Cyclic Codes, CRCs and Burst Errors", "undergraduate", 12),
    ("Part V", "BCH and Reed-Solomon Codes", "undergraduate", 13),
    ("Part VI", "Soft Decisions and Concatenation", "graduate", 9),
    ("Part VII", "Convolutional Codes and Viterbi Decoding", "graduate", 10),
    ("Part VIII", "Turbo Codes and Iterative Decoding", "graduate", 10),
    ("Part IX", "LDPC Codes", "graduate", 15),
    ("Part X", "Polar Codes", "graduate", 8),
    ("Part XI", "ECC in Memory, Storage and Links", "graduate", 14),
    ("Part XII", "Coded Modulation and Wireless", "graduate", 9),
    ("Part XIII", "Advanced Algebraic Coding Theory", "research", 19),
    ("Part XIV", "Codes, Cryptography and Security", "research", 10),
    ("Part XV", "Frontiers", "research", 7),
    ("Reference", "Reference", "back", 3),
]

# Part I as published: the nine teaching chapters, then the case study, the
# project and the lab, numbered 1 to 12.
PART_I = [
    ("why-codes-exist", "teaching"),
    ("probability-for-coding", "teaching"),
    ("channels-signals-and-likelihoods", "teaching"),
    ("information-and-entropy", "teaching"),
    ("source-coding-and-prefix-codes", "teaching"),
    ("compression-in-practice", "teaching"),
    ("capacity-and-shannons-theorem", "teaching"),
    ("distance-detection-and-correction", "teaching"),
    ("error-control-strategies-and-performance", "teaching"),
    ("case-study-codes-in-the-wild", "case study"),
    ("project-channel-simulator", "project"),
    ("lab-first-hardware-checker", "lab"),
]


def half_unit(stated: str) -> float:
    text = stated.replace(",", "").rstrip("%")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", "").rstrip("%"))
    if stated.endswith("%"):
        computed = 100 * float(computed)
    print(f"  {label:<64} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<64} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def xor(a: str, b: str) -> str:
    """Bitwise XOR of two bit strings of equal length."""
    assert len(a) == len(b)
    return "".join("1" if x != y else "0" for x, y in zip(a, b))


def parity(bits: str) -> int:
    """The XOR of every bit of the string: 1 when it holds an odd number of 1s."""
    return bits.count("1") % 2


def strings(n: int) -> list:
    """Every n-bit string, in counting order."""
    return ["".join(t) for t in product("01", repeat=n)]


def chapters(level: str) -> int:
    return sum(c for _, _, lv, c in PARTS if lv == level)


def check_book() -> None:
    print("The book")
    same("chapters in the book", sum(c for *_, c in PARTS), 187)
    same("Parts, front and back matter aside", sum(1 for lb, *_ in PARTS if lb.startswith("Part ")), 15)
    same("chapters inside the fifteen Parts", sum(c for lb, _, _, c in PARTS if lb.startswith("Part ")), 182)
    same("chapters of front matter", chapters("front"), 2)
    same("chapters of back matter", chapters("back"), 3)
    same("front and back matter together", chapters("front") + chapters("back"), 5)
    same("undergraduate chapters, Parts I to V", chapters("undergraduate"), 71)
    same("graduate chapters, Parts VI to XII", chapters("graduate"), 75)
    same("research chapters, Parts XIII to XV", chapters("research"), 36)
    same("undergraduate Parts", [lb for lb, _, lv, _ in PARTS if lv == "undergraduate"],
         ["Part I", "Part II", "Part III", "Part IV", "Part V"])
    same("graduate Parts", [lb for lb, _, lv, _ in PARTS if lv == "graduate"],
         ["Part VI", "Part VII", "Part VIII", "Part IX", "Part X", "Part XI", "Part XII"])
    same("research Parts", [lb for lb, _, lv, _ in PARTS if lv == "research"],
         ["Part XIII", "Part XIV", "Part XV"])
    same("chapters in each section, front matter to Reference", [c for *_, c in PARTS],
         [2, 12, 19, 15, 12, 13, 9, 10, 10, 15, 8, 14, 9, 19, 10, 7, 3])
    parts = [(lb, c) for lb, _, _, c in PARTS if lb.startswith("Part ")]
    most = max(c for _, c in parts)
    same("most chapters in one Part", most, 19)
    same("the Parts with that many", [lb for lb, c in parts if c == most], ["Part II", "Part XIII"])
    fewest = min(c for _, c in parts)
    same("fewest chapters in one Part", fewest, 7)
    same("the Part with that few", [lb for lb, c in parts if c == fewest], ["Part XV"])
    same("chapters in Part I", len(PART_I), dict((lb, c) for lb, _, _, c in PARTS)["Part I"])
    same("teaching chapters in Part I", sum(1 for _, k in PART_I if k == "teaching"), 9)
    same("Part I closing chapters, in order", [k for _, k in PART_I[9:]], ["case study", "project", "lab"])
    same("chapter number of the case study", PART_I.index(("case-study-codes-in-the-wild", "case study")) + 1, 10)
    same("chapter number of the project", PART_I.index(("project-channel-simulator", "project")) + 1, 11)
    same("chapter number of the lab", PART_I.index(("lab-first-hardware-checker", "lab")) + 1, 12)
    same("reading order of this chapter", 1, 1)
    same("reading order of chapter 1, after the two front chapters", chapters("front") + 1, 3)


def check_sums_and_products() -> None:
    print("Sums and products")
    same("1 + 2 + 3 + 4 + 5", sum(range(1, 6)), 15)
    same("1 x 2 x 3 x 4", prod(range(1, 5)), 24)
    same("a sum with no terms", sum([]), 0)
    same("a product with no factors", prod([]), 1)
    p = 0.01
    check("all 8 bits survive at p = 0.01, (0.99)^8", (1 - p) ** 8, "0.9227")
    check("at least one of 8 bits flips at p = 0.01", 1 - (1 - p) ** 8, "0.0773")
    check("the same as a product of eight factors of 0.99", prod([1 - p] * 8), "0.9227")


def check_binomials() -> None:
    print("Factorials and binomial coefficients")
    same("orders of three letters, 3!", factorial(3), 6)
    same("the orders listed", sorted("".join(t) for t in permutations("abc")),
         ["abc", "acb", "bac", "bca", "cab", "cba"])
    same("5!", factorial(5), 120)
    same("0!", factorial(0), 1)
    two_of_four = [s for s in strings(4) if s.count("1") == 2]
    same("4-bit strings with exactly two 1s", sorted(two_of_four, reverse=True),
         ["1100", "1010", "1001", "0110", "0101", "0011"])
    same("their count, 4 choose 2", comb(4, 2), 6)
    same("4! / (2! 2!)", factorial(4) // (factorial(2) * factorial(2)), 6)
    same("4! and 2! 2!", (factorial(4), factorial(2) * factorial(2)), (24, 4))
    same("ordered choices of 2 positions from 4, 4 x 3", 4 * 3, 12)
    same("each unordered pair counted 2! times, 12 / 2", 12 // factorial(2), 6)
    same("ways exactly 2 of 8 bits flip, 8 choose 2", comb(8, 2), 28)
    same("8 x 7 / 2", 8 * 7 // 2, 28)
    same("Pascal: 7 choose 3", comb(7, 3), 35)
    same("Pascal: 6 choose 2 + 6 choose 3", comb(6, 2) + comb(6, 3), 35)
    same("6 choose 2 and 6 choose 3", (comb(6, 2), comb(6, 3)), (15, 20))
    row7 = [comb(7, i) for i in range(8)]
    same("row 7 of the binomial coefficients", row7, [1, 7, 21, 35, 35, 21, 7, 1])
    same("its sum", sum(row7), 128)
    same("2^7", 2 ** 7, 128)
    same("row 7 reads the same backward", row7 == row7[::-1], True)
    same("Pascal's rule holds for every 1 <= k <= n - 1, n <= 40",
         all(comb(n, k) == comb(n - 1, k - 1) + comb(n - 1, k) for n in range(2, 41) for k in range(1, n)), True)
    same("formula n!/(k!(n-k)!) matches a direct count, n <= 12",
         all(sum(1 for s in strings(n) if s.count("1") == k) == factorial(n) // (factorial(k) * factorial(n - k))
             for n in range(13) for k in range(n + 1)), True)
    same("ordered picks of 2 of 4 positions, 4!/2!", factorial(4) // factorial(2), 12)
    same("n(n-1)...(n-k+1) = n!/(n-k)!, n up to 12", all(prod(range(n - k + 1, n + 1)) == factorial(n) // factorial(n - k) for n in range(1, 13) for k in range(n + 1)), True)
    same("first rows of Pascal's triangle", [[comb(n, k) for k in range(n + 1)] for n in range(4)], [[1], [1, 1], [1, 2, 1], [1, 3, 3, 1]])


def check_powers_and_logs() -> None:
    print("Exponents and base-2 logarithms")
    same("2^10", 2 ** 10, 1024)
    same("values of an 8-bit byte, 2^8", 2 ** 8, 256)
    same("2^3 x 2^4 = 2^7", 2 ** 3 * 2 ** 4, 2 ** 7)
    same("2^7", 2 ** 7, 128)
    same("(2^3)^2 = 2^6", (2 ** 3) ** 2, 64)
    same("2^0", 2 ** 0, 1)
    check("2^-1", 2 ** -1, "0.5")
    check("2^-3", 2 ** -3, "0.125")
    same("log2 1024", log2(1024), 10.0)
    check("log2 1000", log2(1000), "9.966")
    check("log2 10", log2(10), "3.322")
    check("3 x log2 10, the exponent that makes 2 reach 1000", 3 * log2(10), "9.966")
    same("log2 (8 x 32) = log2 8 + log2 32", (log2(8 * 32), log2(8) + log2(32)), (8.0, 8.0))
    same("log2 256", log2(256), 8.0)
    same("seconds per error at 10^-12 on a 10^10 bit/s link", 10 ** 12 // 10 ** 10, 100)
    check("2^(1/2), the square root of 2", 2 ** 0.5, "1.414")


def check_floor_and_ceiling() -> None:
    print("Floor and ceiling")
    same("floor 3.7", floor(3.7), 3)
    same("ceiling 3.2", ceil(3.2), 4)
    same("floor and ceiling of -3.7", (floor(-3.7), ceil(-3.7)), (-4, -3))
    same("cutting the fraction off -3.7", int(-3.7), -3)
    same("floor 4 and ceiling 4", (floor(4), ceil(4)), (4, 4))
    check("100 / 8", 100 / 8, "12.5")
    same("bytes for 100 bits, ceiling of 100/8", ceil(100 / 8), 13)
    same("bits to name 8 things, ceiling of log2 8", ceil(log2(8)), 3)
    same("bits to name 9 things, ceiling of log2 9", ceil(log2(9)), 4)
    same("whole 2s that fit in 5, floor of 5/2", floor(5 / 2), 2)


def check_mod2_and_xor() -> None:
    print("Mod-2 arithmetic and XOR")
    same("7 mod 2", 7 % 2, 1)
    same("13 mod 2", 13 % 2, 1)
    same("1 + 1", 1 + 1, 2)
    same("(1 + 1) mod 2", (1 + 1) % 2, 0)
    same("mod-2 addition table", [[(a + b) % 2 for b in (0, 1)] for a in (0, 1)], [[0, 1], [1, 0]])
    same("XOR truth table", [a ^ b for a in (0, 1) for b in (0, 1)], [0, 1, 1, 0])
    same("XOR agrees with addition mod 2 on every pair",
         all(a ^ b == (a + b) % 2 for a in (0, 1) for b in (0, 1)), True)
    same("mod-2 multiplication equals AND", all(a * b % 2 == (a & b) for a in (0, 1) for b in (0, 1)), True)
    same("1011 XOR 0110", xor("1011", "0110"), "1101")
    same("1101 XOR 0110, the change undone", xor("1101", "0110"), "1011")
    same("parity of 1011", parity("1011"), 1)
    same("cases in the commutativity proof", len(list(product((0, 1), repeat=2))), 4)
    same("cases in the associativity proof", len(list(product((0, 1), repeat=3))), 8)
    vs = strings(4)
    zero = "0000"
    same("commutative on all 4-bit vectors", all(xor(u, v) == xor(v, u) for u in vs for v in vs), True)
    same("associative on all 4-bit vectors",
         all(xor(xor(u, v), w) == xor(u, xor(v, w)) for u in vs for v in vs for w in vs), True)
    same("zero is the identity", all(xor(u, zero) == u for u in vs), True)
    same("every vector is its own inverse", all(xor(u, u) == zero for u in vs), True)
    same("r = m XOR e gives back m = r XOR e",
         all(xor(xor(m, e), e) == m for m in vs for e in vs), True)


def check_sets() -> None:
    print("Sets")
    same("the set {0,1}^3", strings(3), ["000", "001", "010", "011", "100", "101", "110", "111"])
    same("its size", len(strings(3)), 8)
    same("{x in {0,1}^4 : exactly two 1s} has", len([s for s in strings(4) if s.count("1") == 2]), 6)
    same("|{0,1}^n| = 2^n for n <= 12", all(len(strings(n)) == 2 ** n for n in range(13)), True)
    same("{0,1}^1", strings(1), ["0", "1"])


def check_proofs() -> None:
    print("The three kinds of proof")
    same("7 x 8", 7 * 8, 56)
    same("n(n + 1) even for every n from -100 to 100", all(n * (n + 1) % 2 == 0 for n in range(-100, 101)), True)
    same("1 + 2 + 4 + ... + 128", sum(2 ** i for i in range(8)), 255)
    same("2^8 - 1", 2 ** 8 - 1, 255)
    same("11111111 in binary", int("11111111", 2), 255)
    same("sum of 2^i for i < n equals 2^n - 1, n <= 64",
         all(sum(2 ** i for i in range(n)) == 2 ** n - 1 for n in range(65)), True)
    same("induction base, n = 1: 2^0 and 2^1 - 1", (2 ** 0, 2 ** 1 - 1), (1, 1))
    same("row 3 of the binomials", [comb(3, i) for i in range(4)], [1, 3, 3, 1])
    same("its sum, 2^3", sum(comb(3, i) for i in range(4)), 8)
    same("3-bit strings grouped by number of 1s",
         [[s for s in strings(3) if s.count("1") == i] for i in range(4)],
         [["000"], ["001", "010", "100"], ["011", "101", "110"], ["111"]])
    same("rows of binomials sum to 2^n, n <= 60",
         all(sum(comb(n, i) for i in range(n + 1)) == 2 ** n for n in range(61)), True)
    same("1101 in binary", int("1101", 2), 13)
    same("8 + 4 + 0 + 1", 8 + 4 + 0 + 1, 13)


def check_worked_examples() -> None:
    print("Worked examples")
    same("at most two flips in a byte, 1 + 8 + 28", sum(comb(8, i) for i in range(3)), 37)
    same("the three terms", [comb(8, i) for i in range(3)], [1, 8, 28])
    check("share of the 256 patterns", 37 / 256, "14.45%")
    m, e = "10110010", "00010100"
    r = xor(m, e)
    same("10110010 XOR 00010100", r, "10100110")
    same("10100110 XOR 00010100", xor(r, e), m)
    same("1s in the error pattern", e.count("1"), 2)
    same("their positions, counted from the left from 1", [i + 1 for i, b in enumerate(e) if b == "1"], [4, 6])
    same("1s in m and in r", (m.count("1"), r.count("1")), (4, 4))
    same("parity of m and of r", (parity(m), parity(r)), (0, 0))
    same("address bits for 1000 sensors, ceiling of log2 1000", ceil(log2(1000)), 10)
    same("2^9 and 2^10", (2 ** 9, 2 ** 10), (512, 1024))
    same("addresses left over", 2 ** 10 - 1000, 24)
    same("address bits for 1025 sensors", ceil(log2(1025)), 11)
    p, n = 1e-6, 8000
    same("bits in a 1,000-byte packet", 1000 * 8, n)
    check("packet arrives with no flip at p = 10^-6", (1 - p) ** n, "0.9920")
    check("packet suffers at least one flip", 1 - (1 - p) ** n, "0.0080")
    check("damaged packets per 1,000", 1000 * (1 - (1 - p) ** n), "8")
    check("exponent of the approximation, -n p", -n * p, "-0.008")
    check("e^(-np) as the approximation", 2.718281828459045 ** (-n * p), "0.9920")
    bits = "11010110"
    same("11010110 in decimal", sum(int(b) * 2 ** i for i, b in enumerate(reversed(bits))), 214)
    same("128 + 64 + 16 + 4 + 2", 128 + 64 + 16 + 4 + 2, 214)
    same("214 mod 2, the last bit", 214 % 2, 0)
    same("1s in 11010110", bits.count("1"), 5)
    same("parity of 11010110, unlike its last bit", parity(bits), 1)
    same("1011 + 0110 with no carries", xor("1011", "0110"), "1101")
    same("1011 + 0110 in ordinary binary", format(0b1011 + 0b0110, "b"), "10001")
    check("e, the base of the approximation", exp(1), "2.718")
    check("(1 + 1/1000)^1000, on its way to e", (1 + 1 / 1000) ** 1000, "2.717")
    check("e^-1, where (1 - 1/m)^m heads", exp(-1), "0.368")


def check_exercises() -> None:
    print("Exercise answers")
    same("sum of i^2 for i = 1 to 6", sum(i * i for i in range(1, 7)), 91)
    same("1 + 4 + 9 + 16 + 25 + 36", 1 + 4 + 9 + 16 + 25 + 36, 91)
    same("product of 2 for i = 1 to 5", prod([2] * 5), 32)
    same("row 4, the starting row", [comb(4, i) for i in range(5)], [1, 4, 6, 4, 1])
    same("row 5", [comb(5, i) for i in range(6)], [1, 5, 10, 10, 5, 1])
    same("row 6 by Pascal's rule", [1] + [comb(5, k - 1) + comb(5, k) for k in range(1, 6)] + [1],
         [1, 6, 15, 20, 15, 6, 1])
    same("row 6 sum, 2^6", sum(comb(6, i) for i in range(7)), 64)
    same("11001010 XOR 01100110", xor("11001010", "01100110"), "10101100")
    same("10101100 XOR 01100110", xor("10101100", "01100110"), "11001010")
    check("1500 / 8", 1500 / 8, "187.5")
    same("bytes for 1500 bits", ceil(1500 / 8), 188)
    same("XOR tree depth over 64 inputs", ceil(log2(64)), 6)
    same("XOR tree depth over 65 inputs", ceil(log2(65)), 7)
    same("log2 4096", log2(4096), 12.0)
    check("log2 of one million", log2(10 ** 6), "19.93")
    same("bits to label a million things", ceil(log2(10 ** 6)), 20)
    same("2^19 and 2^20", (2 ** 19, 2 ** 20), (524288, 1048576))
    same("n^2 is even exactly when n is, n from -100 to 100",
         all((n * n) % 2 == n % 2 for n in range(-100, 101)), True)
    same("1 + 2 + ... + 100", sum(range(1, 101)), 5050)
    same("100 x 101 / 2", 100 * 101 // 2, 5050)
    same("1 + ... + n = n(n + 1)/2 for n <= 500", all(sum(range(1, n + 1)) == n * (n + 1) // 2 for n in range(501)), True)
    same("8 choose 3 and 8 choose 5", (comb(8, 3), comb(8, 5)), (56, 56))
    same("8 x 7 x 6 / (3 x 2 x 1)", 8 * 7 * 6 // (3 * 2 * 1), 56)
    same("n choose k = n choose n-k for n <= 40",
         all(comb(n, k) == comb(n, n - k) for n in range(41) for k in range(n + 1)), True)
    odd3 = [s for s in strings(3) if parity(s) == 1]
    same("3-bit strings with an odd number of 1s", odd3, ["001", "010", "100", "111"])
    same("their count", len(odd3), 4)
    check("1000 bits with no flip at p = 0.001", 0.999 ** 1000, "0.3677")
    check("e^-1, the approximation", 2.718281828459045 ** -1, "0.3679")
    check("log of 0.999 times 1000", 1000 * log(0.999), "-1.0005")
    check("n p for 1000 bits at 0.001", 1000 * 0.001, "1")
    check("gap between e^-1 and the exact value", 2.718281828459045 ** -1 - 0.999 ** 1000, "0.0002")


def check_looking_ahead() -> None:
    print("Looking ahead")
    same("01100100 in decimal", int("01100100", 2), 100)
    same("its high bit flipped, 11100100", int("11100100", 2), 228)
    same("the flip as an XOR with 10000000", xor("01100100", "10000000"), "11100100")


def main() -> None:
    check_book()
    check_sums_and_products()
    check_binomials()
    check_powers_and_logs()
    check_floor_and_ceiling()
    check_mod2_and_xor()
    check_sets()
    check_proofs()
    check_worked_examples()
    check_exercises()
    check_looking_ahead()
    print("\nAll chapter f1 numbers reproduced.")


if __name__ == "__main__":
    main()
