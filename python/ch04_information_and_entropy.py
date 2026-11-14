"""Reproduce every number in chapter 4 of Error Correction Coding, "Information and Entropy".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch04_information_and_entropy.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees by
more than half a unit in the last digit the chapter prints. A run that reaches
the final line has reproduced every number in the chapter.

The module also holds the entropy and mutual-information functions that
ch04_charts.py plots, so the charts and this check share one source.
"""

import sys
from collections import Counter
from math import e, log, log2


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


def require(label: str, condition: bool) -> None:
    print(f"  {label:<60} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH in '{label}'")
        sys.exit(1)


def rnd(value: float, places: int) -> float:
    """The value as the chapter prints it, so sums of printed values can be checked too."""
    return round(value, places)


# Information measures --------------------------------------------------------------

P_FLIP = 0.1          # the flip probability of the chapter's binary symmetric channel
RUNNING_Q = 0.25      # P(X = 1) in the running example


def surprise(p: float) -> float:
    """The surprise of an outcome of probability p, in bits."""
    return log2(1 / p)


def entropy(probs) -> float:
    """H = sum of p log2(1/p), with 0 log 1/0 taken as 0."""
    return sum(p * log2(1 / p) for p in probs if p > 0)


def h(p: float) -> float:
    """The binary entropy function."""
    return entropy([p, 1 - p])


def h_slope(p: float) -> float:
    """dh/dp = log2((1 - p) / p)."""
    return log2((1 - p) / p)


def relative_entropy(P, Q) -> float:
    """D(P || Q) in bits; infinite when Q gives zero to an outcome P does not."""
    total = 0.0
    for p, q in zip(P, Q):
        if p > 0:
            if q == 0:
                return float("inf")
            total += p * log2(p / q)
    return total


def bsc_joint(q: float, p: float) -> dict:
    """Joint distribution of (X, Y) when P(X = 1) = q and X crosses BSC(p)."""
    return {(0, 0): (1 - q) * (1 - p), (0, 1): (1 - q) * p, (1, 0): q * p, (1, 1): q * (1 - p)}


def bsc_output_one(q: float, p: float) -> float:
    """P(Y = 1) for the binary symmetric channel, r = q(1 - p) + (1 - q)p."""
    return q * (1 - p) + (1 - q) * p


def bsc_mi(q: float, p: float) -> float:
    """I(X; Y) = h(r) - h(p) for the binary symmetric channel."""
    return h(bsc_output_one(q, p)) - h(p)


def z_mi(q: float, p: float) -> float:
    """I(X; Y) = h(q(1 - p)) - q h(p) for the Z-channel, where only a 1 can flip, to 0."""
    return h(q * (1 - p)) - q * h(p)


def bec_mi(q: float, eps: float) -> float:
    """I(X; Y) = (1 - eps) h(q) for the binary erasure channel."""
    return (1 - eps) * h(q)


def mi_from_joint(joint: dict) -> float:
    """I = sum of p(x, y) log2(p(x, y) / (p(x) p(y))), straight from a joint table."""
    px, py = Counter(), Counter()
    for (x, y), v in joint.items():
        px[x] += v
        py[y] += v
    return sum(v * log2(v / (px[x] * py[y])) for (x, y), v in joint.items() if v > 0)


def dice_sum_counts() -> Counter:
    return Counter(a + b for a in range(1, 7) for b in range(1, 7))


# The four distributions the chapter draws as bars.
DISTRIBUTIONS = (
    ("fair coin", [0.5, 0.5]),
    ("biased coin", [0.9, 0.1]),
    ("fair die", [1 / 6] * 6),
    ("loaded die", [0.5] + [0.1] * 5),
)

LOADED_DIE = [0.5] + [0.1] * 5
DYADIC_DIE = [1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32, 1 / 32]
LOPSIDED_DIE = [0.3, 0.2, 0.2, 0.1, 0.1, 0.1]


def main() -> None:
    print("Measuring Surprise")
    check("fair coin heads", surprise(1 / 2), "1")
    check("one of four equally likely outcomes", surprise(1 / 4), "2")
    check("one face of a fair die", surprise(1 / 6), "2.585")
    check("an ace, probability 1/13", surprise(1 / 13), "3.700")
    check("a heart, probability 1/4", surprise(1 / 4), "2")
    check("the ace of hearts, probability 1/52", surprise(1 / 52), "5.700")
    check("rank plus suit", surprise(1 / 13) + surprise(1 / 4), "5.700")
    check("biased coin, the common outcome 0.9", surprise(0.9), "0.152")
    check("biased coin, the rare outcome 0.1", surprise(0.1), "3.322")
    check("two dice both showing six, 1/36", surprise(1 / 36), "5.170")
    check("die face, rounded, doubled", 2 * rnd(surprise(1 / 6), 3), "5.170")
    check("bits per nat, log2 e", log2(e), "1.4427")
    require("surprise is additive: s(1/52) = s(1/13) + s(1/4)",
            abs(surprise(1 / 52) - surprise(1 / 13) - surprise(1 / 4)) < 1e-12)

    print("\nEntropy as Average Surprise")
    check("biased coin: 0.9 times its surprise", 0.9 * rnd(surprise(0.9), 3), "0.1368")
    check("biased coin: 0.1 times its surprise", 0.1 * rnd(surprise(0.1), 3), "0.3322")
    check("biased coin entropy", entropy([0.9, 0.1]), "0.4690")
    check("four-sided source 1/2, 1/4, 1/8, 1/8", entropy([0.5, 0.25, 0.125, 0.125]), "1.75")
    check("fair die", entropy([1 / 6] * 6), "2.5850")
    check("loaded die", entropy(LOADED_DIE), "2.1610")
    check("loaded die, rounded surprise 3.3219", 0.5 + 5 * 0.1 * 3.3219, "2.1610")
    check("figure: fair coin", entropy([0.5, 0.5]), "1")
    check("figure: biased coin", entropy([0.9, 0.1]), "0.469")
    check("figure: fair die", entropy([1 / 6] * 6), "2.585")
    check("figure: loaded die", entropy(LOADED_DIE), "2.161")
    for p in (1e-2, 1e-4, 1e-6):
        require(f"p log2(1/p) at p = {p:g} is below 2 sqrt(p) / ln 2",
                p * log2(1 / p) < 2 * p**0.5 / log(2))
    check("p log2(1/p) at p = 0.001", 0.001 * log2(1000), "0.00997")

    print("\nThe Binary Entropy Function")
    check("h(0.01)", h(0.01), "0.0808")
    check("h(0.01) / 0.01, about eight", h(0.01) / 0.01, "8.08")
    check("h(0.1)", h(0.1), "0.4690")
    check("h(0.25)", h(0.25), "0.8113")
    check("h(0.5)", h(0.5), "1")
    check("slope at p = 0.01, log2 99", h_slope(0.01), "6.629")
    check("slope at p = 0.1, log2 9", h_slope(0.1), "3.170")
    check("slope at p = 0.25, log2 3", h_slope(0.25), "1.585")
    require("h(p) = h(1 - p) on a grid", all(abs(h(k / 100) - h(1 - k / 100)) < 1e-12 for k in range(1, 100)))
    require("h rises below 1/2 and falls above it",
            all(h(k / 1000) < h((k + 1) / 1000) for k in range(1, 500))
            and all(h(k / 1000) > h((k + 1) / 1000) for k in range(500, 999)))
    step = 1e-6
    require("finite-difference slope matches log2((1-p)/p) at p = 0.2",
            abs((h(0.2 + step) - h(0.2 - step)) / (2 * step) - h_slope(0.2)) < 1e-6)
    check("1 - h(0.1), chapter 1's capacity", 1 - h(0.1), "0.531")

    print("\nConvexity and Jensen's Inequality")
    check("X in {1, 4}: E[log2 X]", (log2(1) + log2(4)) / 2, "1")
    check("X in {1, 4}: log2 E[X]", log2(2.5), "1.3219")
    check("X in {1, 4}: E[X^2]", (1 + 16) / 2, "8.5")
    check("X in {1, 4}: (E[X])^2", 2.5**2, "6.25")
    check("ln 2", log(2), "0.6931")
    check("ln 0.5", log(0.5), "-0.6931")
    require("ln x <= x - 1 on a grid from 0.01 to 10", all(log(k / 100) <= k / 100 - 1 + 1e-15 for k in range(1, 1001)))
    xs, ws = [0.5, 2.0, 7.0], [0.2, 0.5, 0.3]
    require("Jensen for log2 on three points", sum(w * log2(x) for x, w in zip(xs, ws)) <= log2(sum(w * x for x, w in zip(xs, ws))))

    print("\nRelative Entropy and Gibbs' Inequality")
    check("log2(0.9 / 0.5)", log2(0.9 / 0.5), "0.8480")
    check("0.9 log2(0.9 / 0.5)", 0.9 * log2(0.9 / 0.5), "0.7632")
    check("0.1 log2(0.1 / 0.5)", 0.1 * log2(0.1 / 0.5), "-0.2322")
    check("D(biased || fair)", relative_entropy([0.9, 0.1], [0.5, 0.5]), "0.5310")
    check("cross-entropy of the biased coin under the fair model", 0.9 * 1 + 0.1 * 1, "1")
    check("that minus the true entropy", 1 - entropy([0.9, 0.1]), "0.5310")
    check("0.5 log2(0.5 / 0.9)", 0.5 * log2(0.5 / 0.9), "-0.4240")
    check("0.5 log2(0.5 / 0.1)", 0.5 * log2(0.5 / 0.1), "1.1610")
    check("D(fair || biased)", relative_entropy([0.5, 0.5], [0.9, 0.1]), "0.7370")
    require("D is infinite when Q forbids an outcome P allows",
            relative_entropy([0.5, 0.5], [1.0, 0.0]) == float("inf"))

    print("\nHow Large Entropy Can Be")
    check("log2 6 minus the loaded die's entropy", log2(6) - entropy(LOADED_DIE), "0.4240")
    check("D(loaded || uniform)", relative_entropy(LOADED_DIE, [1 / 6] * 6), "0.4240")
    require("D(P || U) = log2|X| - H(P) for the loaded die",
            abs(relative_entropy(LOADED_DIE, [1 / 6] * 6) - (log2(6) - entropy(LOADED_DIE))) < 1e-12)

    print("\nJoint and Conditional Entropy (running example, q = 0.25, BSC(0.1))")
    q, p = RUNNING_Q, P_FLIP
    J = bsc_joint(q, p)
    check("p(0, 0)", J[(0, 0)], "0.675")
    check("p(0, 1)", J[(0, 1)], "0.075")
    check("p(1, 0)", J[(1, 0)], "0.025")
    check("p(1, 1)", J[(1, 1)], "0.225")
    check("P(X = 0)", J[(0, 0)] + J[(0, 1)], "0.75")
    check("P(Y = 0)", J[(0, 0)] + J[(1, 0)], "0.7")
    check("P(Y = 1)", J[(0, 1)] + J[(1, 1)], "0.3")
    terms = {k: v * log2(1 / v) for k, v in J.items()}
    check("term 0.675 log2(1/0.675)", terms[(0, 0)], "0.3828")
    check("term 0.075 log2(1/0.075)", terms[(0, 1)], "0.2803")
    check("term 0.025 log2(1/0.025)", terms[(1, 0)], "0.1330")
    check("term 0.225 log2(1/0.225)", terms[(1, 1)], "0.4842")
    check("sum of the four rounded terms", sum(rnd(t, 4) for t in terms.values()), "1.2803")
    HXY = entropy(J.values())
    check("H(X, Y)", HXY, "1.2803")
    HX = h(q)
    check("0.75 log2(4/3)", 0.75 * log2(4 / 3), "0.3113")
    check("H(X) = h(0.25)", HX, "0.8113")
    HYgX = (1 - q) * h(p) + q * h(p)
    check("H(Y | X = 0) = h(0.1)", h(p), "0.4690")
    check("H(Y | X)", HYgX, "0.4690")
    check("chain rule: H(X) + H(Y | X)", HX + HYgX, "1.2803")
    check("chain rule, rounded parts", rnd(HX, 4) + rnd(HYgX, 4), "1.2803")
    HY = h(0.3)
    check("0.3 log2(10/3)", 0.3 * log2(10 / 3), "0.5211")
    check("0.7 log2(10/7)", 0.7 * log2(10 / 7), "0.3602")
    check("H(Y) = h(0.3)", HY, "0.8813")
    HXgY = HXY - HY
    check("H(X | Y) = H(X, Y) - H(Y)", HXgY, "0.3990")
    check("H(X | Y), rounded parts", rnd(HXY, 4) - rnd(HY, 4), "0.3990")
    post0 = J[(1, 0)] / 0.7
    post1 = J[(1, 1)] / 0.3
    check("P(X = 1 | Y = 0) = 0.025 / 0.7", post0, "0.0357")
    check("P(X = 1 | Y = 1) = 0.225 / 0.3", post1, "0.75")
    check("h(0.0357)", h(post0), "0.2223")
    check("h(0.75)", h(post1), "0.8113")
    check("H(X | Y) = 0.7 h(1/28) + 0.3 h(0.75)", 0.7 * h(post0) + 0.3 * h(post1), "0.3990")
    check("same from the rounded entropies", 0.7 * 0.2223 + 0.3 * 0.8113, "0.3990")
    check("chain rule the other way: H(Y) + H(X | Y)", HY + HXgY, "1.2803")
    require("H(X | Y = 1) equals H(X) exactly here", abs(h(post1) - HX) < 1e-12)
    # Three-variable chain rule on a small check: two fair bits and their XOR.
    xor = {(a, b, a ^ b): 0.25 for a in (0, 1) for b in (0, 1)}
    check("H(X, Y, Z) for two fair bits and their XOR", entropy(xor.values()), "2")

    print("\nConditioning Reduces Entropy")
    require("H(X | Y) <= H(X) in the running example", HXgY <= HX)
    qa = 0.05
    ya1 = bsc_output_one(qa, p)
    pa1 = qa * (1 - p) / ya1
    check("alarm: P(X = 1)", qa, "0.05")
    check("alarm: H(X) = h(0.05)", h(qa), "0.2864")
    check("alarm: P(Y = 1)", ya1, "0.14")
    check("alarm: P(X = 1 | Y = 1) = 0.045 / 0.14", pa1, "0.3214")
    check("alarm: H(X | Y = 1)", h(pa1), "0.9059")
    require("alarm: a received 1 raises the uncertainty", h(pa1) > h(qa))

    print("\nMutual Information")
    I = HX - HXgY
    check("I(X; Y) = H(X) - H(X | Y)", I, "0.4123")
    check("rounded: 0.8113 - 0.3990", rnd(HX, 4) - rnd(HXgY, 4), "0.4123")
    check("I(X; Y) = H(Y) - H(Y | X)", HY - HYgX, "0.4123")
    check("rounded: 0.8813 - 0.4690", rnd(HY, 4) - rnd(HYgX, 4), "0.4123")
    check("I(X; Y) = H(X) + H(Y) - H(X, Y)", HX + HY - HXY, "0.4123")
    check("rounded: 0.8113 + 0.8813 - 1.2803", rnd(HX, 4) + rnd(HY, 4) - rnd(HXY, 4), "0.4123")
    check("I from the joint table directly", mi_from_joint(J), "0.4123")
    check("share of H(X) that crosses, I / H(X)", I / HX, "0.508")
    indep = {(x, y): px * py for x, px in ((0, 0.75), (1, 0.25)) for y, py in ((0, 0.6), (1, 0.4))}
    check("independent pair: I", mi_from_joint(indep), "0")
    same = {(0, 0): 0.75, (1, 1): 0.25}
    check("Y = X: I(X; X) = H(X)", mi_from_joint(same), "0.8113")

    print("\nThe Information Diagram")
    check("left crescent H(X | Y)", HXgY, "0.3990")
    check("overlap I(X; Y)", I, "0.4123")
    check("right crescent H(Y | X)", HYgX, "0.4690")
    check("union H(X, Y)", HXgY + I + HYgX, "1.2803")
    check("XOR: I(X; Y), X and Y independent", 0, "0")
    # I(X; Y | Z) = H(X | Z) - H(X | Y, Z) for the XOR triple.
    hx_z = entropy(xor.values()) - 1          # H(X, Z) = 2 and H(Z) = 1
    check("XOR: H(X | Z)", hx_z, "1")
    check("XOR: H(X | Y, Z)", 0, "0")
    check("XOR: I(X; Y | Z)", hx_z - 0, "1")
    check("XOR: the centre region I(X; Y) - I(X; Y | Z)", 0 - (hx_z - 0), "-1")

    print("\nMutual Information of the Binary Symmetric Channel")
    require("r - 1/2 = (1 - 2p)(q - 1/2) on a grid",
            all(abs(bsc_output_one(a / 20, b / 20) - 0.5 - (1 - 2 * b / 20) * (a / 20 - 0.5)) < 1e-12
                for a in range(21) for b in range(21)))
    check("uniform input: 1 - h(0.1)", bsc_mi(0.5, p), "0.5310")
    check("uniform input: P(Y = 1)", bsc_output_one(0.5, p), "0.5")
    rows = ((0.5, "0.5", "1", "0.5310"), (0.4, "0.42", "0.9815", "0.5125"), (0.3, "0.34", "0.9248", "0.4558"),
            (0.25, "0.3", "0.8813", "0.4123"), (0.1, "0.18", "0.6801", "0.2111"), (0.05, "0.14", "0.5842", "0.1152"))
    for qq, r_s, hy_s, i_s in rows:
        r = bsc_output_one(qq, p)
        check(f"table q = {qq}: P(Y = 1)", r, r_s)
        check(f"table q = {qq}: H(Y)", h(r), hy_s)
        check(f"table q = {qq}: I", bsc_mi(qq, p), i_s)
        check(f"table q = {qq}: rounded H(Y) - 0.4690", rnd(h(r), 4) - 0.4690, i_s)
    require("uniform input is best on a fine grid for BSC(0.1)",
            max(range(1, 1000), key=lambda k: bsc_mi(k / 1000, p)) == 500)
    require("uniform input: I(1 - p) = I(p)", all(abs(bsc_mi(0.5, k / 100) - bsc_mi(0.5, 1 - k / 100)) < 1e-12
                                                  for k in range(1, 100)))
    check("BSC(0.5): I", bsc_mi(0.5, 0.5), "0")
    check("figure: q = 0.25 at p = 0 carries all of H(X)", h(0.25), "0.811")
    check("BEC(0.1), uniform input", bec_mi(0.5, 0.1), "0.9")
    check("Z(0.1), uniform input: h(0.45)", h(0.45), "0.9928")
    check("Z(0.1), uniform input: 0.5 h(0.1)", 0.5 * h(0.1), "0.2345")
    check("Z(0.1), uniform input: I", z_mi(0.5, p), "0.7583")

    print("\nLooking Ahead")
    check("four-sided source: average questions", 0.5 * 1 + 0.25 * 2 + 0.125 * 3 + 0.125 * 3, "1.75")

    print("\nWorked Examples")
    counts = dice_sum_counts()
    probs = [counts[s] / 36 for s in range(2, 13)]
    check("dice sum: log2 36", log2(36), "5.1699")
    check("dice sum: log2 18", log2(18), "4.1699")
    check("dice sum: log2 12", log2(12), "3.5850")
    check("dice sum: log2 9", log2(9), "3.1699")
    check("dice sum: log2 7.2", log2(7.2), "2.8480")
    check("dice sum: log2 6", log2(6), "2.5850")
    HS = entropy(probs)
    check("dice sum: H(S)", HS, "3.2744")
    check("dice sum: from the rounded surprises",
          sum(c / 36 * rnd(log2(36 / c), 4) for c in counts.values()), "3.2744")
    check("dice sum: H(D1, D2)", log2(36), "5.1699")
    check("dice sum: H(D1, D2 | S)", log2(36) - HS, "1.8955")
    check("dice sum: H(S | D1)", log2(6), "2.5850")
    check("dice sum: I(D1; S)", HS - log2(6), "0.6894")
    check("dice sum: rounded I", 3.2744 - 2.5850, "0.6894")
    check("dice sum: rounded H(D1, D2 | S)", 5.1699 - 3.2744, "1.8955")
    check("dice sum: P(S = 7)", counts[7] / 36, "0.1667")

    ia = bsc_mi(qa, p)
    check("alarm: h(0.14)", h(ya1), "0.5842")
    check("alarm: I = h(0.14) - h(0.1)", ia, "0.1152")
    check("alarm: share of H(X) that crosses", ia / h(qa), "0.40")
    pa0 = qa * p / (1 - ya1)
    check("alarm: P(X = 1 | Y = 0) = 0.005 / 0.86", pa0, "0.0058")
    check("alarm: h(0.0058...)", h(pa0), "0.05154")
    check("alarm: h(0.3214...)", h(pa1), "0.90593")
    check("alarm: H(X | Y) from five-place entropies", 0.86 * 0.05154 + 0.14 * 0.90593, "0.1712")
    check("alarm: H(X | Y) exactly", (1 - ya1) * h(pa0) + ya1 * h(pa1), "0.1712")
    check("alarm: H(X) - H(X | Y)", h(qa) - ((1 - ya1) * h(pa0) + ya1 * h(pa1)), "0.1152")
    check("alarm: rounded 0.2864 - 0.1712", 0.2864 - 0.1712, "0.1152")

    check("Z: uniform I", z_mi(0.5, p), "0.7583")
    check("Z: q = 0.4", z_mi(0.4, p), "0.7551")
    best_q = max(range(1, 1000), key=lambda k: z_mi(k / 1000, p)) / 1000
    check("Z: best q on a grid of step 0.001", best_q, "0.456")
    check("Z: I at that q", z_mi(best_q, p), "0.7628")
    check("Z: P(Y = 1) at q = 0.456", 0.456 * 0.9, "0.4104")
    check("Z: gain over the uniform input", z_mi(best_q, p) - z_mi(0.5, p), "0.0046")
    require("Z: the gain is under 0.005 bits", z_mi(best_q, p) - z_mi(0.5, p) < 0.005)
    # The optimum from setting the derivative to zero: (1 - p) log2((1 - r)/r) = h(p), r = (1 - p) q.
    slope_target = h(p) / (1 - p)
    check("Z optimum: h(0.1) / 0.9", slope_target, "0.5211")
    check("Z optimum: 2 to that power", 2**slope_target, "1.4351")
    r_opt = 1 / (1 + 2**slope_target)
    check("Z optimum: 1/r = 1 + 1.4351", 1 / r_opt, "2.4351")
    check("Z optimum: r = 1 / (1 + 1.4351)", r_opt, "0.4107")
    check("Z optimum: same from the rounded 1.4351", 1 / (1 + 1.4351), "0.4107")
    check("Z optimum: q = r / 0.9", r_opt / (1 - p), "0.4563")
    check("Z optimum: same from the rounded 0.4107", 0.4107 / 0.9, "0.4563")
    check("Z optimum: I there", z_mi(r_opt / (1 - p), p), "0.7628")
    require("Z optimum: the derivative vanishes there",
            abs((1 - p) * log2((1 - r_opt) / r_opt) - h(p)) < 1e-12)
    require("Z optimum: the grid's best q is the optimum to three places", round(r_opt / (1 - p), 3) == best_q)
    check("BSC(0.1) uniform, for comparison", bsc_mi(0.5, p), "0.5310")

    check("dyadic die entropy", entropy(DYADIC_DIE), "1.9375")
    require("dyadic die surprises are 1, 2, 3, 4, 5, 5",
            [round(surprise(x), 12) for x in DYADIC_DIE] == [1, 2, 3, 4, 5, 5])
    check("dyadic die: average questions", sum(x * surprise(x) for x in DYADIC_DIE), "1.9375")
    check("dyadic die: log2 6 - H", log2(6) - entropy(DYADIC_DIE), "0.6475")
    check("dyadic die: D(P || U)", relative_entropy(DYADIC_DIE, [1 / 6] * 6), "0.6475")

    print("\nExercises")
    check("face card, probability 3/13: log2(13/3)", surprise(3 / 13), "2.115")
    check("lopsided die entropy", entropy(LOPSIDED_DIE), "2.4464")
    check("lopsided die: below log2 6 by", log2(6) - entropy(LOPSIDED_DIE), "0.1385")
    check("lopsided die: 0.3 log2(10/3), five places", 0.3 * log2(10 / 3), "0.52109")
    check("lopsided die: 0.2 log2 5, five places", 0.2 * log2(5), "0.46439")
    check("lopsided die: 0.1 log2 10, five places", 0.1 * log2(10), "0.33219")
    check("lopsided die: rounded sum 0.52109 + 2(0.46439) + 3(0.33219)", 0.52109 + 2 * 0.46439 + 3 * 0.33219,
          "2.44644")
    check("lopsided die entropy, five places", entropy(LOPSIDED_DIE), "2.44644")
    check("log2 6, five places", log2(6), "2.58496")
    check("lopsided die: gap, five places", log2(6) - entropy(LOPSIDED_DIE), "0.13852")
    check("lopsided die: rounded 2.58496 - 2.44644", 2.58496 - 2.44644, "0.13852")
    check("h(0.2): log2 5", log2(5), "2.3219")
    check("h(0.2): log2 1.25", log2(1.25), "0.3219")
    check("h(0.2)", h(0.2), "0.7219")
    check("h(0.2) from the rounded logarithms", 0.2 * 2.3219 + 0.8 * 0.3219, "0.7219")
    check("X in {1, 5}: E[X^2]", (1 + 25) / 2, "13")
    check("X in {1, 5}: (E[X])^2", 3**2, "9")
    check("X in {1, 5}: variance", (1 + 25) / 2 - 9, "4")
    check("D(fair || biased) again", relative_entropy([0.5, 0.5], [0.9, 0.1]), "0.7370")
    check("D(uniform || loaded die)", relative_entropy([1 / 6] * 6, LOADED_DIE), "0.3500")
    check("BSC(0.01), uniform input: 1 - h(0.01)", bsc_mi(0.5, 0.01), "0.9192")
    check("BEC(0.1), q = 0.25", bec_mi(0.25, 0.1), "0.7302")
    check("BEC(0.1), q = 0.25 from the rounded H(X)", 0.9 * 0.8113, "0.7302")
    check("BSC(0.1), q = 0.25", bsc_mi(0.25, 0.1), "0.4123")
    check("ln 6, the fair die in nats", log(6), "1.7918")
    check("ln 6 in bits", log(6) * log2(e), "2.5850")
    check("rounded: 1.7918 times 1.4427", 1.7918 * 1.4427, "2.5850")

    print()
    print("All chapter 4 numbers reproduced.")


if __name__ == "__main__":
    main()
