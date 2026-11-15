"""Write the five computed figures of chapter 23 of Error Correction Coding as SVG files.

    python3 python/ch23_charts.py --out figures

Every plotted value comes from the functions in ch23_weight_distributions_and_undetected_errors.py,
so a figure is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet: real
text in IBM Plex, colors from the page's tokens with paper fallbacks, and no
element that the site's figure check would refuse.
"""

import argparse
import os
from math import comb, log10

from ch23_weight_distributions_and_undetected_errors import (EXTENDED, GRID, HAMMING, TRIPLE, code_from_checks,
                                                             columns_dual, curve, dist, from_dual, peak, pud,
                                                             shortened)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
MINUS = "−"


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


class LinAxis:
    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = lo, hi, px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (v - self.lo) / (self.hi - self.lo)


class LogAxis:
    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = log10(lo), log10(hi), px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (log10(v) - self.lo) / (self.hi - self.lo)


def text(x: float, y: float, s: str, anchor: str = "start", fill: str = "", mono: bool = False,
         size: int = 0) -> str:
    attrs = [f'x="{fmt(x)}"', f'y="{fmt(y)}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if fill:
        attrs.append(f'fill="{fill}"')
    if mono:
        attrs.append(f'font-family="{MONO}"')
    if size:
        attrs.append(f'font-size="{size}"')
    return f'<text {" ".join(attrs)}>{s}</text>'


def power(base: str, exponent: str) -> str:
    """A base with a raised exponent, as text content."""
    return f'{base}<tspan dy="-5" font-size="9">{exponent}</tspan>'


def sci(v: float, digits: int = 3) -> str:
    """A number as mantissa times a power of ten, with the exponent raised."""
    mantissa, exponent = f"{v:.{digits - 1}e}".split("e")
    return f"{mantissa} \u00d7 " + power("10", MINUS + str(-int(exponent)) if int(exponent) < 0 else str(int(exponent)))


def spoken(v: float, digits: int = 3) -> str:
    """A small number in words for alt text, such as 6.79 times 10 to the minus 6."""
    mantissa, exponent = f"{v:.{digits - 1}e}".split("e")
    return f"{mantissa} times 10 to the minus {-int(exponent)}"


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "", opacity: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, width: float = 1.5, rx: float = 0) -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"/>')


def frame(xlabel: str, ylabel: str, x1: float = X1, h_: int = H) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {fmt(x1)},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + x1) / 2, h_ - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def x_ticks(ax, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        out.append(text(px, Y1 + 20, lab, "middle", mono=True, size=11))
    return out


def y_ticks(ax, values: list, labels: list, x1: float = X1) -> list:
    out = []
    for v, lab in zip(values, labels):
        py = ax(v)
        out.append(line(X0 - 4, py, X0, py))
        out.append(line(X0, py, x1, py, TEXT3, 1, opacity="0.35"))
        out.append(text(X0 - 8, py + 4, lab, "end", mono=True, size=11))
    return out


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_line(x: float, y: float, color: str, label: str, dash: str = "") -> list:
    return [line(x, y - 4, x + 24, y - 4, color, 2, dash), text(x + 32, y, label, size=12)]


def legend_box(x: float, y: float, fill: str, stroke: str, label: str) -> list:
    return [rect(x + 6, y - 10, 12, 12, fill, stroke, 1.5), text(x + 32, y, label, size=12)]


def grid_p() -> list:
    return [j / GRID for j in range(GRID // 2 + 1)]


def linear_ticks_p(x) -> list:
    values = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    return x_ticks(x, values, ["0", "0.1", "0.2", "0.3", "0.4", "0.5"])


# The charts --------------------------------------------------------------------------

def caught_and_missed() -> str:
    """Error patterns of each weight in 7 bits, split into those the (7,4) code catches and those it misses."""
    A = dist(HAMMING)
    x = LinAxis(0.5, 7.5, X0, X1)
    y = LinAxis(0, 40, Y1, Y0)
    body = frame("weight of the error pattern, the number of flipped bits", "number of error patterns")
    body += y_ticks(y, list(range(0, 41, 10)), [str(v) for v in range(0, 41, 10)])
    body += x_ticks(x, list(range(1, 8)), [str(w) for w in range(1, 8)])
    for w in range(1, 8):
        total, missed = comb(7, w), A[w]
        caught = total - missed
        left = x(w) - 28
        if missed:
            body.append(rect(left, y(missed), 56, y(0) - y(missed), BG2, RED, 1.5))
        if caught:
            body.append(rect(left, y(total), 56, y(missed) - y(total), BG2, GREEN, 1.5))
        if caught and missed:
            body.append(text(x(w), y(total) - 8, str(total), "middle", mono=True, size=11))
            body.append(text(x(w), y(missed / 2) + 4, f"{missed} missed", "middle", RED, mono=True, size=10))
        elif missed:
            body.append(text(x(w), y(total) - 8, f"{missed} missed", "middle", RED, mono=True, size=11))
        else:
            body.append(text(x(w), y(total) - 8, str(total), "middle", mono=True, size=11))
    body += legend_box(400, 44, BG2, GREEN, "caught, syndrome not 000")
    body += legend_box(400, 66, BG2, RED, "missed, the pattern is a codeword")
    counts = ", ".join(f"{comb(7, w)} of weight {w} with {A[w]} missed" for w in range(1, 8))
    return svg("Stacked bar chart of the 127 nonzero error patterns of 7 bits, grouped by weight, each bar split into "
               f"the patterns the (7,4) Hamming code catches and those it misses. The counts are {counts}. Every "
               "pattern of weight 1 or 2 is caught, and the 15 missed patterns are exactly the nonzero codewords, "
               "7 of weight 3, 7 of weight 4 and 1 of weight 7.", body)


def undetected_probability_curves() -> str:
    """P_ud of the (7,4) and extended (8,4) codes against p on log axes, with the levels 2^-3 and 2^-4."""
    x = LogAxis(1e-3, 0.5, X0, X1)
    y = LogAxis(1e-12, 1.0, Y1, Y0)
    body = frame("bit error probability p, log scale", "probability of an undetected error")
    body += y_ticks(y, [1e-12, 1e-9, 1e-6, 1e-3, 1.0],
                    [power("10", MINUS + "12"), power("10", MINUS + "9"), power("10", MINUS + "6"),
                     power("10", MINUS + "3"), "1"])
    body += x_ticks(x, [1e-3, 1e-2, 1e-1, 0.5], ["0.001", "0.01", "0.1", "0.5"])
    ps = [10 ** (-3 + i * (log10(0.5) + 3) / 200) for i in range(201)]
    codes = [(dist(HAMMING), GREEN, 3, "(7,4) Hamming"), (dist(EXTENDED), CYAN, 4, "extended (8,4)")]
    for A, color, r, _ in codes:
        body.append(line(X0, y(2.0 ** -r), X1, y(2.0 ** -r), color, 1, "3 3"))
        body.append(polyline([(x(p), y(pud(A, p))) for p in ps], color, 2))
        for p in (0.01, 0.1):
            body.append(dot(x(p), y(pud(A, p)), color, 4))
    A7, A8 = dist(HAMMING), dist(EXTENDED)
    body.append(text(x(0.01) - 8, y(pud(A7, 0.01)) - 10, sci(pud(A7, 0.01)), "end", GREEN, mono=True, size=11))
    body.append(text(x(0.01) + 10, y(pud(A8, 0.01)) + 18, sci(pud(A8, 0.01)), fill=CYAN, mono=True, size=11))
    body.append(text(X0 + 8, y(2.0 ** -3) - 6, power("2", MINUS + "3") + " = 0.125", fill=GREEN, mono=True, size=11))
    body.append(text(X0 + 8, y(2.0 ** -4) + 16, power("2", MINUS + "4") + " = 0.0625", fill=CYAN, mono=True, size=11))
    body += legend_line(360, 240, GREEN, "(7,4) Hamming code")
    body += legend_line(360, 262, CYAN, "extended (8,4) code")
    body += legend_line(360, 284, "currentColor", "rule of thumb " + power("2", MINUS + "r"), "3 3")
    return svg("Probability of an undetected error against the bit error probability p, both on log scales, for the "
               "(7,4) Hamming code and the extended (8,4) Hamming code. Both curves rise steadily from p = 0.001 to "
               f"p = 0.5. At p = 0.01 they are {spoken(pud(A7, 0.01))} and {spoken(pud(A8, 0.01))}, at p = 0.1 they are "
               f"{pud(A7, 0.1):.5g} and {pud(A8, 0.1):.5g}, and at p = 1/2 they reach 15/128 and 15/256, just under "
               "the dashed levels 2 to the minus 3 and 2 to the minus 4.", body)


def code_that_is_not_good() -> str:
    """P_ud of the triplicated-parity (8,5) code against the (7,4) code, on a linear p axis, with the level 1/8."""
    x1 = 540
    x = LinAxis(0, 0.5, X0, x1)
    y = LinAxis(0, 0.14, Y1, Y0)
    body = frame("bit error probability p", "probability of an undetected error", x1)
    body += y_ticks(y, [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14],
                    ["0", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12", "0.14"], x1)
    body += linear_ticks_p(x)
    body.append(line(X0, y(0.125), x1, y(0.125), "currentColor", 1, "3 3"))
    ps = grid_p()
    triple, hamming = curve(dist(TRIPLE)), curve(dist(HAMMING))
    body.append(polyline([(x(p), y(v)) for p, v in zip(ps, hamming)], GREEN, 2))
    body.append(polyline([(x(p), y(v)) for p, v in zip(ps, triple)], ORANGE, 2))
    top, j = peak(dist(TRIPLE))
    body.append(dot(x(j / GRID), y(top), ORANGE, 4))
    body.append(line(x(j / GRID), y(top) - 6, x(j / GRID), y(0.1385) + 4, ORANGE, 1))
    body.append(text(x(j / GRID), y(0.1385), f"{top:.5f} at p = 0.{j:03d}", "middle", ORANGE, mono=True, size=11))
    body.append(text(X0 + 8, y(0.125) - 6, "1/8 = 0.125", mono=True, size=11))
    body.append(text(x1 + 8, y(triple[-1]) - 2, "31/256", fill=ORANGE, mono=True, size=11))
    body.append(text(x1 + 8, y(hamming[-1]) + 12, "15/128", fill=GREEN, mono=True, size=11))
    body += legend_line(330, 252, ORANGE, "triplicated parity (8,5)")
    body += legend_line(330, 274, GREEN, "(7,4) Hamming code")
    body += legend_line(330, 296, "currentColor", "the level 1/8", "3 3")
    return svg("Probability of an undetected error against the bit error probability p from 0 to 1/2 for two codes "
               "with three parity checks. The (7,4) Hamming code rises steadily to 15/128 at p = 1/2 and stays under "
               "the dashed level 1/8. The (8,5) code that sends five data bits and their parity three times rises "
               f"above 1/8, peaks at {top:.5f} near p = 0.{j:03d}, and falls back to 31/256 at p = 1/2.", body)


def shortened_hamming_lengths() -> str:
    """P_ud against p for the (7,4) code cut to n = 4, 5, 6 and 7 bits, with the level 1/8."""
    x1 = 516
    x = LinAxis(0, 0.5, X0, x1)
    y = LinAxis(0, 0.14, Y1, Y0)
    body = frame("bit error probability p", "probability of an undetected error", x1)
    body += y_ticks(y, [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14],
                    ["0", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12", "0.14"], x1)
    body += linear_ticks_p(x)
    body.append(line(X0, y(0.125), x1, y(0.125), "currentColor", 1, "3 3"))
    body.append(text(x1 + 8, y(0.125) - 2, "1/8", mono=True, size=11))
    styles = {4: (VIOLET, "3 3", "1/16"), 5: (ORANGE, "", "3/32"), 6: (CYAN, "3 3", "7/64"), 7: (GREEN, "", "15/128")}
    ps = grid_p()
    finals = []
    for n, (color, dash, frac) in styles.items():
        A = dist(code_from_checks(shortened(n)))
        values = curve(A)
        body.append(polyline([(x(p), y(v)) for p, v in zip(ps, values)], color, 2, dash))
        body.append(text(x1 + 8, y(values[-1]) + 4, f"n = {n}, {frac}", fill=color, mono=True, size=11))
        finals.append(f"{frac} for n = {n}")
    return svg("Probability of an undetected error against the bit error probability p from 0 to 1/2 for the "
               "(7,4) Hamming code and for the codes made by keeping only its first 4, 5 or 6 positions, all with "
               "the same three parity checks. Every curve rises steadily and ends at p = 1/2 below the dashed level "
               "1/8, at " + ", ".join(finals) + ".", body)


def five_check_family() -> str:
    """2^5 times the largest P_ud on the grid for the codes of length 17 to 31 with five parity checks."""
    x = LinAxis(16.5, 31.5, X0, X1)
    y = LinAxis(0.96, 1.12, Y1, Y0)
    body = frame("length n of the code, all with five parity checks", "largest value on the grid, times 32")
    body += y_ticks(y, [0.96, 1.0, 1.04, 1.08, 1.12], ["0.96", "1.00", "1.04", "1.08", "1.12"])
    body += x_ticks(x, list(range(17, 32, 2)), [str(n) for n in range(17, 32, 2)])
    body.append(line(X0, y(1.0), X1, y(1.0), "currentColor", 1, "3 3"))
    body.append(text(X1 - 4, y(1.0) + 18, "1, the level " + power("2", MINUS + "5"), "end", mono=True, size=11))
    marks = []
    for n in range(17, 32):
        A = from_dual(columns_dual(list(range(1, n + 1)), 5), 5)
        top, j = peak(A)
        scaled = 32 * top
        bad = top > 2.0 ** -5
        body.append(dot(x(n), y(scaled), ORANGE if bad else GREEN, 5 if bad else 4))
        if bad:
            body.append(text(x(n) + 10, y(scaled) + 4, f"n = {n}: {scaled:.4f}, peak at p = 0.{j:03d}", fill=ORANGE,
                             mono=True, size=11))
            marks.append(f"n = {n} to {scaled:.4f} with its peak at p = 0.{j:03d}")
    body.append(text(x(31), y(1.0) - 12, "(31,26)", "middle", GREEN, mono=True, size=11))
    body += legend_box(100, 290, ORANGE, ORANGE, "rises above the level, not good")
    body += legend_box(100, 312, GREEN, GREEN, "largest value at p = 1/2, good")
    return svg("Dot plot of 32 times the largest probability of an undetected error on the grid p = 0 to 1/2, for the "
               "codes of length 17 to 31 whose parity-check columns are the numbers 1 to n written in five bits. Two "
               "lengths rise above the dashed level 1, " + ", and ".join(marks) + ". From n = 19 to 31 the largest "
               "value sits at p = 1/2 and equals 1 minus 2 raised to the power 5 minus n, which is below 1 by less "
               "than 0.0001. The full (31,26) Hamming code is the last dot.", body)


CHARTS = {
    "fig-23-3-caught-and-missed-by-weight.svg": caught_and_missed,
    "fig-23-5-undetected-probability-curves.svg": undetected_probability_curves,
    "fig-23-6-a-code-that-is-not-good.svg": code_that_is_not_good,
    "fig-23-7-shortened-hamming-lengths.svg": shortened_hamming_lengths,
    "fig-23-8-five-check-family.svg": five_check_family,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="figures", help="folder to write the SVG files into")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
