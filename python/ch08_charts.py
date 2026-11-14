"""Write the four computed figures of chapter 8 of Error Correction Coding as SVG files.

    python3 python/ch08_charts.py --out figures

Every plotted value comes from the functions in ch08_distance_detection_and_correction.py,
so a figure is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet: real
text in IBM Plex, colors from the page's tokens with paper fallbacks, and no
element that the site's figure check would refuse.
"""

import argparse
import os
from fractions import Fraction
from math import comb, log10

from ch08_distance_detection_and_correction import (ball, cube_edges, digits, distance_distribution, isbn10_valid,
                                                    isbn13_valid, luhn_valid, make_valid, parity_coverage,
                                                    random_code, scheme_rates, weight)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"


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


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def x_ticks(ax, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        out.append(text(px, Y1 + 20, lab, "middle", mono=True, size=11))
    return out


def y_ticks(ax, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        py = ax(v)
        out.append(line(X0 - 4, py, X0, py))
        out.append(line(X0, py, X1, py, TEXT3, 1, opacity="0.35"))
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


# The charts --------------------------------------------------------------------------

def random_code_distances() -> str:
    """The 120 pairwise distances of the chapter 7 random code against the binomial expectation."""
    code = random_code(1, 16, 12)
    counts = distance_distribution(code)
    pairs = sum(counts.values())
    expected = [pairs * comb(12, k) / 4096 for k in range(13)]
    x = LinAxis(-0.5, 12.5, X0, X1)
    y = LinAxis(0, 30, Y1, Y0)
    body = frame("distance between the two codewords of a pair", "number of pairs")
    body += y_ticks(y, list(range(0, 31, 5)), [str(v) for v in range(0, 31, 5)])
    body += x_ticks(x, list(range(13)), [str(k) for k in range(13)])
    for k, c in counts.items():
        color = ORANGE if k == min(counts) else GREEN
        body.append(rect(x(k) - 14, y(c), 28, y(0) - y(c), BG2, color, 1.5))
    body.append(polyline([(x(k), y(e)) for k, e in enumerate(expected)], "currentColor", 1, "3 3"))
    body += [dot(x(k), y(e), "currentColor", 3) for k, e in enumerate(expected)]
    d = min(counts)
    body.append(text(x(d), y(11.6), "two pairs", "middle", ORANGE, size=12))
    body.append(text(x(d), y(10.1), f"at distance {d}", "middle", ORANGE, size=12))
    body.append(line(x(d), y(9.4), x(d), y(counts[d]) - 4, ORANGE, 1))
    body += legend_box(410, 44, BG2, GREEN, "pairs in the code")
    body += legend_line(410, 66, "currentColor", "binomial expectation", "3 3")
    listing = ", ".join(f"{counts.get(k, 0)} at {k}" for k in range(2, 11))
    return svg("Bar chart of the distances between all 120 pairs of codewords in a random code of 16 codewords of "
               f"12 bits. The counts are {listing}. Dots joined by a dashed line show the binomial expectation, "
               f"120 times C(12, k) over 4096, which peaks at {expected[6]:.1f} pairs at distance 6. The bar at "
               f"distance {d} is marked, since two pairs of codewords differ in only {d} positions.", body)


def four_cube_layers() -> str:
    """The 4-cube in layers of weight, with the radius-1 balls of the repetition code {0000, 1111}."""
    n = 4
    cols = [70, 195, 320, 445, 570]
    layers = [[v for v in range(1 << n) if weight(v) == k] for k in range(n + 1)]
    pos = {}
    for k, layer in enumerate(layers):
        for i, v in enumerate(layer):
            pos[v] = (cols[k], 190 + (i - (len(layer) - 1) / 2) * 48)
    near0, near1 = ball(0, n, 1), ball((1 << n) - 1, n, 1)
    body = []
    for a, b in cube_edges(n):
        lo, hi = (a, b) if weight(a) < weight(b) else (b, a)
        (xa, ya), (xb, yb) = pos[lo], pos[hi]
        body.append(line(xa + 24, ya, xb - 24, yb, TEXT3, 1, opacity="0.6"))
    for v, (px, py) in pos.items():
        color = GREEN if v in near0 else CYAN if v in near1 else ORANGE
        body.append(rect(px - 24, py - 12, 48, 24, BG2, color, 1.5, 3))
        body.append(text(px, py + 4, format(v, "04b"), "middle", color, mono=True, size=12))
    for k, x in enumerate(cols):
        body.append(text(x, 26, f"weight {k}", "middle", TEXT3, size=12))
    spans = [(cols[0] - 24, cols[1] + 24, GREEN, "within distance 1 of 0000"),
             (cols[2] - 24, cols[2] + 24, ORANGE, "distance 2 from both"),
             (cols[3] - 24, cols[4] + 24, CYAN, "within distance 1 of 1111")]
    for x_lo, x_hi, color, label in spans:
        body.append(line(x_lo, 354, x_hi, 354, color, 2))
        body.append(text((x_lo + x_hi) / 2, 376, label, "middle", color, size=12))
    covered = len(near0 | near1)
    return svg("The 16 words of 4 bits arranged in five columns by weight, 1, 4, 6, 4 and 1 words, with a line "
               "joining every two words that differ in one bit. The radius-1 ball around 0000 holds 0000 and the "
               "four words of weight 1, and the ball around 1111 holds 1111 and the four words of weight 3. The six "
               f"words of weight 2 lie in neither ball, so the two balls cover {covered} of the 16 words.", body)


def parity_coverage_chart() -> str:
    """Fraction of corrupted words that a single parity bit flags, against p, for three word lengths."""
    assert parity_coverage(8, Fraction(1, 2)) == Fraction(128, 255)
    assert all(abs(parity_coverage(n, 0.5) - 0.5) < 1e-9 for n in (32, 64))
    x = LogAxis(1e-4, 0.5, X0, X1)
    y = LinAxis(0.4, 1.0, Y1, Y0)
    body = frame("bit error probability p, log scale", "fraction of corrupted words flagged")
    body += y_ticks(y, [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], ["0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])
    body += x_ticks(x, [1e-4, 1e-3, 1e-2, 1e-1, 0.5], ["0.0001", "0.001", "0.01", "0.1", "0.5"])
    body.append(line(x(0.01), Y0, x(0.01), Y1, TEXT3, 1, "3 3"))
    grid = [10 ** (-4 + i * (log10(0.5) + 4) / 200) for i in range(201)]
    styles = [(8, GREEN, ""), (32, CYAN, "3 3"), (64, ORANGE, "")]
    for n, color, dash in styles:
        body.append(polyline([(x(p), y(parity_coverage(n, p))) for p in grid], color, 2, dash))
        body.append(dot(x(0.01), y(parity_coverage(n, 0.01)), color, 4))
    for i, (n, color, dash) in enumerate(styles):
        body += legend_line(100, 240 + 22 * i, color, f"n = {n}, {parity_coverage(n, 0.01):.3f} at p = 0.01", dash)
    return svg("Fraction of corrupted words that a single parity bit flags, plotted against the bit error "
               "probability p from 0.0001 to 0.5 on a log scale, for words of n = 8, 32 and 64 bits. Each curve "
               f"starts near 1 and falls toward one half. At p = 0.01, marked by a dashed vertical line, the "
               f"fractions are {parity_coverage(8, 0.01):.3f}, {parity_coverage(32, 0.01):.3f} and "
               f"{parity_coverage(64, 0.01):.3f}. At p = 0.5 the fraction is 128/255, or {parity_coverage(8, 0.5):.4f}, "
               "for n = 8, and within a billionth of one half for the longer words.", body)


def check_digit_coverage() -> str:
    """Percentage of single substitutions, adjacent swaps and jump swaps caught by three check-digit schemes."""
    schemes = [("ISBN-10", isbn10_valid, digits("0-13-042672-5"), GREEN),
               ("ISBN-13", isbn13_valid, digits("978-0-13-042672-7"), CYAN),
               ("Luhn", luhn_valid, make_valid(luhn_valid, digits("3141592653"), {}), ORANGE)]
    rates = [(name, scheme_rates(valid, base), color) for name, valid, base, color in schemes]
    groups = [("substitution", "single substitution"), ("adjacent", "adjacent swap"), ("jump", "swap two apart")]
    y = LinAxis(0, 100, Y1, Y0)
    body = frame("kind of error, tested at every position and every pair of digits", "errors caught, percent")
    body += y_ticks(y, list(range(0, 101, 20)), [str(v) for v in range(0, 101, 20)])
    for g, (key, group_label) in enumerate(groups):
        center = X0 + (X1 - X0) * (2 * g + 1) / 6
        for s, (name, rate, color) in enumerate(rates):
            cx = center + (s - 1) * 54
            value = 100 * float(rate[key])
            if value > 0:
                body.append(rect(cx - 18, y(value), 36, y(0) - y(value), BG2, color, 1.5))
            body.append(text(cx, y(value) - 6, f"{value:.1f}".rstrip("0").rstrip("."), "middle", mono=True, size=11))
            body.append(text(cx, Y1 + 16, name, "middle", color, size=10))
        body.append(text(center, Y1 + 36, group_label, "middle", size=12))
    summary = []
    for key, group_label in groups:
        part = ", ".join(f"{100 * float(rate[key]):.1f}" for _, rate, _ in rates)
        article = "an" if group_label[0] in "aeiou" else "a"
        summary.append(f"For {article} {group_label} they catch {part} percent")
    return svg("Grouped bar chart of the percentage of errors caught by the ISBN-10, ISBN-13 and Luhn check digits, "
               "in that order, for three kinds of error. " + ". ".join(summary) + ".", body)


CHARTS = {
    "fig-8-2-random-code-distances.svg": random_code_distances,
    "fig-8-4-four-cube-layers.svg": four_cube_layers,
    "fig-8-6-parity-coverage.svg": parity_coverage_chart,
    "fig-8-7-check-digit-coverage.svg": check_digit_coverage,
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
