"""Write the three charts of chapter 6 of Error Correction Coding as SVG files.

    python3 python/ch06_charts.py --out figures

Every plotted value comes from the functions in ch06_compression_in_practice.py
and ch05_source_coding_and_prefix_codes.py, so a chart is corrected by fixing
those functions and running this script again, never by editing its output.
The SVG follows the book's figure style sheet: real text in IBM Plex, colors
from the page's tokens with paper fallbacks, and no element that the site's
figure check would refuse.
"""

import argparse
import os
from math import log10

from ch05_source_coding_and_prefix_codes import block_rate
from ch06_compression_in_practice import (arith_expected_rate, best_uniform, binary_entropy, binary_source, db,
                                          entropy_coded_uniform, lz78_bits)

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


def square(x: float, y: float, fill: str, half: float = 4) -> str:
    return (f'<rect x="{fmt(x - half)}" y="{fmt(y - half)}" width="{fmt(2 * half)}" height="{fmt(2 * half)}" '
            f'fill="{fill}" stroke="{fill}" stroke-width="1"/>')


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


def legend_mark(x: float, y: float, color: str, label: str, shape: str = "dot") -> list:
    mark = dot(x + 12, y - 4, color, 4) if shape == "dot" else square(x + 12, y - 4, color)
    return [mark, text(x + 32, y, label, size=12)]


def geometric_grid(lo: int, hi: int, count: int) -> list:
    """About count integers spread evenly on a log scale from lo to hi, without repeats."""
    out = []
    for i in range(count):
        n = round(lo * (hi / lo) ** (i / (count - 1)))
        if not out or n > out[-1]:
            out.append(n)
    return out


# The charts --------------------------------------------------------------------------

def arithmetic_rate() -> str:
    """Expected bits per letter of the arithmetic code on n letters of the 0.9 and 0.1 source."""
    hb = binary_entropy(0.1)
    x = LogAxis(1, 1000, X0, X1)
    y = LinAxis(0.4, 1.4, Y1, Y0)
    body = frame("message length n in letters (log scale)", "bits per source letter")
    body += y_ticks(y, [0.4, 0.6, 0.8, 1.0, 1.2, 1.4], ["0.4", "0.6", "0.8", "1.0", "1.2", "1.4"])
    body += x_ticks(x, [1, 10, 100, 1000], ["1", "10", "100", "1000"])
    body.append(line(X0, y(hb), X1, y(hb), TEXT3, 1.5, "3 3"))
    ns = geometric_grid(2, 1000, 140)
    bound = [(x(n), y(hb + 2 / n)) for n in ns if hb + 2 / n <= 1.4]
    body.append(polyline(bound, ORANGE, 1.5, "3 3"))
    body.append(polyline([(x(n), y(arith_expected_rate(0.1, n))) for n in ns], GREEN, 2))
    for n in range(1, 9):
        body.append(square(x(n), y(block_rate(0.1, n)), CYAN, 3.5))
    body += legend_line(330, 50, GREEN, "arithmetic code, expected")
    body += legend_line(330, 72, ORANGE, "upper bound h(0.1) + 2/n", "3 3")
    body += legend_mark(330, 94, CYAN, "Huffman code on blocks of n", "square")
    body += legend_line(330, 116, TEXT3, f"entropy h(0.1) = {hb:.4f}", "3 3")
    return svg("Expected bits per source letter against message length n from 2 to 1000 on a log scale, for a "
               "binary source with probabilities 0.9 and 0.1. The arithmetic code costs 1.3 bits per letter at "
               "n = 2, 0.6751 at n = 8, 0.4841 at n = 100 and 0.4705 at n = 1000, always below the dashed bound "
               "h + 2/n. Squares mark Huffman codes on blocks of 1 to 8 letters, falling from 1 to 0.4758. A dashed "
               "line marks the entropy 0.4690.", body)


LZ78_LENGTHS = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]


def lz78_rate() -> str:
    """Bits per letter of LZ78 on one sequence of the 0.9 and 0.1 source, read to growing lengths."""
    hb = binary_entropy(0.1)
    bits = binary_source(0.1, LZ78_LENGTHS[-1])
    rates = [lz78_bits(bits[:n]) for n in LZ78_LENGTHS]
    x = LogAxis(100, 1000000, X0, X1)
    y = LinAxis(0.4, 1.05, Y1, Y0)
    body = frame("letters read, n (log scale)", "bits per source letter")
    body += y_ticks(y, [0.4, 0.6, 0.8, 1.0], ["0.4", "0.6", "0.8", "1.0"])
    body += x_ticks(x, [100, 1000, 10000, 100000, 1000000], ["100", "1000", "10^4", "10^5", "10^6"])
    body.append(line(X0, y(hb), X1, y(hb), TEXT3, 1.5, "3 3"))
    body.append(line(X0, y(1), X1, y(1), ORANGE, 1.5, "3 3"))
    body.append(polyline([(x(n), y(r)) for n, r in zip(LZ78_LENGTHS, rates)], GREEN, 2))
    for n, r in zip(LZ78_LENGTHS, rates):
        body.append(dot(x(n), y(r), GREEN, 3.5))
    body += legend_line(330, 72, GREEN, "LZ78, one sequence")
    body += legend_line(330, 94, ORANGE, "Huffman on single letters", "3 3")
    body += legend_line(330, 116, TEXT3, f"entropy h(0.1) = {hb:.4f}", "3 3")
    return svg("Bits per source letter of LZ78 against the number of letters read, from 100 to one million on a "
               f"log scale, for one sequence from a binary source with probabilities 0.9 and 0.1. The rate falls "
               f"from {rates[0]:.4f} at 100 letters to {rates[6]:.4f} at ten thousand and {rates[-1]:.4f} at one "
               f"million, still well above the dashed entropy line at {hb:.4f}. A second dashed line at 1 marks "
               "Huffman coding of single letters.", body)


def rate_distortion() -> str:
    """Distortion against rate for a unit-variance Gaussian: the bound, fixed-rate and entropy-coded quantizers."""
    x = LinAxis(0, 5, X0, X1)
    y = LinAxis(-32, 0, Y1, Y0)
    body = frame("rate R in bits per sample", "mean squared error in dB")
    body += y_ticks(y, [0, -5, -10, -15, -20, -25, -30], ["0", "-5", "-10", "-15", "-20", "-25", "-30"])
    body += x_ticks(x, [0, 1, 2, 3, 4, 5], ["0", "1", "2", "3", "4", "5"])
    body.append(polyline([(x(0), y(0)), (x(5), y(db(2 ** -10)))], "currentColor", 2))
    curve = []
    step = 3.0
    while step > 0.1:
        r, m = entropy_coded_uniform(step)
        if r <= 5:
            curve.append((x(r), y(db(m))))
        step *= 0.94
    body.append(polyline(curve, CYAN, 1.5, "3 3"))
    for rate in range(1, 6):
        body.append(dot(x(rate), y(db(best_uniform(2 ** rate)[1])), GREEN, 4))
    body += legend_line(330, 50, "currentColor", "bound D(R), 6.02 dB per bit")
    body += legend_mark(330, 72, GREEN, "best uniform quantizer, 2^R levels")
    body += legend_line(330, 94, CYAN, "uniform quantizer, entropy coded", "3 3")
    pts = ", ".join(f"{db(best_uniform(2 ** r)[1]):.2f}" for r in range(1, 6))
    return svg("Mean squared error in decibels against rate in bits per sample for a Gaussian source of unit "
               "variance. The rate-distortion bound is a straight line falling 6.02 dB per bit to -30.1 dB at 5 "
               f"bits. The best uniform quantizers with 2, 4, 8, 16 and 32 levels reach {pts} dB, drifting away "
               "from the bound. A dashed curve for entropy-coded uniform quantizers runs about 1.5 dB above the "
               "bound at high rates.", body)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="figures", help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    charts = (("fig-6-2-arithmetic-rate.svg", arithmetic_rate()),
              ("fig-6-5-lz78-rate.svg", lz78_rate()),
              ("fig-6-7-rate-distortion.svg", rate_distortion()))
    for name, content in charts:
        with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
        print(f"wrote {os.path.join(out, name)}")


if __name__ == "__main__":
    main()
