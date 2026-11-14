"""Write the three charts of chapter 5 of Error Correction Coding as SVG files.

    python3 python/ch05_charts.py --out figures

Every plotted value comes from the functions in ch05_source_coding_and_prefix_codes.py,
so a chart is corrected by fixing those functions and running this script
again, never by editing its output. The SVG follows the book's figure style
sheet: real text in IBM Plex, colors from the page's tokens with paper
fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import os
from math import log2

from ch05_source_coding_and_prefix_codes import (LETTERS6, SIX, average_length, binary_entropy, block_rate,
                                                 entropy, huffman, shannon_lengths)

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


def rect(x: float, y: float, w: float, h_: float, fill: str = "none", stroke: str = "currentColor",
         width: float = 1, rx: float = 2) -> str:
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h_)}" rx="{fmt(rx)}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def lin_x_ticks(ax: LinAxis, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        out.append(text(px, Y1 + 20, lab, "middle", mono=True, size=11))
    return out


def lin_y_ticks(ax: LinAxis, values: list, labels: list) -> list:
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


def legend_row(x: float, y: float, kind: str, color: str, label: str, dash: str = "") -> list:
    """One legend entry: a swatch (a filled box or a short line) and its label."""
    if kind == "box":
        swatch = [rect(x, y - 9, 14, 10, color, TEXT3 if color == BG2 else color, 1)]
    else:
        swatch = [line(x, y - 4, x + 24, y - 4, color, 2, dash)]
    return swatch + [text(x + 32, y, label, size=12)]


# The charts --------------------------------------------------------------------------

def code_lengths() -> str:
    """Per letter of the six-letter source: the ideal length log2(1/p), Shannon's length and Huffman's."""
    probs = [float(p) for p in SIX]
    ideal = [log2(1 / p) for p in probs]
    shannon = shannon_lengths(SIX)
    huff = huffman(SIX, merged_high=True)[0]
    y = LinAxis(0, 5.5, Y1, Y0)
    group = (X1 - X0) / len(probs)
    bar, gap = 18, 4
    body = frame("letter and its probability", "codeword length in bits")
    body += lin_y_ticks(y, [0, 1, 2, 3, 4, 5], ["0", "1", "2", "3", "4", "5"])
    for i, (letter, p) in enumerate(zip(LETTERS6, probs)):
        centre = X0 + group * (i + 0.5)
        left = centre - (3 * bar + 2 * gap) / 2
        for k, (val, fill, label) in enumerate(((ideal[i], BG2, f"{ideal[i]:.2f}"),
                                                (shannon[i], ORANGE, str(shannon[i])),
                                                (huff[i], GREEN, str(huff[i])))):
            bx = left + k * (bar + gap)
            body.append(rect(bx, y(val), bar, Y1 - y(val), fill, TEXT3 if fill == BG2 else fill, 1, 1))
            body.append(text(bx + bar / 2, y(val) - 5, label, "middle", mono=True, size=9 if k == 0 else 11))
        body.append(text(centre, Y1 + 18, letter, "middle", mono=True, size=12))
        body.append(text(centre, Y1 + 33, f"{p:.2f}", "middle", mono=True, size=11, fill=TEXT3))
    h6 = entropy(SIX)
    body += legend_row(X0 + 16, 44, "box", BG2, f"ideal length log2(1/p), average H = {h6:.4f}")
    body += legend_row(X0 + 16, 64, "box", ORANGE, f"Shannon code, average {float(average_length(SIX, shannon)):.2f}")
    body += legend_row(X0 + 16, 84, "box", GREEN, f"Huffman code, average {float(average_length(SIX, huff)):.2f}")
    return svg("Grouped bars for the six letters a to f, with probabilities 0.45, 0.20, 0.10, 0.10, 0.10 and 0.05. "
               "For each letter three bars show the ideal length log2 of 1 over p, the Shannon code length and the "
               "Huffman code length. Ideal lengths are 1.15, 2.32, 3.32, 3.32, 3.32 and 4.32. Shannon lengths are "
               "2, 3, 4, 4, 4 and 5. Huffman lengths are 1, 3, 3, 3, 4 and 4. Averages: entropy 2.1955, Shannon "
               "2.95, Huffman 2.25.", body)


def block_coding() -> str:
    """Huffman bits per letter on blocks of n letters of the 0.9 and 0.1 source, n = 1 to 8."""
    hb = binary_entropy(0.1)
    ns = list(range(1, 9))
    rates = [block_rate(0.1, n) for n in ns]
    x = LinAxis(0.5, 8.5, X0, X1)
    y = LinAxis(0.3, 1.5, Y1, Y0)
    body = frame("block length n, in letters", "bits per source letter")
    body += lin_y_ticks(y, [0.4, 0.6, 0.8, 1.0, 1.2, 1.4], ["0.4", "0.6", "0.8", "1.0", "1.2", "1.4"])
    body += lin_x_ticks(x, ns, [str(n) for n in ns])
    body.append(line(X0, y(hb), X1, y(hb), TEXT3, 1.5, "3 3"))
    curve = [(x(n / 10), y(hb + 10 / n)) for n in range(10, 81)]
    body.append(polyline(curve, ORANGE, 1.5, "3 3"))
    body.append(polyline([(x(n), y(r)) for n, r in zip(ns, rates)], GREEN, 2))
    for n, r in zip(ns, rates):
        body.append(dot(x(n), y(r), GREEN))
    body += legend_row(380, 50, "line", GREEN, "Huffman code, bits per letter")
    body += legend_row(380, 72, "line", ORANGE, "upper bound h(0.1) + 1/n", "3 3")
    body += legend_row(380, 94, "line", TEXT3, f"entropy h(0.1) = {hb:.4f}", "3 3")
    return svg("Bits per source letter against block length n from 1 to 8, for a binary source with probabilities "
               "0.9 and 0.1. The Huffman points are 1, 0.645, 0.5327, 0.4926, 0.4802, 0.4702, 0.4743 and 0.4758. "
               "A dashed upper bound h + 1/n falls from 1.469 at n = 1 to 0.594 at n = 8, and a dashed horizontal "
               "line marks the entropy 0.4690. The Huffman points stay between the two lines and rise slightly at "
               "n = 7.", body)


def rate_against_p() -> str:
    """Huffman bits per letter on single letters, pairs and triples of a binary source, against h(p)."""
    x = LinAxis(0, 0.5, X0, X1)
    y = LinAxis(0, 1.05, Y1, Y0)
    body = frame("probability p of the rarer letter", "bits per source letter")
    body += lin_y_ticks(y, [0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
    body += lin_x_ticks(x, [0, 0.1, 0.2, 0.3, 0.4, 0.5], ["0", "0.1", "0.2", "0.3", "0.4", "0.5"])
    ps = [i / 1000 for i in range(10, 501)]
    body.append(line(x(0.01), y(1), x(0.5), y(1), TEXT3, 1.5, "3 3"))
    body.append(polyline([(x(p), y(block_rate(p, 2))) for p in ps], CYAN, 1.5))
    body.append(polyline([(x(p), y(block_rate(p, 3))) for p in ps], GREEN, 1.5))
    body.append(polyline([(x(p), y(binary_entropy(p))) for p in ps], "currentColor", 2))
    for val, fill in ((block_rate(0.1, 2), CYAN), (block_rate(0.1, 3), GREEN), (binary_entropy(0.1), "currentColor")):
        body.append(dot(x(0.1), y(val), fill, 3.5))
    body += legend_row(392, 176, "line", TEXT3, "Huffman, single letters", "3 3")
    body += legend_row(392, 198, "line", CYAN, "Huffman, pairs")
    body += legend_row(392, 220, "line", GREEN, "Huffman, triples")
    body += legend_row(392, 242, "line", "currentColor", "entropy h(p)")
    return svg("Bits per source letter against the probability p of the rarer letter of a binary source, p from "
               "0.01 to 0.5. Coding single letters always costs 1 bit. Coding pairs costs from 0.515 at p = 0.01 "
               "up to 1, coding triples from 0.353 up to 1, and the entropy h(p) rises from 0.081 to 1. The gap "
               "to the entropy is widest for small p. Dots mark p = 0.1: entropy 0.469, triples 0.533, pairs "
               "0.645.", body)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="figures", help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    charts = (("fig-5-6-code-lengths.svg", code_lengths()),
              ("fig-5-7-block-coding.svg", block_coding()),
              ("fig-5-8-rate-against-p.svg", rate_against_p()))
    for name, content in charts:
        with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
        print(f"wrote {os.path.join(out, name)}")


if __name__ == "__main__":
    main()
