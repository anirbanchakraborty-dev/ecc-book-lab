"""Write the two computed figures of chapter 18 of Error Correction Coding as SVG files.

    python3 python/ch18_charts.py --out figures

Every drawn value comes from the functions in
ch18_syndrome_decoding_and_the_standard_array.py, so a figure is corrected by
fixing those functions and running this script again, never by editing its
output. The SVG follows the book's figure style sheet: real text in IBM Plex,
colors from the page's tokens with paper fallbacks, and no element that the
site's figure check would refuse.
"""

import argparse
import math
import os

from ch18_syndrome_decoding_and_the_standard_array import (
    C7, leader_counts, likelihood, p_correct, standard_array)

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
    def __init__(self, lo_exp: int, hi_exp: int, px_lo: float, px_hi: float):
        self.lin = LinAxis(lo_exp, hi_exp, px_lo, px_hi)

    def __call__(self, v: float) -> float:
        return self.lin(math.log10(v))


def text(x: float, y: float, s: str, anchor: str = "start", fill: str = "", mono: bool = False,
         size: int = 0, bold: bool = False) -> str:
    attrs = [f'x="{fmt(x)}"', f'y="{fmt(y)}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if fill:
        attrs.append(f'fill="{fill}"')
    if mono:
        attrs.append(f'font-family="{MONO}"')
    if size:
        attrs.append(f'font-size="{size}"')
    if bold:
        attrs.append('font-weight="600"')
    return f'<text {" ".join(attrs)}>{s}</text>'


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "", opacity: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}/>')


def polyline(points: list, stroke: str, width: float = 2, dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def frame(x1: float, xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {fmt(x1)},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + x1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


LIKELIHOOD_P = [(0.1, GREEN, ""), (0.5, TEXT3, "3 3"), (0.8, ORANGE, "")]


def likelihood_chart() -> str:
    """The likelihood p^w (1-p)^(6-w) of an error pattern of weight w, for w = 0 to 6,
    at three crossover probabilities, on a log scale. The flat line at p = 1/2 carries
    no dots, so that no two dots ever touch."""
    x1 = 530
    x = LinAxis(0, 6, X0 + 24, x1 - 24)
    y = LogAxis(-6, 0, Y1, Y0)
    body = frame(x1, "weight w of the error pattern, which is the distance to the codeword",
                 "likelihood, log scale")
    for k in range(-6, 1):
        py = y(10.0 ** k)
        body.append(line(X0, py, x1, py, TEXT3, 1, opacity="0.35"))
        label = "1" if k == 0 else f'10<tspan dy="-5" font-size="8">\u2212{-k}</tspan>'
        body.append(text(X0 - 8, py + 4, label, "end", mono=True, size=11))
    for w in range(7):
        body.append(line(x(w), Y1, x(w), Y1 + 4))
        body.append(text(x(w), Y1 + 18, str(w), "middle", mono=True, size=11))
    for p, color, dash in LIKELIHOOD_P:
        pts = [(x(w), y(likelihood(p, 6, w))) for w in range(7)]
        body.append(polyline(pts, color, 2, dash))
        if not dash:
            body.extend(dot(px, py, color, 3.5) for px, py in pts)
        body.append(text(x1 + 10, pts[-1][1] + 4, f"p = {p}", fill=color, mono=True, size=12))
    return svg("Likelihood of an error pattern against its weight for three crossover probabilities", body)


def correct_decoding_chart() -> str:
    """The probability of correct decoding of the (7,4) Hamming code with its syndrome
    table, against the probability that four bits sent with no code all arrive intact."""
    _, leaders = standard_array(C7, 7)
    counts = leader_counts(leaders, 7)
    x = LinAxis(0, 0.5, X0, X1)
    y = LinAxis(0, 1, Y1, Y0)
    body = frame(X1, "crossover probability p", "probability of correct decoding")
    for k in range(6):
        v = k / 10
        body.append(line(x(v), Y1, x(v), Y1 + 4))
        body.append(text(x(v), Y1 + 18, f"{v:.1f}", "middle", mono=True, size=11))
    for k in range(6):
        v = k / 5
        py = y(v)
        body.append(line(X0, py, X1, py, TEXT3, 1, opacity="0.35"))
        body.append(text(X0 - 8, py + 4, f"{v:.1f}", "end", mono=True, size=11))
    ps = [i / 200 for i in range(101)]
    body.append(polyline([(x(p), y((1 - p) ** 4)) for p in ps], ORANGE, 2, "3 3"))
    body.append(polyline([(x(p), y(p_correct(counts, 7, p))) for p in ps], GREEN, 2))
    for p in (0.01, 0.05, 0.1):
        body.append(dot(x(p), y(p_correct(counts, 7, p)), GREEN, 4))
        body.append(dot(x(p), y((1 - p) ** 4), ORANGE, 4))
    body.append(text(x(0.1) + 10, y(p_correct(counts, 7, 0.1)) - 6, f"{p_correct(counts, 7, 0.1):.3f}",
                     fill=GREEN, mono=True, size=12))
    body.append(text(x(0.1) - 10, y(0.9 ** 4) + 18, f"{0.9 ** 4:.3f}", "end", ORANGE, mono=True, size=12))
    meet = p_correct(counts, 7, 0.5)
    body.append(dot(x(0.5), y(meet), "currentColor", 4))
    lx, ly = 392, 54
    body.append(line(lx, ly - 4, lx + 28, ly - 4, GREEN, 2))
    body.append(text(lx + 36, ly, "(7,4) code, syndrome table", size=12))
    body.append(line(lx, ly + 18, lx + 28, ly + 18, ORANGE, 2, "3 3"))
    body.append(text(lx + 36, ly + 22, "four bits sent bare", size=12))
    body.append(dot(lx + 14, ly + 40, "currentColor", 4))
    body.append(text(lx + 36, ly + 44, f"both {meet:.4f} at p = 0.5", size=12))
    return svg("Probability of correct decoding for the (7,4) Hamming code and for four bits sent bare", body)


CHARTS = {
    "fig-18-2-likelihood-by-weight.svg": likelihood_chart,
    "fig-18-8-correct-decoding-probability.svg": correct_decoding_chart,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="figures", help="folder to write the SVG files into")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, build in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(build())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
