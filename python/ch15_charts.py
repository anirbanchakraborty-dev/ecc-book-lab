"""Write the two computed figures of chapter 15 of Error Correction Coding as SVG files.

    python3 python/ch15_charts.py --out figures

Every drawn value comes from the functions in ch15_vector_spaces_over_finite_fields.py,
so a figure is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet: real
text in IBM Plex, colors from the page's tokens with paper fallbacks, and no
element that the site's figure check would refuse.
"""

import argparse
import math
import os
from math import comb

from ch15_vector_spaces_over_finite_fields import basis_fraction, gaussian

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
    """A base-10 axis whose ends are exponents, so 10^lo sits at px_lo."""

    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lin = LinAxis(lo, hi, px_lo, px_hi)

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


def polyline(points: list, stroke: str, width: float = 1.5) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"/>'


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def ring(x: float, y: float, stroke: str, r: float = 6.5) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="none" stroke="{stroke}" stroke-width="1.5"/>'


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def x_ticks(ax, values: list) -> list:
    out = []
    for v in values:
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        out.append(text(px, Y1 + 20, str(v), "middle", mono=True, size=11))
    return out


def y_ticks(ax, values: list, labels: list) -> list:
    out = []
    for v, label in zip(values, labels):
        py = ax(v)
        out.append(line(X0 - 4, py, X0, py))
        out.append(line(X0, py, X1, py, TEXT3, 1, opacity="0.35"))
        out.append(text(X0 - 8, py + 4, label, "end", mono=True, size=11))
    return out


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


# The charts --------------------------------------------------------------------------

def subspaces_of_gf2_7() -> str:
    """Subspaces of GF(2)^7 of each dimension k, beside the number of ways to pick k of 7 positions."""
    n = 7
    ks = list(range(n + 1))
    subspaces = [gaussian(n, k, 2) for k in ks]
    subsets = [comb(n, k) for k in ks]
    ax = LinAxis(-0.5, 7.5, X0, X1)
    ay = LogAxis(-0.5, 5, Y1, Y0)
    body = frame("dimension k", "count, on a logarithmic scale")
    body += x_ticks(ax, ks)
    decades = [1, 10, 100, 1000, 10000, 100000]
    body += y_ticks(ay, decades, [str(v) for v in decades])
    body.append(polyline([(ax(k), ay(v)) for k, v in zip(ks, subsets)], CYAN))
    body.append(polyline([(ax(k), ay(v)) for k, v in zip(ks, subspaces)], GREEN))
    for k, v in zip(ks, subsets):
        body.append(ring(ax(k), ay(v), CYAN))
        body.append(text(ax(k), ay(v) + 22, str(v), "middle", CYAN, mono=True, size=11))
    for k, v in zip(ks, subspaces):
        body.append(dot(ax(k), ay(v), GREEN))
        body.append(text(ax(k), ay(v) - 14, str(v), "middle", GREEN, mono=True, size=11))
    lx, ly = X0 + 150, Y1 - 40
    body.append(dot(lx + 12, ly - 4, GREEN))
    body.append(text(lx + 32, ly, "subspaces of dimension k", size=12))
    body.append(ring(lx + 12, ly + 16, CYAN))
    body.append(text(lx + 32, ly + 20, "sets of k of the 7 positions", size=12))
    return svg("Subspaces of GF(2)^7 of each dimension beside the sets of k of 7 positions", body)


def basis_fraction_chart() -> str:
    """The fraction of lists of n vectors of GF(q)^n that are bases, for q = 2, 3, 5."""
    ns = list(range(1, 11))
    ax = LinAxis(0.5, 10.5, X0, X1)
    ay = LinAxis(0, 1, Y1, Y0)
    body = frame("length n, with n vectors chosen at random", "fraction that form a basis")
    body += x_ticks(ax, ns)
    ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    body += y_ticks(ay, ticks, [f"{v:.1f}" for v in ticks])
    for q, color in ((5, ORANGE), (3, CYAN), (2, GREEN)):
        values = [basis_fraction(q, n) for n in ns]
        body.append(polyline([(ax(n), ay(v)) for n, v in zip(ns, values)], color))
        for n, v in zip(ns, values):
            body.append(dot(ax(n), ay(v), color))
        body.append(text(ax(6.5), ay(values[-1]) + 22, f"q = {q}, settling near {values[-1]:.3f}",
                         fill=color, size=12))
    return svg("The fraction of lists of n vectors of GF(q)^n that form a basis, for q = 2, 3 and 5", body)


CHARTS = {
    "fig-15-4-subspace-counts.svg": subspaces_of_gf2_7,
    "fig-15-5-basis-fractions.svg": basis_fraction_chart,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="figures", help="directory for the SVG files")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
