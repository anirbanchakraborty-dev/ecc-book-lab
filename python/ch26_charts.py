"""Write the computed figure of chapter 26 of Error Correction Coding as an SVG file.

    python3 python/ch26_charts.py --out figures

Every drawn value comes from the functions in ch26_reed_muller_codes.py, so
the figure is corrected by fixing those functions and running this script
again, never by editing its output. The SVG follows the book's figure style
sheet: real text in IBM Plex, colors from the page's tokens with paper
fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import math
import os

from ch26_reed_muller_codes import direct_additions, fht_additions

W, H = 640, 400
X0, X1, Y0, Y1 = 92, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
TEXT3 = "var(--text-3,#5a5448)"
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
         dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 3.5) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def polyline(points: list, stroke: str, dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return ('<polyline points="' + " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
            + f'" fill="none" stroke="{stroke}" stroke-width="1.5"{d}/>')


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def fht_operation_count() -> str:
    """Additions to correlate a word against all 2^m first-order patterns, fast and direct."""
    ms = range(2, 11)
    ax = LinAxis(2, 10, X0 + 16, X1 - 64)
    ay = LinAxis(0, 7, Y1, Y0)
    body = [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 22, "number of variables m, code length 2^m", "middle"),
        f'<text transform="translate(22 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">'
        'additions and subtractions</text>',
    ]
    for e in range(0, 8):
        py = ay(e)
        body.append(line(X0 - 4, py, X0, py))
        body.append(line(X0, py, X1, py, TEXT3, 0.5, "3 3"))
        body.append(text(X0 - 8, py + 4, f"{10 ** e:,}", "end", mono=True, size=10))
    for m in ms:
        px = ax(m)
        body.append(line(px, Y1, px, Y1 + 4))
        body.append(text(px, Y1 + 18, str(m), "middle", mono=True, size=11))
    direct = [(ax(m), ay(math.log10(direct_additions(m)))) for m in ms]
    fast = [(ax(m), ay(math.log10(fht_additions(m)))) for m in ms]
    body.append(polyline(direct, "currentColor", "3 3"))
    body.append(polyline(fast, GREEN))
    for x, y in direct:
        body.append(dot(x, y, "currentColor"))
    for x, y in fast:
        body.append(dot(x, y, GREEN))
    xe, yd = direct[-1]
    _, yf = fast[-1]
    # The end labels sit above their points, clear of the gridline each point lands near.
    body.append(text(xe + 6, yd - 9, f"{direct_additions(10):,}", mono=True, size=11))
    body.append(text(xe + 6, yf - 9, f"{fht_additions(10):,}", fill=GREEN, mono=True, size=11))
    # The legend sits in the empty band between the gridlines for 1 and 10, right of m = 6.
    lx, ly = ax(6) + 20, ay(0.5) - 10
    body.append(line(lx, ly, lx + 24, ly, "currentColor", 1.5, "3 3"))
    body.append(text(lx + 32, ly + 4, "direct correlation, 2^m (2^m - 1)", size=12))
    body.append(line(lx, ly + 22, lx + 24, ly + 22, GREEN, 1.5))
    body.append(text(lx + 32, ly + 26, "fast Hadamard transform, m 2^m", size=12))
    return svg("Additions needed to correlate a received word with every first-order pattern, "
               "direct against fast, for m from 2 to 10", body)


CHARTS = {
    "fig-26-5-fht-operation-count.svg": fht_operation_count,
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
